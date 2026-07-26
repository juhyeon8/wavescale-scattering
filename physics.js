// physics.js — RGD(Rayleigh–Gans–Debye) 소입자 산란 물리 코어
//
// 순수 계산만 담는다: DOM API도, Node 전용 API(console/fs/process 등)도 쓰지 않는다.
// Node에서는 require로, 브라우저에서는 <script src>로 그대로 재사용한다.
//
// 규약: 입사파 +z 진행, 전기장 x축 편광, 시간 규약 e^{-iωt}.
// 무차원화: 반지름 a=1 고정, λ가 조작 변인(x = ka = 2π/λ, a=1이므로 λ 값이 곧 λ/a).
// 총 분극률 α_total=1 고정, 개별 쌍극자 α = α_total/N.
//   (개별 α를 상수로 고정하면 N이 λ에 따라 바뀌면서 총 분극률이 격자 해상도에
//   오염되므로 반드시 α_total을 고정하고 N으로 나눈다.)
//
// α를 x(=ka)에 무관한 상수로 두는 물리적 근거: ω ≪ ω₀ 영역에서는 복원력이
// 지배하여 변위 진폭이 거의 일정하고, 따라서 분극률도 거의 일정하다.
// 이 시뮬레이션은 분광학적 축(ω/ω₀)이 아니라 기하학적 축(2πa/λ)만을 다룬다.
//
// 모형의 한계:
//   - Born 근사: 쌍극자 간 상호작용을 무시한다.
//   - RGD 적용 조건: |m−1| ≪ 1 이고 2x|m−1| ≪ 1.
//   - 그림자 효과와 흡수가 없으므로 x ≳ 3 에서는 물리적으로 신뢰할 수 없다.
//     기본 탐색 범위는 x ∈ [0.0063, 2] 로 제한한다.
//   - 목적은 정량적 예측이 아니라 경향과 기제의 확인이다.

var RGD = (function () {
  'use strict';

  // ---------------------------------------------------------------------
  // 이산화
  // ---------------------------------------------------------------------

  // 격자 간격 h가 아니라 지름당 격자점 수 n_d(정수)로 이산화를 지정한다.
  // 제약: k·h ≤ 0.3 ⟺ n_d ≥ 20x/3. 기하학적 충실도를 위해 n_d ≥ 20도 강제한다.
  function chooseNd(x) {
    return Math.max(20, Math.ceil((20 * x) / 3));
  }

  // 반지름 a의 구 안에 들어오는 정육면체 격자점만 취한다(난수 배치 금지 —
  // 무작위 배치는 |F|²에 샷 노이즈를 남긴다).
  //
  // 격자점은 셀 중심(cell-centered)으로 배치한다:
  //   x_i = −a + (i + 0.5)·h,  i = 0 … n_d−1  (y, z 동일)
  // 이렇게 하면 n_d의 홀짝과 무관하게 격자가 항상 원점 대칭이 되고,
  // Im(F)/N = 0 (검증 6)의 전제가 구현자의 임의 선택에 좌우되지 않는다.
  //
  // columns: [{x, z, m}] — 같은 (x_i, z_k) 열에 속한 y 격자점 개수 m.
  // φ=0 슬라이스에서는 q_y = 0 이므로 이 열 구조로 sigma1D를 가속할 수 있다.
  function buildDipoleGrid(a, nd) {
    var h = (2 * a) / nd;
    var coords = new Array(nd);
    for (var i = 0; i < nd; i++) {
      coords[i] = -a + (i + 0.5) * h;
    }

    var points = [];
    var columns = [];
    var a2 = a * a;

    for (var ix = 0; ix < nd; ix++) {
      var x = coords[ix];
      for (var iz = 0; iz < nd; iz++) {
        var z = coords[iz];
        var m = 0;
        for (var iy = 0; iy < nd; iy++) {
          var y = coords[iy];
          if (x * x + y * y + z * z <= a2) {
            points.push({ x: x, y: y, z: z });
            m++;
          }
        }
        if (m > 0) {
          columns.push({ x: x, z: z, m: m });
        }
      }
    }

    return { points: points, N: points.length, h: h, columns: columns };
  }

  // ---------------------------------------------------------------------
  // 산란 기하
  // ---------------------------------------------------------------------

  // q = k(ẑ − n̂),  n̂ = (sinθcosφ, sinθsinφ, cosθ)
  function qVector(k, theta, phi) {
    var st = Math.sin(theta);
    var ct = Math.cos(theta);
    var nx = st * Math.cos(phi);
    var ny = st * Math.sin(phi);
    var nz = ct;
    return { x: k * -nx, y: k * -ny, z: k * (1 - nz) };
  }

  // F(q) = Σ_j exp(i q·r_j)
  function formFactor(points, qVec) {
    var re = 0;
    var im = 0;
    for (var j = 0; j < points.length; j++) {
      var p = points[j];
      var phase = qVec.x * p.x + qVec.y * p.y + qVec.z * p.z;
      re += Math.cos(phase);
      im += Math.sin(phase);
    }
    return { re: re, im: im };
  }

  // φ=0 슬라이스 전용 최적화: 이 슬라이스에서는 항상 q_y=0 이므로
  // q·r_j에 y_j가 들어가지 않는다(근사가 아니라 정확한 항등식). 따라서
  // 격자점을 (x_i, z_k) 열로 묶고 그 열의 y 격자점 개수 m을 무게로 쓰면
  // formFactor(points, q)와 정확히 같은 값을 훨씬 적은 항으로 계산할 수
  // 있다 (n_d=20에서 항 수가 약 4200 → 314(≈π·n_d²/4)로 13배 감소).
  //   F = Σ_columns m(x_i, z_k) · exp(i(q_x·x_i + q_z·z_k))
  // sigma2D는 φ≠0을 다루므로 이 최적화를 쓸 수 없다.
  function formFactorPhi0Columns(columns, qx, qz) {
    var re = 0;
    var im = 0;
    for (var c = 0; c < columns.length; c++) {
      var col = columns[c];
      var phase = qx * col.x + qz * col.z;
      re += col.m * Math.cos(phase);
      im += col.m * Math.sin(phase);
    }
    return { re: re, im: im };
  }

  // 균일한 구의 해석적 형상 인자: f(u) = 3(sin u − u·cos u)/u³
  // u→0에서 f→1이지만 그 식은 0/0이므로, |u|가 작을 때는 급수 전개를 쓴다:
  //   f(u) = 1 − u²/10 + u⁴/280 − u⁶/15120 + ...
  function analyticShapeFactor(u) {
    if (Math.abs(u) < 0.1) {
      var u2 = u * u;
      return 1 - u2 / 10 + (u2 * u2) / 280 - (u2 * u2 * u2) / 15120;
    }
    return (3 * (Math.sin(u) - u * Math.cos(u))) / (u * u * u);
  }

  // ---------------------------------------------------------------------
  // Gauss–Legendre 구적 (Newton법, 순수 JS 구현, n별로 캐싱)
  // ---------------------------------------------------------------------

  var glCache = {};

  function gaussLegendreNodes(n) {
    if (glCache[n]) {
      return glCache[n];
    }

    var nodes = new Array(n);
    var weights = new Array(n);
    var m = Math.floor((n + 1) / 2);

    for (var i = 0; i < m; i++) {
      var z = Math.cos((Math.PI * (i + 0.75)) / (n + 0.5));
      var z1, pp, p1, p2, p3;
      do {
        p1 = 1;
        p2 = 0;
        for (var j = 0; j < n; j++) {
          p3 = p2;
          p2 = p1;
          p1 = ((2 * j + 1) * z * p2 - j * p3) / (j + 1);
        }
        pp = (n * (z * p1 - p2)) / (z * z - 1);
        z1 = z;
        z = z1 - p1 / pp;
      } while (Math.abs(z - z1) > 1e-15);

      nodes[i] = -z;
      nodes[n - 1 - i] = z;
      weights[i] = 2 / ((1 - z * z) * pp * pp);
      weights[n - 1 - i] = weights[i];
    }

    var result = { nodes: nodes, weights: weights };
    glCache[n] = result;
    return result;
  }

  // ---------------------------------------------------------------------
  // 산란 단면적
  // ---------------------------------------------------------------------

  // dσ/dΩ ∝ sin²Θ·|F(q)|²  — 편광 인자를 명시적으로 쓴다.
  // 극좌표 산란 패턴(편광면 φ=0, 수직면 φ=90°)의 기반 함수.
  //
  // 편광 인자 이중 계산 금지: 이 함수는 sin²Θ = 1 − (n̂·x̂)² 를 명시적으로
  // 쓴다. sigma1D의 (1+μ²) 인자는 이미 ∫sin²Θ dφ 를 적분해 둔 것이므로,
  // 이 함수와 sigma1D의 인자를 동시에 곱하면 안 된다.
  //
  // 예외: 이 함수는 각도마다 반복 호출되는 저수준 함수이므로, sigma1D/sigma2D와
  // 달리 k·h 자기 방어를 넣지 않는다(넣으면 경고가 폭주한다). 이는 의도된 설계다.
  function dSigmaDOmega(points, k, theta, phi) {
    var st = Math.sin(theta);
    var ndotx = st * Math.cos(phi);
    var polFactor = 1 - ndotx * ndotx;
    var q = qVector(k, theta, phi);
    var F = formFactor(points, q);
    var F2 = F.re * F.re + F.im * F.im;
    return polFactor * F2;
  }

  // φ 적분을 해석적으로 처리한 1차원(생산) 경로.
  //   σ = (8π/3)·k⁴·α_total²·G(x)
  //   G(x) = (3/8)·∫_{−1}^{1} (1+μ²)·|F(μ)/N|² dμ,   μ = cosθ
  // F(μ)는 실제 격자 합을 φ=0 슬라이스로 고정해 계산한다(등방성 가정).
  //
  // 편광 인자 이중 계산 금지: (1+μ²)를 쓴다 — 이미 ∫sin²Θ dφ 가 적분된
  // 결과이므로 별도로 sin²Θ를 곱하면 안 된다.
  //
  // grid는 buildDipoleGrid가 반환한 {points, N, h, columns} 객체를 받는다
  // (원시 points 배열이 아님) — h를 별도로 다시 추정하면 그 계산 자체가
  // 어긋날 수 있으므로, buildDipoleGrid가 반환한 h를 그대로 실어 나른다.
  //
  // 자기 방어: console을 쓸 수 없으므로(physics.js는 순수 계산만) 직접
  // 경고를 출력하는 대신 kh/khOk를 반환값에 실어, 호출자(verify.js)가
  // 판단해 출력하게 한다.
  function sigma1D(grid, k, alphaTotal, opts) {
    opts = opts || {};
    var n = opts.n || 32;
    var gl = gaussLegendreNodes(n);
    var N = grid.N;
    var sum = 0;

    for (var i = 0; i < n; i++) {
      var mu = gl.nodes[i];
      var w = gl.weights[i];
      var theta = Math.acos(mu);
      var q = qVector(k, theta, 0);
      // φ=0 슬라이스이므로 q.y는 항상 0 — 열 축약 최적화를 쓴다(근사 아님).
      var F = formFactorPhi0Columns(grid.columns, q.x, q.z);
      var Fmag2 = (F.re * F.re + F.im * F.im) / (N * N);
      sum += w * (1 + mu * mu) * Fmag2;
    }

    var G = (3 / 8) * sum;
    var sigma = ((8 * Math.PI) / 3) * Math.pow(k, 4) * alphaTotal * alphaTotal * G;
    var kh = k * grid.h;

    return { sigma: sigma, G: G, kh: kh, khOk: kh <= 0.3 + 1e-9 };
  }

  // 검증 전용 완전한 2차원 적분 (θ: Gauss–Legendre, φ: 균일 격자).
  // 1D 경로의 등방성 가정(φ=0 슬라이스로 F를 대표시키는 것)이 실제로
  // 성립하는지 검증하는 목적이므로, sin²Θ|F/N|²를 φ에 대해서도 그대로
  // 수치 적분한다(해석적 φ-적분을 쓰지 않는다).
  //
  // σ = k⁴·α_total² · ∫∫ sin²Θ·|F(q)/N|² dΩ
  //
  // 편광 인자 이중 계산 금지: dSigmaDOmega가 이미 sin²Θ = 1−(n̂·x̂)² 를
  // 명시적으로 곱하므로, 여기서 다시 곱하지 않는다.
  //
  // sigma1D와 시그니처를 통일한다: grid 객체를 받고, 동일한 k·h≤0.3 자기
  // 방어를 갖는다.
  function sigma2D(grid, k, alphaTotal, opts) {
    opts = opts || {};
    var nTheta = opts.nTheta || 32;
    var nPhi = opts.nPhi || 16;
    var gl = gaussLegendreNodes(nTheta);
    var N = grid.N;
    var integral = 0;

    for (var i = 0; i < nTheta; i++) {
      var mu = gl.nodes[i];
      var w = gl.weights[i];
      var theta = Math.acos(mu);
      var phiSum = 0;
      for (var j = 0; j < nPhi; j++) {
        var phi = (2 * Math.PI * j) / nPhi;
        phiSum += dSigmaDOmega(grid.points, k, theta, phi);
      }
      var phiIntegral = ((2 * Math.PI) / nPhi) * phiSum;
      integral += w * phiIntegral;
    }

    integral /= N * N;

    var sigma = Math.pow(k, 4) * alphaTotal * alphaTotal * integral;
    var kh = k * grid.h;

    return { sigma: sigma, kh: kh, khOk: kh <= 0.3 + 1e-9 };
  }

  // ---------------------------------------------------------------------
  // 근접장 (완전한 쌍극자 장 — 원거리 근사를 쓰지 않는다)
  // ---------------------------------------------------------------------

  // Jackson 9.18. 단위 x-편극 쌍극자(p̂ = x̂) 기준 커널.
  // 항등식 (R̂×p̂)×R̂ = p̂ − R̂(R̂·p̂) 로 전개하면:
  //
  //   E(R) = e^{ikR}/R · { k²[p̂ − R̂(R̂·p̂)] + [3R̂(R̂·p̂) − p̂]·(1/R² − ik/R) }
  //
  // 이것을 2갈래(복사/준정전기)가 아니라 반드시 3갈래로 분리한다.
  // (1/R² − ik/R) 묶음 안에 k 의존성이 서로 다른 두 항이 섞여 있기 때문이다 —
  // 묶음의 진폭은 (1/R³)·√(1+(kR)²) 이므로 "준정전기항은 λ 무관"은 틀렸다.
  //
  //   복사 radiation : k²[p̂ − R̂(R̂·p̂)]        · e^{ikR}/R    → 1/R,  k²  → λ⁻²
  //   유도 induction : [3R̂(R̂·p̂) − p̂]·(−ik/R) · e^{ikR}/R    → 1/R², k¹  → λ⁻¹
  //   정전 static    : [3R̂(R̂·p̂) − p̂]·(1/R²)  · e^{ikR}/R    → 1/R³, k⁰  → λ 무관
  //
  // 거리 차수와 k 차수가 1:1로 맞물리므로("1/R³↔k⁰, 1/R²↔k¹, 1/R↔k²"),
  // "파장이 길어질수록 먼 곳까지 가는 성분부터 순서대로 사라진다"가 그대로 읽힌다.
  //
  // 반환: { radiation, induction, static } — 각각 {ex, ey, ez}, 각 성분 {re, im}.
  // 세 항을 합쳐 전체 장을 만드는 것은 호출자의 책임이다(항별 표시가 UI 요구사항).
  function dipoleFieldKernel(dx, dy, dz, k) {
    var R2 = dx * dx + dy * dy + dz * dz;
    var R = Math.sqrt(R2);
    var R3 = R2 * R;

    // R̂ 성분과 c = R̂·p̂ = R̂_x
    var ux = dx / R;
    var uy = dy / R;
    var uz = dz / R;
    var c = ux;

    // V1 = p̂ − R̂(R̂·p̂)   (복사항의 벡터 인자)
    var v1x = 1 - ux * c;
    var v1y = -uy * c;
    var v1z = -uz * c;

    // V2 = 3R̂(R̂·p̂) − p̂   (유도·정전항이 공유하는 벡터 인자)
    var v2x = 3 * ux * c - 1;
    var v2y = 3 * uy * c;
    var v2z = 3 * uz * c;

    var kR = k * R;
    var cosKR = Math.cos(kR);
    var sinKR = Math.sin(kR);

    // 복사: 계수 = k²/R · e^{ikR}
    var radA = (k * k) / R;
    var radRe = radA * cosKR;
    var radIm = radA * sinKR;

    // 유도: 계수 = (−ik/R)·e^{ikR}/R = (k/R²)·(sin kR − i·cos kR)
    var indA = k / R2;
    var indRe = indA * sinKR;
    var indIm = -indA * cosKR;

    // 정전: 계수 = (1/R²)·e^{ikR}/R = (1/R³)·e^{ikR}
    var staRe = cosKR / R3;
    var staIm = sinKR / R3;

    return {
      radiation: {
        ex: { re: v1x * radRe, im: v1x * radIm },
        ey: { re: v1y * radRe, im: v1y * radIm },
        ez: { re: v1z * radRe, im: v1z * radIm }
      },
      induction: {
        ex: { re: v2x * indRe, im: v2x * indIm },
        ey: { re: v2y * indRe, im: v2y * indIm },
        ez: { re: v2z * indRe, im: v2z * indIm }
      },
      static: {
        ex: { re: v2x * staRe, im: v2x * staIm },
        ey: { re: v2y * staRe, im: v2y * staIm },
        ez: { re: v2z * staRe, im: v2z * staIm }
      }
    };
  }

  // 전체 격자에 대한 산란 근접장 — 브루트포스 합(대칭 최적화 없음).
  // verify.js 항목 7·8 전용이며, 임의의 3D 관측점을 지원한다.
  // UI 최적화 경로(scatteredFieldXZPlane)와 분리해 두어야 검증이 커널 자체를
  // 검증하는 의미를 갖는다.
  //
  // 개별 쌍극자: p_j = α·E₀·x̂·exp(i k z_j),  **α = alphaTotal/N** (개별 분극률).
  // alphaTotal을 그대로 곱하면 N²배(n_d=20에서 약 1.8e7배) 어긋난다.
  //
  // 반환: { radiation, induction, static } — 각각 {ex, ey, ez}, 각 성분 {re, im}.
  // 항목 8이 항별로 검증하므로 합치지 않고 분리해 반환한다.
  function scatteredFieldAt(obsPoint, grid, k, alphaTotal) {
    var points = grid.points;
    var N = grid.N;
    var alpha = alphaTotal / N;

    var acc = {
      radiation: { ex: { re: 0, im: 0 }, ey: { re: 0, im: 0 }, ez: { re: 0, im: 0 } },
      induction: { ex: { re: 0, im: 0 }, ey: { re: 0, im: 0 }, ez: { re: 0, im: 0 } },
      static: { ex: { re: 0, im: 0 }, ey: { re: 0, im: 0 }, ez: { re: 0, im: 0 } }
    };
    var termNames = ['radiation', 'induction', 'static'];
    var compNames = ['ex', 'ey', 'ez'];

    for (var j = 0; j < points.length; j++) {
      var p = points[j];
      var kern = dipoleFieldKernel(obsPoint.x - p.x, obsPoint.y - p.y, obsPoint.z - p.z, k);

      // 구동 진폭 α·exp(i k z_j)
      var da = alpha * Math.cos(k * p.z);
      var db = alpha * Math.sin(k * p.z);

      for (var t = 0; t < 3; t++) {
        var tn = termNames[t];
        for (var m = 0; m < 3; m++) {
          var cn = compNames[m];
          var v = kern[tn][cn];
          // 복소수 곱: (re + i·im)·(da + i·db)
          acc[tn][cn].re += v.re * da - v.im * db;
          acc[tn][cn].im += v.re * db + v.im * da;
        }
      }
    }

    return acc;
  }

  // UI 전용 최적화: 관측점이 xz 평면(y_obs = 0) 위에 있을 때의 산란 근접장.
  //
  // [y 대칭 — 이 함수 안에서 쓰는 유일한 대칭]
  //   dy = −y_j 이고, 커널의 E_x·E_z 성분은 dy에 대해 짝함수(V1x, V1z, V2x, V2z가
  //   uy를 짝수 번만 포함), E_y 성분은 홀함수다. 격자가 y_j ↔ −y_j로 짝지어져
  //   있고 구동 위상 exp(i k z_j)는 z에만 의존하므로 짝의 진폭이 같다.
  //   따라서 y_j > 0인 쌍극자만 돌고 E_x, E_z에 2배 하면 정확하다(근사 아님).
  //   E_y는 상쇄되어 항상 0이므로 계산도 반환도 하지 않는다.
  //   셀 중심 배치라 y=0 격자점이 존재하지 않으므로 이중 계산 위험이 없다.
  //
  // [x 대칭은 여기에 넣지 않는다 — 넣으면 틀린다]
  //   관측점 x_o가 고정된 상태에서 dx = x_o − x_j 와 x_o − (−x_j) = x_o + x_j 는
  //   부호만 다른 관계가 아니다(x_o = 0에서만 성립). 실제로 성립하는 것은
  //   전체 합에 대한 관계이며,
  //     E_x(−x_o, z_o) = +E_x(x_o, z_o),  E_z(−x_o, z_o) = −E_z(x_o, z_o)
  //   커널이 dx에 대해 E_x는 짝함수·E_z는 홀함수라는 성질과, 격자가 x→−x
  //   대칭이라는 조건이 함께 있어야 성립한다. 단일 관측점을 받는 이 함수는
  //   관측점 쌍을 만들 수 없으므로, x 대칭은 호출자(렌더 루프)가 적용한다.
  //   (비구형 산란체로 확장하면 후자가 깨지므로 이 최적화도 무효가 된다.)
  //
  // α = alphaTotal/N (개별 분극률) — scatteredFieldAt과 동일.
  //
  // 반환: { radiation:{ex,ez}, induction:{ex,ez}, static:{ex,ez} } (각 {re,im}).
  //   세 항을 항상 함께 반환한다 — 화면 점마다 세 항을 캐싱해야 (B)의
  //   [전체]/[복사]/[유도]/[정전] 토글이 재계산 없이 표시만 바꾸는 일이 된다.
  //   opts.term('full'|'radiation'|'induction'|'static')을 주면 그 항의 합을
  //   최상위 {ex, ez}로도 함께 실어 준다(단일 지점 조회 편의용).
  function scatteredFieldXZPlane(xObs, zObs, grid, k, alphaTotal, opts) {
    opts = opts || {};
    var points = grid.points;
    var N = grid.N;
    var alpha = alphaTotal / N;

    var radExRe = 0, radExIm = 0, radEzRe = 0, radEzIm = 0;
    var indExRe = 0, indExIm = 0, indEzRe = 0, indEzIm = 0;
    var staExRe = 0, staExIm = 0, staEzRe = 0, staEzIm = 0;

    for (var j = 0; j < points.length; j++) {
      var p = points[j];
      if (p.y <= 0) continue; // y>0 절반만 순회하고 마지막에 2배

      var dx = xObs - p.x;
      var dy = -p.y;
      var dz = zObs - p.z;

      var kern = dipoleFieldKernel(dx, dy, dz, k);

      var da = alpha * Math.cos(k * p.z);
      var db = alpha * Math.sin(k * p.z);

      var v;
      v = kern.radiation.ex;
      radExRe += v.re * da - v.im * db;
      radExIm += v.re * db + v.im * da;
      v = kern.radiation.ez;
      radEzRe += v.re * da - v.im * db;
      radEzIm += v.re * db + v.im * da;

      v = kern.induction.ex;
      indExRe += v.re * da - v.im * db;
      indExIm += v.re * db + v.im * da;
      v = kern.induction.ez;
      indEzRe += v.re * da - v.im * db;
      indEzIm += v.re * db + v.im * da;

      v = kern.static.ex;
      staExRe += v.re * da - v.im * db;
      staExIm += v.re * db + v.im * da;
      v = kern.static.ez;
      staEzRe += v.re * da - v.im * db;
      staEzIm += v.re * db + v.im * da;
    }

    var result = {
      radiation: {
        ex: { re: 2 * radExRe, im: 2 * radExIm },
        ez: { re: 2 * radEzRe, im: 2 * radEzIm }
      },
      induction: {
        ex: { re: 2 * indExRe, im: 2 * indExIm },
        ez: { re: 2 * indEzRe, im: 2 * indEzIm }
      },
      static: {
        ex: { re: 2 * staExRe, im: 2 * staExIm },
        ez: { re: 2 * staEzRe, im: 2 * staEzIm }
      }
    };

    if (opts.term) {
      if (opts.term === 'full') {
        result.ex = {
          re: result.radiation.ex.re + result.induction.ex.re + result.static.ex.re,
          im: result.radiation.ex.im + result.induction.ex.im + result.static.ex.im
        };
        result.ez = {
          re: result.radiation.ez.re + result.induction.ez.re + result.static.ez.re,
          im: result.radiation.ez.im + result.induction.ez.im + result.static.ez.im
        };
      } else {
        result.ex = result[opts.term].ex;
        result.ez = result[opts.term].ez;
      }
    }

    return result;
  }

  return {
    chooseNd: chooseNd,
    buildDipoleGrid: buildDipoleGrid,
    qVector: qVector,
    formFactor: formFactor,
    formFactorPhi0Columns: formFactorPhi0Columns,
    analyticShapeFactor: analyticShapeFactor,
    gaussLegendreNodes: gaussLegendreNodes,
    dSigmaDOmega: dSigmaDOmega,
    sigma1D: sigma1D,
    sigma2D: sigma2D,
    dipoleFieldKernel: dipoleFieldKernel,
    scatteredFieldAt: scatteredFieldAt,
    scatteredFieldXZPlane: scatteredFieldXZPlane
  };
})();

if (typeof module !== 'undefined' && module.exports) module.exports = RGD;
