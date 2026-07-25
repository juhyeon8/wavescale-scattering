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
    sigma2D: sigma2D
  };
})();

if (typeof module !== 'undefined' && module.exports) module.exports = RGD;
