// script.js — RGD 소입자 산란 시각화 (탐색 화면)
//
// 1행 (A)(B)(C): xz 평면 단면의 공간 장 분포. physics.js의 근접장 API
// (dipoleFieldKernel / scatteredFieldXZPlane)를 쓴다.
//
// 계산 격자와 표시 격자는 별개다:
//   - σ/G 계산 경로(sigma1D/sigma2D, nd=20)는 절대 건드리지 않는다.
//   - 장 지도 렌더링용으로만 별도의 coarse 격자(nd_display=14)를 쓴다.

(function () {
  'use strict';

  // ---------------------------------------------------------------------
  // 설정
  // ---------------------------------------------------------------------

  // 표시 전용 격자. k·h = 0.286 이므로 x=2 에서도 k·h ≤ 0.3 제약을 만족한다.
  // 2행 (D)(E)는 이 값이 아니라 물리 격자 nd=20 을 쓴다.
  var ND_DISPLAY = 14;

  // 표시 범위는 λ 기준이 기본이다: z 범위 = 3λ 고정.
  // 그러면 λ를 바꿔도 (A) 입사파 줄이 언제나 똑같이 생기고(마루 세 개),
  // 변하는 것이 (B) 하나뿐이라 원인과 결과가 분리되어 보인다. 동시에 입자는
  // λ에 비해 점점 작은 점이 되어 λ/a 가 화면에서 직접 읽힌다.
  var RANGE_LAMBDA_Z = 3;  // z 범위 = 3λ
  var RANGE_A_Z = 16;      // 'a 기준' 토글: z ∈ [−8a, 8a]

  // 화면 격자는 패널 가로세로비에 맞춰 자동 산출한다(픽셀 정사각형 유지).
  //   NZ * NX ≈ TOTAL_SAMPLES,  NZ / NX ≈ aspect
  var TOTAL_SAMPLES = 8000;
  var Z_HALF = 7.5; // a 단위. updateGridDims/λ 에 따라 다시 잡힌다
  var X_HALF = 1.1;
  var NZ = 216;
  var NX = 36;      // 짝수여야 x 거울 대칭(절반만 계산)이 성립한다

  // 마스크 반경은 1.5a로 고정한다 — 축소하지 않는다.
  // nd_display=14, 셀 중심 배치이므로 최외곽 쌍극자는 r ≈ 0.93~1.00a 에 있다.
  // 1.2a로 줄이면 최근접 쌍극자까지 0.21a 밖에 안 되어, 1/R³ 발산 때문에
  // 단일 쌍극자 하나가 전체 장의 최대 13%를 차지한다 → 격자점 위치가
  // 알갱이 무늬로 드러나 논문 그림으로 쓸 수 없다. 1.5a에서는 1~3%다.
  var MASK_R = 1.5;

  // --- 실제 단위 ---------------------------------------------------------
  // 길이의 기준 단위는 μm 다. a = 1 고정을 해제했으므로 a 와 λ 가 각각
  // 독립 슬라이더이고, 물리량(σ, Q)이 실제 단위로 나온다.
  //
  // 그림 자체는 a 에 의존하지 않는다: 격자를 a 배로 늘리면 k 는 1/a 배가 되고
  // α ∝ a³ 이므로 세 항(k²/R, k/R², 1/R³)에서 a 가 정확히 상쇄된다.
  // 즉 λ/a 가 같으면 그림이 같다 — "λ/a 고정" 체크박스가 그것을 보여준다.
  var A_MIN = 0.01, A_MAX = 10, A_DEFAULT = 0.1;         // μm
  var LAM_MIN = 0.1, LAM_MAX = 1e6, LAM_DEFAULT = 0.5;   // μm (100nm ~ 1m)

  // 규산염 m ≈ 1.7 의 (m²−1)/(m²+2). α_total = 0.39·a³ 이 된다.
  // 계수 없이 a³ 만 쓰면 완전도체 구의 상한이라 λ/a = 5 에서 Q_sca 가 3.7까지
  // 올라간다. UI 컨트롤로 노출하지 않는 상수다.
  var POLARIZABILITY_FACTOR = 0.39;

  // RGD 적용 하한. λ/a < π (x > 2) 이면 |m−1| ≪ 1 조건이 깨진다.
  // 계산은 계속하되 화면을 흐리게 하고 숫자를 회색으로 둔다.
  var LAMBDA_OVER_A_MIN = Math.PI;

  // 개발용 대칭 어서션. 기본값 false — true 로 두면 미러링으로 채운 열 하나를
  // scatteredFieldXZPlane 으로 직접 계산해 대조한다.
  //
  // 왜 필요한가: verify.js 항목 9b는 x 거울 "관계식"이 성립함을 확인할 뿐,
  // 이 렌더 루프가 그 관계를 올바로 적용했는지는 검증하지 못한다
  // (x 대칭은 physics.js가 아니라 여기에 있다).
  var DEBUG_SYMMETRY = false;

  var TERM_NAMES = ['radiation', 'induction', 'static'];
  var TERM_LABELS = { full: '전체', radiation: '복사', induction: '유도', static: '정전' };

  // ---------------------------------------------------------------------
  // 상태
  // ---------------------------------------------------------------------

  // 2행 (D)(E)는 표시 격자가 아니라 물리 격자 nd=20 을 쓴다 —
  // (E)가 nd=14를 쓰면 (F) 산란 패턴과 다른 격자의 결과를 나란히 놓게 된다.
  var ND_PHYS = 20;
  var THETA_SAMPLES = 181; // (F) 곡선: 1° 간격

  var aUm = A_DEFAULT;
  var lambdaUm = LAM_DEFAULT;
  var lockRatio = false;                        // λ/a 고정 체크박스
  var lockedRatio = LAM_DEFAULT / A_DEFAULT;    // 고정할 당시의 λ/a
  var atLimit = false;                          // 고정 상태에서 슬라이더 끝에 닿음

  // 격자는 a 에 의존하므로 a 가 바뀔 때마다 다시 만든다.
  var grid = null;      // 표시용 nd = 14
  var physGrid = null;  // 물리 nd = 20 (σ/G 경로와 동일)

  function rebuildGrids() {
    grid = RGD.buildDipoleGrid(aUm, ND_DISPLAY);
    physGrid = RGD.buildDipoleGrid(aUm, ND_PHYS);
  }
  rebuildGrids();

  function lambdaOverA() { return lambdaUm / aUm; }
  function waveK() { return (2 * Math.PI) / lambdaUm; }   // μm⁻¹
  function sizeX() { return (2 * Math.PI * aUm) / lambdaUm; } // x = ka
  function alphaTotal() { return POLARIZABILITY_FACTOR * aUm * aUm * aUm; }

  var theta = Math.PI / 2; // (E)(F) 전용 θ. 탭 2 진입 시 기본 90°
  var activeTab = 1;       // a·λ 는 두 탭이 공유한다(탭을 바꿔도 유지)
  var term = 'full';
  // 배율 기본값은 ×1 이다. 자동 배율이면 λ를 늘려도 (B)가 늘 비슷하게 보여서
  // "산란파가 약해진다"는 메시지가 통째로 사라진다.
  var scaleMode = '1';
  var rangeMode = 'lambda';
  var showArrows = true;
  // 지표는 기본 숨김. 켠 상태는 세션 동안 유지된다(λ를 바꿔도 꺼지지 않는다).
  var showMetrics = false;

  // --- 애니메이션 --------------------------------------------------------
  // 화면 주기를 λ와 무관하게 고정한다. 실제 ω를 쓰면 21cm 전파에서는 멈춰
  // 보이고 가시광선에서는 폭주한다.
  var SCREEN_PERIOD = 2.0; // 초
  var animPhase = 0;       // ωt
  var animating = true;
  var animSpeed = 1;
  var lastFrameT = 0;
  var rafId = null;

  // 화면 점마다 세 항의 복소 진폭을 캐싱한다 — 프레임마다 하는 일은
  // val = Re·cos(ωt) + Im·sin(ωt) 곱셈 두 번뿐이고, 무거운 쌍극자 합은
  // a 또는 λ가 바뀔 때만 돈다. (시간 규약 e^{-iωt})
  // full 은 세 항의 합을 미리 더해 둔 것이다(프레임마다 더하지 않기 위해).
  var FIELD_KEYS = ['radiation', 'induction', 'static', 'full'];
  var fieldRe = {}, fieldIm = {};
  var incidentRe = null, incidentIm = null, maskFlag = null;
  var metrics = { radiation: 0, induction: 0, static: 0, maxScaAmp: 0 };

  function allocBuffers() {
    var n = NZ * NX;
    for (var i = 0; i < FIELD_KEYS.length; i++) {
      fieldRe[FIELD_KEYS[i]] = new Float64Array(n);
      fieldIm[FIELD_KEYS[i]] = new Float64Array(n);
    }
    incidentRe = new Float64Array(n);
    incidentIm = new Float64Array(n);
    maskFlag = new Uint8Array(n);
    offscreen = {}; // 크기가 바뀌었으므로 오프스크린 캐시를 버린다
  }
  allocBuffers();

  // ---------------------------------------------------------------------
  // 좌표
  // ---------------------------------------------------------------------

  var hz = 0, hx = 0;

  // 표시 범위와 화면 격자를 현재 패널 크기·λ 에 맞춰 다시 잡는다.
  // 반환값은 "격자 칸 수가 바뀌었는가" — 바뀌었으면 장을 다시 계산해야 한다.
  function updateGeometry() {
    Z_HALF = (rangeMode === 'lambda' ? RANGE_LAMBDA_Z * lambdaOverA() : RANGE_A_Z) / 2;

    var body = document.querySelector('#panelA .body');
    var w = body ? body.clientWidth : 0;
    var h = body ? body.clientHeight : 0;
    var changed = false;

    if (w > 0 && h > 0) {
      var aspect = w / h;
      var nx = Math.round(Math.sqrt(TOTAL_SAMPLES / aspect) / 2) * 2; // 짝수
      if (nx < 8) nx = 8;
      var nz = Math.max(16, Math.round(nx * aspect));
      if (nz !== NZ || nx !== NX) {
        NZ = nz; NX = nx;
        allocBuffers();
        changed = true;
      }
    }

    // 픽셀이 정사각형이 되도록 x 범위를 z 범위에서 유도한다
    X_HALF = (Z_HALF * NX) / NZ;
    hz = (2 * Z_HALF) / NZ;
    hx = (2 * X_HALF) / NX;
    return changed;
  }

  // 셀 중심 배치 — x 표본이 x=0 에 대해 대칭이어야 좌우 반전이 정확해진다
  // (1단계 물리 격자가 셀 중심으로 Im(F)/N=0 을 보장했던 것과 같은 논리).
  // x=0 표본은 존재하지 않고, 행 j 와 NX−1−j 가 정확한 거울쌍이 된다.
  // 화면 좌표는 계속 "a 단위"다(Z_HALF = 5 는 5a 라는 뜻). 물리 함수를 부를
  // 때만 aUm 을 곱해 μm 로 바꾼다 — 눈금 막대·경계원 같은 오버레이가 전부
  // a 단위로 그려지므로, 화면 쪽을 실단위로 바꾸면 그쪽이 전부 흔들린다.
  function zAt(i) { return -Z_HALF + (i + 0.5) * hz; }
  function xAt(j) { return X_HALF - (j + 0.5) * hx; } // 행 0 이 화면 위쪽(+x)

  // ---------------------------------------------------------------------
  // 장 계산
  // ---------------------------------------------------------------------

  function compute() {
    var k = waveK();
    var alpha = alphaTotal();
    var halfRows = NX / 2; // 위쪽 절반이 x>0, 나머지는 거울상
    var maxAmp = 0;

    for (var j = 0; j < halfRows; j++) {
      var x = xAt(j);
      var jm = NX - 1 - j; // 거울 행
      for (var i = 0; i < NZ; i++) {
        var z = zAt(i);
        var idx = j * NZ + i;
        var idxm = jm * NZ + i;

        // 입사파 e^{ikz} — 복사파와 같은 규약으로 실·허수부를 함께 둔다.
        var ph = k * z * aUm;
        var cr = Math.cos(ph), ci = Math.sin(ph);
        incidentRe[idx] = cr; incidentRe[idxm] = cr; // 입사파는 x에 무관
        incidentIm[idx] = ci; incidentIm[idxm] = ci;

        // 마스크 안쪽은 계산 자체를 생략한다(Born 근사 적용 범위 밖 + 성능).
        if (x * x + z * z < MASK_R * MASK_R) {
          maskFlag[idx] = 1;
          maskFlag[idxm] = 1;
          continue;
        }
        maskFlag[idx] = 0;
        maskFlag[idxm] = 0;

        var f = RGD.scatteredFieldXZPlane(x * aUm, z * aUm, grid, k, alpha);
        var fr = 0, fi = 0;
        for (var t = 0; t < 3; t++) {
          var tn = TERM_NAMES[t];
          var vr = f[tn].ex.re, vi = f[tn].ex.im;
          fr += vr; fi += vi;
          // x 거울 관계: E_x(−x_o, z_o) = +E_x(x_o, z_o) (짝함수).
          // 커널이 dx에 대해 E_x는 짝·E_z는 홀이라는 성질과, 격자가 x→−x
          // 대칭이라는 조건이 함께 있어야 성립한다. 표시량이 E_x이므로
          // 부호 그대로 복사한다(E_z를 쓰게 되면 부호를 뒤집어야 한다).
          fieldRe[tn][idx] = vr; fieldRe[tn][idxm] = vr;
          fieldIm[tn][idx] = vi; fieldIm[tn][idxm] = vi;
        }
        fieldRe.full[idx] = fr; fieldRe.full[idxm] = fr;
        fieldIm.full[idx] = fi; fieldIm.full[idxm] = fi;

        var amp = Math.sqrt(fr * fr + fi * fi);
        if (amp > maxAmp) maxAmp = amp;
      }
    }

    metrics.maxScaAmp = maxAmp; // 입사파 진폭 1 대비
    metrics.ratioAtLambda = scatteredRatioAtLambda();
    if (DEBUG_SYMMETRY) assertMirror(k);
    computeMetrics(k);
    computeArrows(k, alpha);
  }

  // --- 축 위 화살표 ------------------------------------------------------
  // 색 지도와 반드시 같은 격자(ND_DISPLAY)로 계산한다 — physGrid(nd=20)로
  // 계산하면 밑에 깔린 색과 값이 미세하게 어긋나 겹쳐 보였을 때 어색해진다.
  //
  // 화면 x 표본은 셀 중심이라 x=0 표본이 존재하지 않는다. 그래서 격자에서
  // 읽지 않고 x=0 에서 scatteredFieldXZPlane 을 직접 부른다
  // (40점 × 열 축약 ≈ 3만 — 무시할 수준).
  var ARROW_N = 40;
  var arrows = null;

  function computeArrows(k, alpha) {
    var zh = Z_HALF;
    var step = (2 * zh) / ARROW_N;
    var z = [], iRe = [], iIm = [];
    var sRe = {}, sIm = {};
    for (var q = 0; q < FIELD_KEYS.length; q++) { sRe[FIELD_KEYS[q]] = []; sIm[FIELD_KEYS[q]] = []; }

    for (var i = 0; i < ARROW_N; i++) {
      var zz = -zh + (i + 0.5) * step;
      if (Math.abs(zz) < MASK_R) continue; // 마스크 안은 생략
      var f = RGD.scatteredFieldXZPlane(0, zz * aUm, grid, k, alpha);
      var ph = k * zz * aUm;
      z.push(zz);
      iRe.push(Math.cos(ph));
      iIm.push(Math.sin(ph));
      var fr = 0, fi = 0;
      for (var t = 0; t < 3; t++) {
        var tn = TERM_NAMES[t];
        sRe[tn].push(f[tn].ex.re);
        sIm[tn].push(f[tn].ex.im);
        fr += f[tn].ex.re; fi += f[tn].ex.im;
      }
      sRe.full.push(fr); sIm.full.push(fi);
    }
    arrows = { z: z, iRe: iRe, iIm: iIm, sRe: sRe, sIm: sIm };
  }

  // 미러링으로 채운 행 하나를 직접 계산해 대조한다.
  function assertMirror(k) {
    var j = NX - 1 - 8; // 거울로 채워진 행
    var x = xAt(j);
    var worst = 0;
    for (var i = 0; i < NZ; i += 7) {
      var idx = j * NZ + i;
      if (maskFlag[idx]) continue;
      var f = RGD.scatteredFieldXZPlane(x * aUm, zAt(i) * aUm, grid, k, alphaTotal());
      for (var t = 0; t < 3; t++) {
        var tn = TERM_NAMES[t];
        var direct = f[tn].ex.re;
        var scale = Math.abs(direct);
        if (scale < 1e-300) continue;
        var dev = Math.abs(direct - fieldRe[tn][idx]) / scale;
        if (dev > worst) worst = dev;
      }
    }
    if (worst > 1e-10) {
      console.warn('DEBUG_SYMMETRY: 미러링 불일치 상대편차 ' + worst.toExponential(3));
    }
  }

  // (C) 고정 기준점 지표 — +z축 r=2a 에서 세 항의 |E|/|E_inc|.
  //
  // 반드시 전체 다이폴 합으로 계산한다. dipoleFieldKernel 은 쌍극자 하나의
  // 커널이므로 그것을 한 번 호출해 지표를 만들면 (C)가 화면에 그리는 전체
  // 합과 어긋난다(x≪1에서는 우연히 비슷하지만 x=2에서는 exp(ikz_j) 위상
  // 변화 때문에 명확히 달라진다). 단일 점이므로 비용은 무시할 수준이다.
  function computeMetrics(k) {
    for (var t = 0; t < 3; t++) {
      var tn = TERM_NAMES[t];
      var f = RGD.scatteredFieldXZPlane(0, 2 * aUm, grid, k, alphaTotal(), { term: tn });
      metrics[tn] = Math.sqrt(
        f.ex.re * f.ex.re + f.ex.im * f.ex.im + f.ez.re * f.ez.re + f.ez.im * f.ez.im
      );
    }
  }

  // ---------------------------------------------------------------------
  // 색
  // ---------------------------------------------------------------------

  // diverging (파랑 ← 흰색 → 빨강). 색 스케일은 입사파 진폭 1로 고정하고
  // (A)(B)(C) 세 패널이 공유한다 — (B)만 따로 자동 스케일하면 진폭 대비가
  // 사라지기 때문이다.
  var NEG = [33, 102, 172];
  var MID = [247, 247, 247];
  var POS = [178, 24, 43];
  var MASK_RGB = [176, 176, 176];

  function colorInto(out, o, v) {
    if (v > 1) v = 1; else if (v < -1) v = -1;
    var a, b, f;
    if (v >= 0) { a = MID; b = POS; f = v; } else { a = MID; b = NEG; f = -v; }
    out[o] = a[0] + (b[0] - a[0]) * f;
    out[o + 1] = a[1] + (b[1] - a[1]) * f;
    out[o + 2] = a[2] + (b[2] - a[2]) * f;
    out[o + 3] = 255;
  }

  // ---------------------------------------------------------------------
  // 렌더링
  // ---------------------------------------------------------------------

  var offscreen = {};
  function getOffscreen(key) {
    if (!offscreen[key]) {
      var c = document.createElement('canvas');
      c.width = NZ;
      c.height = NX;
      offscreen[key] = { canvas: c, ctx: c.getContext('2d'), image: c.getContext('2d').createImageData(NZ, NX) };
    }
    return offscreen[key];
  }

  // valueAt(idx) → 표시값. 마스크는 회색으로 덮는다.
  function paint(key, valueAt) {
    var off = getOffscreen(key);
    var data = off.image.data;
    for (var p = 0; p < NZ * NX; p++) {
      var o = p * 4;
      if (maskFlag[p]) {
        data[o] = MASK_RGB[0]; data[o + 1] = MASK_RGB[1]; data[o + 2] = MASK_RGB[2]; data[o + 3] = 255;
      } else {
        colorInto(data, o, valueAt(p));
      }
    }
    off.ctx.putImageData(off.image, 0, 0);
    return off.canvas;
  }

  // 표시 캔버스를 패널 크기에 맞추고, 종횡비를 유지한 레터박스로 그린다.
  // (CSS 로 늘리면 10a×6a 비율이 깨져 장 분포가 왜곡된다.)
  function blit(canvas, source) {
    var body = canvas.parentNode;
    var dpr = window.devicePixelRatio || 1;
    var cw = body.clientWidth;
    var ch = body.clientHeight;
    if (cw <= 0 || ch <= 0) return null;

    canvas.width = Math.round(cw * dpr);
    canvas.height = Math.round(ch * dpr);
    canvas.style.width = cw + 'px';
    canvas.style.height = ch + 'px';

    var ctx = canvas.getContext('2d');
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, cw, ch);

    var scale = Math.min(cw / NZ, ch / NX);
    var dw = NZ * scale;
    var dh = NX * scale;
    var dx = (cw - dw) / 2;
    var dy = (ch - dh) / 2;

    ctx.imageSmoothingEnabled = true;
    ctx.drawImage(source, dx, dy, dw, dh);

    ctx.strokeStyle = '#bbb';
    ctx.lineWidth = 1;
    ctx.strokeRect(dx + 0.5, dy + 0.5, dw - 1, dh - 1);

    // 월드 좌표 → 화면 좌표 변환기를 함께 돌려준다(눈금 막대·기준점 표시용)
    return {
      ctx: ctx,
      zToPx: function (z) { return dx + ((z + Z_HALF) / (2 * Z_HALF)) * dw; },
      xToPx: function (x) { return dy + ((X_HALF - x) / (2 * X_HALF)) * dh; },
      box: { x: dx, y: dy, w: dw, h: dh }
    };
  }

  // (A) 파장 눈금 막대. λ 기준 축척에서는 화면 폭의 정확히 1/3 이 된다.
  // 'a 기준'으로 바꾸면 λ가 창보다 길 때 화살표가 패널 밖으로 잘려 나가는데,
  // 그 잘림 자체가 "파장이 이 화면보다 크다"는 메시지다.
  function drawWavelengthBar(v) {
    var ctx = v.ctx;
    var loa = lambdaOverA(); // 화면 좌표가 a 단위이므로 막대 길이도 λ/a
    var y = v.box.y + v.box.h * 0.87;
    var x0 = v.zToPx(-loa / 2);
    var x1 = v.zToPx(loa / 2);

    ctx.save();
    ctx.beginPath();
    ctx.rect(v.box.x, v.box.y, v.box.w, v.box.h);
    ctx.clip(); // 패널 경계에서 잘린다

    ctx.strokeStyle = '#111';
    ctx.fillStyle = '#111';
    ctx.lineWidth = 1.5;

    ctx.beginPath();
    ctx.moveTo(x0, y);
    ctx.lineTo(x1, y);
    ctx.stroke();

    var head = 6;
    [[x0, 1], [x1, -1]].forEach(function (e) {
      ctx.beginPath();
      ctx.moveTo(e[0], y);
      ctx.lineTo(e[0] + e[1] * head, y - head * 0.6);
      ctx.lineTo(e[0] + e[1] * head, y + head * 0.6);
      ctx.closePath();
      ctx.fill();
    });

    // 세로 끝단 표시
    ctx.beginPath();
    ctx.moveTo(x0, y - head); ctx.lineTo(x0, y + head);
    ctx.moveTo(x1, y - head); ctx.lineTo(x1, y + head);
    ctx.stroke();

    var label = 'λ = ' + fmtLen(lambdaUm);
    ctx.font = '600 12px Consolas, monospace';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'bottom';
    var tw = ctx.measureText(label).width;
    var cx = Math.min(Math.max((x0 + x1) / 2, v.box.x + tw / 2 + 4), v.box.x + v.box.w - tw / 2 - 4);
    ctx.fillStyle = 'rgba(255,255,255,0.85)';
    ctx.fillRect(cx - tw / 2 - 3, y - 20, tw + 6, 15);
    ctx.fillStyle = '#111';
    ctx.fillText(label, cx, y - 6);

    ctx.restore();
  }

  // (B) kR = 1 경계원 — 반지름 R = 1/k = λ/2π.
  // 원 안은 반응성 근접 영역(kR<1), 밖은 복사 영역(kR>1).
  // 원이 화면 밖이거나 마스크 안이면 원 대신 배지로 알린다.
  function cornerR() { return Math.sqrt(Z_HALF * Z_HALF + X_HALF * X_HALF); }

  function drawBoundaryCircle(v, k) {
    var R = 1 / (k * aUm); // a 단위. 1/k [μm] 를 a 로 나눈 값 = 1/x
    var badge = document.getElementById('zoneBadge');

    if (R > cornerR()) {
      badge.textContent = '이 화면 전체가 반응성 근접 영역 (kR < 1)';
      badge.hidden = false;
      return;
    }
    if (R < MASK_R) {
      badge.textContent = '이 화면 전체가 복사 영역 (kR > 1)';
      badge.hidden = false;
      return;
    }
    badge.hidden = true;

    var ctx = v.ctx;
    var cx = v.zToPx(0);
    var cy = v.xToPx(0);
    var rpx = v.zToPx(R) - cx; // z·x 스케일이 같으므로(픽셀 정사각형) 한쪽만 재면 된다

    ctx.save();
    ctx.beginPath();
    ctx.rect(v.box.x, v.box.y, v.box.w, v.box.h);
    ctx.clip();

    ctx.strokeStyle = '#2f6f9f';
    ctx.lineWidth = 1.5;
    ctx.setLineDash([5, 4]);
    ctx.beginPath();
    ctx.arc(cx, cy, rpx, 0, 2 * Math.PI);
    ctx.stroke();
    ctx.setLineDash([]);

    var label = 'kR = 1';
    ctx.font = '600 11px Consolas, monospace';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    var tw = ctx.measureText(label).width;
    var ly = cy - rpx;
    ctx.fillStyle = 'rgba(255,255,255,0.9)';
    ctx.fillRect(cx - tw / 2 - 3, ly - 7, tw + 6, 14);
    ctx.fillStyle = '#2f6f9f';
    ctx.fillText(label, cx, ly);
    ctx.restore();
  }

  // 기준점 r=2a (+z축) 표시 — 지표가 어느 지점 값인지 화면에 명기한다.
  // 지표를 (B)로 옮겼으므로 이 마커도 (B) 이미지에 그린다.
  function drawReferencePoint(v) {
    var ctx = v.ctx;
    var px = v.zToPx(2);
    var py = v.xToPx(0);
    ctx.save();
    ctx.strokeStyle = '#111';
    ctx.fillStyle = '#111';
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    ctx.arc(px, py, 4, 0, 2 * Math.PI);
    ctx.stroke();
    ctx.font = '600 11px Consolas, monospace';
    ctx.textAlign = 'left';
    ctx.textBaseline = 'bottom';
    ctx.fillStyle = 'rgba(255,255,255,0.85)';
    ctx.fillRect(px + 5, py - 15, 62, 13);
    ctx.fillStyle = '#111';
    ctx.fillText('기준 r=2a', px + 7, py - 4);
    ctx.restore();
  }

  // 입자 표시. λ 기준 축척에서는 λ/a 가 커질수록 몇 픽셀짜리 점이 되므로
  // 최소 2px(반지름 1px)는 그려서 입자 위치가 늘 보이게 한다.
  function drawParticle(v) {
    var ctx = v.ctx;
    var px = v.zToPx(0);
    var py = v.xToPx(0);
    var rpx = Math.max(v.zToPx(MASK_R) - px, 1);

    ctx.save();
    ctx.beginPath();
    ctx.arc(px, py, rpx, 0, 2 * Math.PI);
    ctx.fillStyle = 'rgb(176,176,176)';
    ctx.fill();
    if (rpx >= 3) {
      ctx.strokeStyle = '#8c8c8c';
      ctx.lineWidth = 1;
      ctx.stroke();
    }
    // 안내 문구는 원이 충분히 클 때만 (작으면 글자가 원 밖으로 삐져나온다)
    if (rpx >= 34) {
      ctx.font = '10px "Malgun Gothic", sans-serif';
      ctx.textAlign = 'center';
      ctx.textBaseline = 'middle';
      ctx.fillStyle = '#fff';
      ctx.fillText('모형 적용', px, py - 6);
      ctx.fillText('범위 밖', px, py + 6);
    }
    ctx.restore();
  }

  // 축 위 화살표 — Matter & Interactions 24.38 과 같은 읽기.
  // 세 줄이 같은 스케일을 쓴다: |E_입사| = 1 이 행 높이의 35%.
  // (B)만은 색 지도와 어긋나지 않도록 표시 배율을 함께 건다.
  var ARROW_UNIT_FRAC = 0.35;

  function drawArrows(v, kind, mag) {
    if (!showArrows || !arrows) return;
    var ctx = v.ctx;
    var c = Math.cos(animPhase), s = Math.sin(animPhase);
    var unit = ARROW_UNIT_FRAC * v.box.h * (mag || 1);
    var y0 = v.xToPx(0);
    var re, im;

    if (kind === 'inc') { re = arrows.iRe; im = arrows.iIm; }
    else { re = arrows.sRe[kind === 'sca' ? term : 'full']; im = arrows.sIm[kind === 'sca' ? term : 'full']; }

    ctx.save();
    ctx.beginPath();
    ctx.rect(v.box.x, v.box.y, v.box.w, v.box.h);
    ctx.clip();
    ctx.strokeStyle = 'rgba(60,60,60,0.8)';
    ctx.fillStyle = 'rgba(60,60,60,0.8)';
    ctx.lineWidth = 1;

    for (var i = 0; i < arrows.z.length; i++) {
      var val = re[i] * c + im[i] * s;
      if (kind === 'tot') val += arrows.iRe[i] * c + arrows.iIm[i] * s;
      var px = v.zToPx(arrows.z[i]);
      var len = val * unit;
      if (Math.abs(len) < 0.6) continue; // 눈에 안 보이는 것은 그리지 않는다
      var y1 = y0 - len;
      ctx.beginPath();
      ctx.moveTo(px, y0);
      ctx.lineTo(px, y1);
      ctx.stroke();
      var hd = Math.min(4, Math.abs(len) * 0.5);
      var dir = len > 0 ? 1 : -1;
      ctx.beginPath();
      ctx.moveTo(px, y1);
      ctx.lineTo(px - hd * 0.7, y1 + dir * hd);
      ctx.lineTo(px + hd * 0.7, y1 + dir * hd);
      ctx.closePath();
      ctx.fill();
    }
    ctx.restore();
  }

  // 색 위에 얹는 작은 라벨. 흰 배경을 깔지 않으면 빨강·파랑 위에서 읽히지 않는다.
  function labelBox(ctx, text, x, y, align, color) {
    var w = ctx.measureText(text).width;
    ctx.fillStyle = 'rgba(255,255,255,0.82)';
    ctx.fillRect(align === 'right' ? x - w - 3 : x - 3, y - 10, w + 6, 12);
    ctx.fillStyle = color;
    ctx.textAlign = align;
    ctx.textBaseline = 'bottom';
    ctx.fillText(text, x, y);
  }

  // 좌하단 좌표 아이콘 — 세 줄 모두에 같은 자리에 둔다
  function drawAxisIcon(v) {
    var ctx = v.ctx;
    ctx.save();
    ctx.font = '10px Consolas, monospace';
    labelBox(ctx, '↑x  →z  (y는 화면 안쪽)', v.box.x + 7, v.box.y + v.box.h - 5, 'left', '#555');
    ctx.restore();
  }

  // 색 규약 — (A)에만 한 번 둔다
  function drawColorNote(v) {
    var ctx = v.ctx;
    ctx.save();
    ctx.font = '10px "Malgun Gothic", sans-serif';
    var rx = v.box.x + v.box.w - 7;
    labelBox(ctx, '빨강 = +x 전기장', rx, v.box.y + v.box.h - 18, 'right', '#c0392b');
    labelBox(ctx, '파랑 = −x 전기장', rx, v.box.y + v.box.h - 5, 'right', '#2166ac');
    ctx.restore();
  }

  // ---------------------------------------------------------------------
  // 표시 배율 (B 전용)
  // ---------------------------------------------------------------------

  // 기본은 ×1 이다. 자동 배율이 켜져 있으면 λ를 늘려도 (B)가 늘 비슷하게
  // 보여서 "산란파가 약해진다"는 메시지가 통째로 사라진다.
  function currentScale() {
    if (scaleMode !== 'auto') return parseFloat(scaleMode);
    var max = 0;
    for (var p = 0; p < NZ * NX; p++) {
      if (maskFlag[p]) continue;
      var re = fieldRe[term][p], im = fieldIm[term][p];
      var v = Math.sqrt(re * re + im * im);
      if (v > max) max = v;
    }
    if (max <= 0) return 1;
    // 최댓값이 색 스케일 상한(1)의 80% 정도가 되게 하되, 1,10,100… 눈금으로 맞춘다
    var raw = 0.8 / max;
    if (raw <= 1) return 1;
    return Math.pow(10, Math.floor(Math.log(raw) / Math.LN10));
  }

  // ---------------------------------------------------------------------
  // 그리기
  // ---------------------------------------------------------------------

  // 숨은 탭의 캔버스는 clientWidth = 0 이라 그릴 수 없다. 활성 탭만 그리고,
  // 탭을 전환할 때 그 탭을 다시 그린다.
  function render() {
    if (activeTab === 1) renderTab1();
    else renderTab2();
  }

  function renderTab2() {
    var k = waveK();
    renderPhasorSum(k);
    renderScatterPattern();
  }

  // 세 줄은 같은 z축·같은 색 스케일·같은 픽셀 크기를 공유한다.
  // 프레임마다 하는 일은 Re·cos(ωt) + Im·sin(ωt) 곱셈뿐이다.
  function renderTab1() {
    var k = waveK();
    var c = Math.cos(animPhase), s = Math.sin(animPhase);

    var srcA = paint('A', function (p) { return incidentRe[p] * c + incidentIm[p] * s; });
    var vA = blit(document.getElementById('canvasA'), srcA);
    if (vA) {
      drawParticle(vA);
      drawWavelengthBar(vA);
      drawArrows(vA, 'inc', 1);
      drawAxisIcon(vA);
      drawColorNote(vA);
    }

    var mag = currentScale();
    var fr = fieldRe[term], fi = fieldIm[term];
    var srcB = paint('B', function (p) { return (fr[p] * c + fi[p] * s) * mag; });
    var vB = blit(document.getElementById('canvasB'), srcB);
    if (vB) {
      drawParticle(vB);
      drawBoundaryCircle(vB, k);
      drawArrows(vB, 'sca', mag);
      drawAxisIcon(vB);
      if (showMetrics) drawReferencePoint(vB);
    }

    // (C)는 항상 세 항을 모두 합친 원본값(배율 1)으로 계산한다.
    // 배율이 걸린 산란파를 더하면 (C)가 물리적으로 무의미해진다.
    var cr = fieldRe.full, ci = fieldIm.full;
    var srcC = paint('C', function (p) {
      return (incidentRe[p] + cr[p]) * c + (incidentIm[p] + ci[p]) * s;
    });
    var vC = blit(document.getElementById('canvasC'), srcC);
    if (vC) {
      drawParticle(vC);
      drawArrows(vC, 'tot', 1);
      drawAxisIcon(vC);
    }

    document.getElementById('scaleBadge').textContent = '×' + fmtInt(mag);
    document.getElementById('scaleBadge').hidden = (mag === 1);
    document.getElementById('termLabel').textContent = TERM_LABELS[term];

    // 산란파가 사실상 백지로 보이면 고장으로 오인되므로 알린다.
    // 배지가 뜬다는 사실 자체가 전달하려는 내용이다.
    //
    // 기준은 화면 최댓값이 아니라 z = λ 지점의 비다. 화면 최댓값은 마스크
    // 가장자리의 1/R³ 근접장이 늘 0.1 수준이라(λ와 거의 무관) 백지 여부를
    // 전혀 가려내지 못한다. z = λ 는 항상 kr = 2π 라 λ와 무관하게 비교된다.
    var weak = document.getElementById('weakBadge');
    var rel = metrics.ratioAtLambda;
    if (rel * mag < 0.02) {
      weak.textContent = '산란파 = 입사파의 ' + (rel * 100).toPrecision(2) +
        '% (z = λ) — 배율을 올려 보세요';
      weak.hidden = false;
    } else {
      weak.hidden = true;
    }

    renderMetrics();
    renderD(k);
  }

  // 지표 줄 — (B) 이미지 위의 오버레이. 세 비는 모두 "산란장의 구성"에 관한
  // 값이므로 (B)에 속한다((B)의 항별 토글을 숫자가 설명해 준다).
  //   기준 r=2a │ 머무는 장: 정전 … · 유도 … │ 나가는 장: 복사 …
  function renderMetrics() {
    var bar = document.getElementById('metricsBar');
    bar.hidden = !showMetrics;
    if (!showMetrics) return;
    bar.innerHTML =
      '<span class="m-ref">기준 r=2a</span><span class="m-sep">│</span>' +
      '<span class="m-stay">머무는 장: 정전 <b>' + metrics.static.toExponential(2) +
      '</b> · 유도 <b>' + metrics.induction.toExponential(2) + '</b></span>' +
      '<span class="m-sep">│</span>' +
      '<span class="m-go">나가는 장: 복사 <b>' + metrics.radiation.toExponential(2) + '</b></span>';
  }

  // ---------------------------------------------------------------------
  // 2행 — 기제 패널 (D)(E)(F)
  //
  // 정사각형 좌표계를 셀 안에 중앙 정렬해 letterbox 한다. 종횡비가 1:1이
  // 아니면 (D)의 원, (E)의 복소평면, (F)의 극좌표가 왜곡되어 논문 표 1의
  // 코르뉘 나선과 대응이 깨진다.
  // ---------------------------------------------------------------------

  // 데이터 범위 [-h, h] 를 정사각형으로 그린다. 중심은 (cx0, cy0).
  function blitSquare(canvas, cx0, cy0, halfSpan) {
    var body = canvas.parentNode;
    var dpr = window.devicePixelRatio || 1;
    var cw = body.clientWidth;
    var ch = body.clientHeight;
    if (cw <= 0 || ch <= 0) return null;

    canvas.width = Math.round(cw * dpr);
    canvas.height = Math.round(ch * dpr);
    canvas.style.width = cw + 'px';
    canvas.style.height = ch + 'px';

    var ctx = canvas.getContext('2d');
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, cw, ch);

    var side = Math.min(cw, ch);       // 셀보다 작아도 무방 — 정사각형 유지가 우선
    var ox = (cw - side) / 2;
    var oy = (ch - side) / 2;
    var s = side / (2 * halfSpan);

    ctx.strokeStyle = '#e4e4e4';
    ctx.lineWidth = 1;
    ctx.strokeRect(ox + 0.5, oy + 0.5, side - 1, side - 1);

    return {
      ctx: ctx,
      side: side,
      s: s,
      X: function (u) { return ox + side / 2 + (u - cx0) * s; },
      Y: function (v) { return oy + side / 2 - (v - cy0) * s; }
    };
  }

  // --- (D) 내부 위상 지도 -------------------------------------------------
  // columns[{x, z, m}] 산점도. 색 = cos(k·z), 크기 = m.
  // 색 줄무늬 개수 — (D)의 핵심 메시지를 숫자로 못박는다.
  // 구동 위상은 cos(k z) 이고 입자 안에서 z ∈ [−a, a] 이므로 위상은 2ka 만큼
  // 변한다. 그 구간에 들어가는 반주기(색이 한 번 뒤집히는 단위) 수를 센다.
  function renderStripes(k) {
    var el = document.getElementById('stripeReadout');
    if (!el) return;
    var halfPeriods = (2 * k * aUm) / Math.PI; // 2ka / π
    var txt, note;
    if (halfPeriods < 0.15) {
      txt = '거의 단색';
      note = '입자 전체가 같은 위상으로 진동 — 논문 Ⅲ.2가 성립하는 영역';
    } else if (halfPeriods < 1) {
      txt = '줄무늬 1개 미만';
      note = '위상 균일성이 흔들리기 시작하는 지점';
    } else {
      txt = '색 줄무늬 약 ' + halfPeriods.toFixed(1) + '개';
      note = '입자 안에서 위상이 뒤집힌다 — "전체가 같은 위상"이 무너진 영역';
    }
    el.innerHTML = '<b>' + txt + '</b><span class="note">' + note + '</span>';
  }

  function renderD(k) {
    renderStripes(k);
    var v = blitSquare(document.getElementById('canvasD'), 0, 0, 1.15);
    if (!v) return;
    var ctx = v.ctx;

    // 입자 경계
    ctx.strokeStyle = '#bbb';
    ctx.setLineDash([4, 3]);
    ctx.beginPath();
    ctx.arc(v.X(0), v.Y(0), v.s, 0, 2 * Math.PI);
    ctx.stroke();
    ctx.setLineDash([]);

    var cols = physGrid.columns;
    var rgb = [0, 0, 0, 0];
    // 마커 반지름: 격자 간격(2a/nd = 0.1a)의 절반을 넘지 않아야 이웃과 겹치지
    // 않는다. 겹치면 m(열의 쌍극자 수) 크기 부호화가 사라지고 통짜 원이 된다.
    var rMax = 0.9 * (1 / ND_PHYS) * v.s;
    for (var c = 0; c < cols.length; c++) {
      var col = cols[c];
      // (A)(B)(C)와 같은 위상으로 움직인다 — 입사파 마루가 입자를 지날 때
      // 그 자리 전하가 함께 변하는 것이 보여야 한다.
      colorInto(rgb, 0, Math.cos(k * col.z - animPhase));
      ctx.fillStyle = 'rgb(' + Math.round(rgb[0]) + ',' + Math.round(rgb[1]) + ',' + Math.round(rgb[2]) + ')';
      ctx.beginPath();
      // 넓이가 m에 비례하도록 반지름은 √m
      // col 좌표는 μm 다. 이 패널의 좌표계는 a 단위(halfSpan = 1.15a)이므로 나눈다.
      ctx.arc(v.X(col.z / aUm), v.Y(col.x / aUm), rMax * Math.sqrt(col.m / 20), 0, 2 * Math.PI);
      ctx.fill();
    }

    // 축 라벨 (z = 진행 방향, 1행과 같은 방향으로 맞춘다)
    ctx.fillStyle = '#888';
    ctx.font = '10px Consolas, monospace';
    ctx.textAlign = 'right';
    ctx.textBaseline = 'bottom';
    ctx.fillText('+z →', v.X(1.1), v.Y(-1.02));
    ctx.textAlign = 'left';
    ctx.fillText('↑ x (편광)', v.X(-1.12), v.Y(0.95));
  }

  // --- (E) 위상자 합 ------------------------------------------------------
  // exp(i q·r_j) 를 이어 붙인 누적합 경로.
  // 그리는 순서는 q̂·r_j 오름차순 — 벡터 합은 교환법칙이 성립하므로 결과는
  // 순서와 무관하지만, 순서가 없으면 그림이 낙서가 된다. 위치 순으로 이어야
  // 위상이 단조 증가해 나선이 된다.
  var walkOrder = null;      // 정렬된 인덱스
  var walkOrderKey = '';     // 정렬이 유효한 (k, θ) 서명

  function phasorWalk(k, th) {
    var q = RGD.qVector(k, th, 0); // φ=0 고정
    var qn = Math.sqrt(q.x * q.x + q.z * q.z);
    var pts = physGrid.points;
    var N = physGrid.N;

    // a 가 바뀌면 physGrid 가 새로 만들어지므로 정렬 캐시도 무효가 된다.
    var key = k.toExponential(12) + '|' + th.toExponential(12) + '|' + aUm.toExponential(12);
    if (walkOrderKey !== key) {
      var idx = new Array(pts.length);
      var proj = new Float64Array(pts.length);
      for (var i = 0; i < pts.length; i++) {
        idx[i] = i;
        // θ→0 에서는 q=0 이라 방향이 정의되지 않으므로 z_j 순서로 대체한다
        proj[i] = qn < 1e-12 ? pts[i].z : (q.x * pts[i].x + q.z * pts[i].z) / qn;
      }
      idx.sort(function (a, b) { return proj[a] - proj[b]; });
      walkOrder = idx;
      walkOrderKey = key;
    }

    // 각 위상자의 길이를 1/N 로 두어 walk 전체 경로 길이가 1이 되게 한다.
    var re = 0, im = 0;
    var path = new Float64Array((pts.length + 1) * 2);
    var o = 2;
    for (var j = 0; j < walkOrder.length; j++) {
      var p = pts[walkOrder[j]];
      var ph = q.x * p.x + q.z * p.z;
      re += Math.cos(ph) / N;
      im += Math.sin(ph) / N;
      path[o] = re; path[o + 1] = im;
      o += 2;
    }
    return { path: path, re: re, im: im };
  }

  function renderPhasorSum(k) {
    // 정사각 데이터 범위: re ∈ [−0.2, 1.1], im ∈ [−0.65, 0.65] (둘 다 span 1.3)
    var v = blitSquare(document.getElementById('canvasPhasor'), 0.45, 0, 0.65);
    if (!v) return;
    var ctx = v.ctx;

    // 축
    ctx.strokeStyle = '#eee';
    ctx.beginPath();
    ctx.moveTo(v.X(-0.2), v.Y(0)); ctx.lineTo(v.X(1.1), v.Y(0));
    ctx.moveTo(v.X(0), v.Y(-0.65)); ctx.lineTo(v.X(0), v.Y(0.65));
    ctx.stroke();

    // θ=0 기준선 — 모든 위상자가 정렬되었을 때의 길이 1 직선
    ctx.strokeStyle = '#ccc';
    ctx.lineWidth = 3;
    ctx.beginPath();
    ctx.moveTo(v.X(0), v.Y(0)); ctx.lineTo(v.X(1), v.Y(0));
    ctx.stroke();

    var w = phasorWalk(k, theta);

    // walk 경로
    ctx.strokeStyle = '#2f6f9f';
    ctx.lineWidth = 1.4;
    ctx.beginPath();
    ctx.moveTo(v.X(w.path[0]), v.Y(w.path[1]));
    for (var i = 2; i < w.path.length; i += 2) {
      ctx.lineTo(v.X(w.path[i]), v.Y(w.path[i + 1]));
    }
    ctx.stroke();

    // 합 벡터
    ctx.strokeStyle = '#c0392b';
    ctx.fillStyle = '#c0392b';
    ctx.lineWidth = 2;
    ctx.beginPath();
    ctx.moveTo(v.X(0), v.Y(0));
    ctx.lineTo(v.X(w.re), v.Y(w.im));
    ctx.stroke();
    ctx.beginPath();
    ctx.arc(v.X(w.re), v.Y(w.im), 3.5, 0, 2 * Math.PI);
    ctx.fill();

    ctx.fillStyle = '#888';
    ctx.font = '10px Consolas, monospace';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'top';
    ctx.fillText('θ=0 기준선 (길이 1)', v.X(0.5), v.Y(0) + 4);

    document.getElementById('phasorTheta').textContent =
      '원거리 방향 θ = ' + Math.round((theta * 180) / Math.PI) + '° (φ=0) · |F|/N = ' +
      Math.sqrt(w.re * w.re + w.im * w.im).toFixed(4);
  }

  // --- (F) 산란 패턴 ------------------------------------------------------
  // |F|는 θ에만 의존하고 φ에는 무관하다(|q| = 2k·sin(θ/2)). 두 곡선은 여기에
  // 편광 인자만 다르게 곱한 것이다:
  //   수직면 φ=90° : sin²Θ = 1     → dσ/dΩ = |F(θ)|²        (위상자 크기의 제곱 그 자체)
  //   편광면 φ=0   : sin²Θ = cos²θ → dσ/dΩ = |F(θ)|²·cos²θ  (θ=90°에서 0)
  // 따라서 |F(θ)|² 를 한 번만 계산하고 두 인자를 곱해 쓴다. φ=0 슬라이스이므로
  // 열 축약(316항)이 그대로 쓰이며, 이는 근사가 아니라 정확한 항등식이다.
  var fCurve = new Float64Array(THETA_SAMPLES); // |F(θ)/N|²

  function computeFCurve(k) {
    var N = physGrid.N;
    for (var i = 0; i < THETA_SAMPLES; i++) {
      var th = (Math.PI * i) / (THETA_SAMPLES - 1);
      var q = RGD.qVector(k, th, 0);
      var F = RGD.formFactorPhi0Columns(physGrid.columns, q.x, q.z);
      fCurve[i] = (F.re * F.re + F.im * F.im) / (N * N);
    }
  }

  // 편광 인자 — 두 곡선은 |F(θ)|² 에 이것만 다르게 곱한 것이다
  function polPerp() { return 1; }                                  // φ=90°: sin²Θ = 1
  function polPara(th) { var c = Math.cos(th); return c * c; }      // φ=0  : sin²Θ = cos²θ

  var polarMode = false; // 기본은 직교 좌표

  // 클릭·드래그로 θ를 고르기 위해, 마지막으로 그린 플롯 기하를 기억해 둔다
  var scatterGeom = null;

  function renderScatterPattern() {
    if (polarMode) renderScatterPolar();
    else renderScatterCartesian();
  }

  // 직교 좌표 (기본) — 세로축이 위상자 합의 |F|/N 제곱과 직접 대응하므로
  // 두 패널의 연결이 분명해진다.
  function renderScatterCartesian() {
    var canvas = document.getElementById('canvasScatter');
    var body = canvas.parentNode;
    var dpr = window.devicePixelRatio || 1;
    var cw = body.clientWidth, ch = body.clientHeight;
    if (cw <= 0 || ch <= 0) return;

    canvas.width = Math.round(cw * dpr);
    canvas.height = Math.round(ch * dpr);
    canvas.style.width = cw + 'px';
    canvas.style.height = ch + 'px';
    var ctx = canvas.getContext('2d');
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, cw, ch);

    var mL = 46, mR = 12, mT = 12, mB = 36;
    var pw = cw - mL - mR, phh = ch - mT - mB;
    if (pw <= 10 || phh <= 10) return;

    var TX = function (deg) { return mL + (deg / 180) * pw; };
    var TY = function (val) { return mT + (1 - val) * phh; };
    scatterGeom = { mode: 'cart', TX: TX, mL: mL, pw: pw, mT: mT, ph: phh };

    // 눈금
    ctx.strokeStyle = '#eee';
    ctx.fillStyle = '#999';
    ctx.lineWidth = 1;
    ctx.font = '10px Consolas, monospace';
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    [0, 0.25, 0.5, 0.75, 1].forEach(function (val) {
      ctx.beginPath();
      ctx.moveTo(mL, TY(val)); ctx.lineTo(mL + pw, TY(val));
      ctx.stroke();
      ctx.fillText(val.toFixed(2), mL - 5, TY(val));
    });
    ctx.textAlign = 'center';
    ctx.textBaseline = 'top';
    [0, 30, 60, 90, 120, 150, 180].forEach(function (deg) {
      ctx.beginPath();
      ctx.moveTo(TX(deg), mT); ctx.lineTo(TX(deg), mT + phh);
      ctx.stroke();
      ctx.fillText(deg + '°', TX(deg), mT + phh + 5);
    });

    // 축 이름 — θ는 공간의 점이 아니라 원거리 관측 방향이다
    ctx.fillStyle = '#666';
    ctx.font = '10.5px "Malgun Gothic", sans-serif';
    ctx.fillText('산란 방향 θ (원거리)', mL + pw / 2, mT + phh + 19);
    ctx.save();
    ctx.translate(11, mT + phh / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.textBaseline = 'top';
    ctx.fillText('dσ/dΩ (θ=0 규격화)', 0, 0);
    ctx.restore();

    function curve(pol, color, width) {
      ctx.strokeStyle = color;
      ctx.lineWidth = width;
      ctx.beginPath();
      for (var i = 0; i < THETA_SAMPLES; i++) {
        var th = (Math.PI * i) / (THETA_SAMPLES - 1);
        var px = TX((th * 180) / Math.PI), py = TY(fCurve[i] * pol(th));
        if (i === 0) ctx.moveTo(px, py); else ctx.lineTo(px, py);
      }
      ctx.stroke();
    }
    curve(polPara, '#c0392b', 1.6);
    curve(polPerp, '#2f6f9f', 2);

    // 현재 방향 마커 — 수직면 곡선 위에 둔다(그 곡선이 |F(θ)|² 그 자체이므로
    // 위상자 크기의 제곱과 직접 대응한다)
    var idx = Math.round((theta / Math.PI) * (THETA_SAMPLES - 1));
    var mx = TX((theta * 180) / Math.PI);
    ctx.strokeStyle = '#111';
    ctx.lineWidth = 1;
    ctx.setLineDash([3, 3]);
    ctx.beginPath();
    ctx.moveTo(mx, mT); ctx.lineTo(mx, mT + phh);
    ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle = '#111';
    ctx.beginPath();
    ctx.arc(mx, TY(fCurve[idx]), 4.5, 0, 2 * Math.PI);
    ctx.fill();

    drawLegend(ctx, mL + 8, mT + 6);
  }

  // 극좌표 ("방향 분포로 보기") — 정사각형이 강제된다
  function renderScatterPolar() {
    var v = blitSquare(document.getElementById('canvasScatter'), 0, 0, 1.16);
    if (!v) return;
    var ctx = v.ctx;
    var cx = v.X(0), cy = v.Y(0), R = v.s;
    scatterGeom = { mode: 'polar', cx: cx, cy: cy, R: R };

    ctx.strokeStyle = '#eee';
    ctx.lineWidth = 1;
    [0.25, 0.5, 0.75, 1].forEach(function (r) {
      ctx.beginPath();
      ctx.arc(cx, cy, R * r, 0, 2 * Math.PI);
      ctx.stroke();
    });
    for (var a = 0; a < 180; a += 30) {
      var t = (a * Math.PI) / 180;
      ctx.beginPath();
      ctx.moveTo(cx - R * Math.cos(t), cy - R * Math.sin(t));
      ctx.lineTo(cx + R * Math.cos(t), cy + R * Math.sin(t));
      ctx.stroke();
    }

    // θ는 +z(오른쪽)에서 잰다. 위/아래로 대칭 복사해 닫힌 로브를 만든다.
    function curve(pol, color, width) {
      ctx.strokeStyle = color;
      ctx.lineWidth = width;
      ctx.beginPath();
      var i, th, r;
      for (i = 0; i < THETA_SAMPLES; i++) {
        th = (Math.PI * i) / (THETA_SAMPLES - 1);
        r = R * fCurve[i] * pol(th);
        if (i === 0) ctx.moveTo(cx + r, cy); else ctx.lineTo(cx + r * Math.cos(th), cy - r * Math.sin(th));
      }
      for (i = THETA_SAMPLES - 1; i >= 0; i--) {
        th = (Math.PI * i) / (THETA_SAMPLES - 1);
        r = R * fCurve[i] * pol(th);
        ctx.lineTo(cx + r * Math.cos(th), cy + r * Math.sin(th));
      }
      ctx.closePath();
      ctx.stroke();
    }
    curve(polPara, '#c0392b', 1.6);
    curve(polPerp, '#2f6f9f', 2);

    var idx = Math.round((theta / Math.PI) * (THETA_SAMPLES - 1));
    var rm = R * fCurve[idx];
    ctx.strokeStyle = '#111';
    ctx.lineWidth = 1;
    ctx.setLineDash([3, 3]);
    ctx.beginPath();
    ctx.moveTo(cx, cy);
    ctx.lineTo(cx + R * Math.cos(theta), cy - R * Math.sin(theta));
    ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle = '#111';
    ctx.beginPath();
    ctx.arc(cx + rm * Math.cos(theta), cy - rm * Math.sin(theta), 4.5, 0, 2 * Math.PI);
    ctx.fill();

    ctx.fillStyle = '#888';
    ctx.font = '10px Consolas, monospace';
    ctx.textAlign = 'right';
    ctx.textBaseline = 'bottom';
    ctx.fillText('+z → (θ=0, 원거리 방향)', v.X(1.13), v.Y(-1.04));
    drawLegend(ctx, v.X(-1.13), v.Y(1.13));
  }

  function drawLegend(ctx, x, y) {
    ctx.font = '10px "Malgun Gothic", sans-serif';
    ctx.textAlign = 'left';
    ctx.textBaseline = 'top';
    ctx.fillStyle = '#2f6f9f';
    ctx.fillText('— 수직면 φ=90°  |F(θ)|²', x, y);
    ctx.fillStyle = '#c0392b';
    ctx.fillText('— 편광면 φ=0    |F(θ)|²·cos²θ', x, y + 13);
  }

  // ---------------------------------------------------------------------
  // 서식
  // ---------------------------------------------------------------------

  function fmtLambda(v) {
    if (v >= 1e5) return v.toExponential(2);
    if (v >= 100) return v.toFixed(0);
    return v.toFixed(2);
  }

  function fmtInt(v) {
    return v >= 1000 ? v.toExponential(0) : String(v);
  }

  // 길이는 μm 로 저장하고 표시 단위만 자동 전환한다 (nm → μm → mm → cm → m)
  function fmtLen(um) {
    var v, u;
    if (um < 1) { v = um * 1e3; u = 'nm'; }
    else if (um < 1e3) { v = um; u = 'μm'; }
    else if (um < 1e4) { v = um / 1e3; u = 'mm'; }
    else if (um < 1e6) { v = um / 1e4; u = 'cm'; }
    else { v = um / 1e6; u = 'm'; }
    var s = v >= 100 ? v.toFixed(0) : (v >= 10 ? v.toFixed(1) : v.toFixed(2));
    return s + ' ' + u;
  }

  function fmtQ(q) {
    return q >= 0.01 ? q.toFixed(2) : q.toExponential(1);
  }

  // Q_sca = σ/(πa²). σ 는 물리 격자(nd=20) 경로인 sigma1D 를 그대로 쓴다.
  // 316열 × 32노드라 드래그 중 매 입력마다 불러도 부담이 없다.
  function computeDerived() {
    var r = RGD.sigma1D(physGrid, waveK(), alphaTotal());
    return { G: r.G, Q: r.sigma / (Math.PI * aUm * aUm) };
  }

  // z = λ 지점(+z축)의 |E_산란|/|E_입사|. 그 자리는 항상 kr = 2π 라
  // λ와 무관하게 원거리 조건이 보장된다. 입사파 진폭은 1이다.
  function scatteredRatioAtLambda() {
    var f = RGD.scatteredFieldXZPlane(0, lambdaUm, grid, waveK(), alphaTotal());
    var re = 0, im = 0;
    for (var t = 0; t < 3; t++) {
      var tn = TERM_NAMES[t];
      re += f[tn].ex.re; im += f[tn].ex.im;
    }
    var ez = f.radiation.ez, ez2 = f.induction.ez, ez3 = f.static.ez;
    var zr = ez.re + ez2.re + ez3.re, zi = ez.im + ez2.im + ez3.im;
    return Math.sqrt(re * re + im * im + zr * zr + zi * zi);
  }

  function updateReadout() {
    var x = sizeX();
    var d = computeDerived();
    var ratio = scatteredRatioAtLambda();
    // Q > 1 은 버그가 아니라 RGD 근사(|m−1| ≪ 1)가 무너진다는 신호다.
    // 값을 자르거나 계수를 조정하지 않고 빨간색으로 표시만 한다.
    document.getElementById('lambdaReadout').innerHTML =
      'a = ' + fmtLen(aUm) + '<span class="sep">·</span>' +
      'λ = ' + fmtLen(lambdaUm) + '<span class="sep">·</span>' +
      'λ/a = ' + fmtLambda(lambdaOverA()) + '<span class="sep">·</span>' +
      'x = ' + (x >= 0.001 ? x.toFixed(4) : x.toExponential(3)) + '<span class="sep">·</span>' +
      '|E_산란|/|E_입사| = ' + fmtQ(ratio) + ' (z = λ)' + '<span class="sep">·</span>' +
      'Q = <b class="' + (d.Q > 1 ? 'q-hot' : 'q-ok') + '">' + fmtQ(d.Q) + '</b>' +
      '<span class="sep">·</span>G = ' + d.G.toFixed(3);
    updateNotices();
  }

  // 상태 배지 — λ/a 고정, 적용 범위 밖
  function updateNotices() {
    var lock = document.getElementById('lockNotice');
    lock.hidden = !lockRatio;
    lock.textContent = atLimit
      ? 'λ/a 고정 — 슬라이더 끝에 닿아 함께 멈춤'
      : 'λ/a 고정 — 그림이 변하지 않습니다';
    lock.classList.toggle('limit', atLimit);

    var out = lambdaOverA() < LAMBDA_OVER_A_MIN;
    document.getElementById('rangeNotice').hidden = !out;
    document.body.classList.toggle('out-of-range', out);

    document.getElementById('aSlider').classList.toggle('at-limit', atLimit);
    document.getElementById('lambdaSlider').classList.toggle('at-limit', atLimit);
  }

  // ---------------------------------------------------------------------
  // 재계산 (계산 중 표시 → 한 틱 뒤 실행)
  // ---------------------------------------------------------------------

  var status = null;

  function recompute() {
    if (!status) status = document.getElementById('status');
    status.textContent = '계산 중…';
    updateReadout();
    // file:// 에서는 Web Worker가 막히는 경우가 많아 쓰지 않는다. 동기 계산이므로
    // 브라우저가 "계산 중…"을 페인트할 틈을 주기 위해 한 틱 미룬다.
    setTimeout(function () {
      updateGeometry(); // 표시 범위는 λ에 딸려 움직인다
      compute();
      computeFCurve(waveK()); // (F) 곡선도 a·λ 에만 의존한다
      render();
      status.textContent = '';
    }, 0);
  }

  // ---------------------------------------------------------------------
  // 애니메이션
  // ---------------------------------------------------------------------

  function frame(t) {
    rafId = requestAnimationFrame(frame);
    if (!animating) { lastFrameT = t; return; }
    var dt = lastFrameT ? (t - lastFrameT) / 1000 : 0;
    lastFrameT = t;
    if (dt > 0.5) dt = 0; // 탭 복귀 등으로 크게 튄 프레임은 버린다
    animPhase = (animPhase + (2 * Math.PI) * (dt / SCREEN_PERIOD) * animSpeed) % (2 * Math.PI);
    render();
  }

  function setAnimating(on) {
    animating = on;
    lastFrameT = 0;
    var btn = document.getElementById('playBtn');
    btn.textContent = on ? '⏸ 정지' : '▶ 재생';
    btn.classList.toggle('on', !on);
    if (!on) render();
  }

  // 검증용 훅 — Playwright 에서 위상을 고정해 프레임을 비교할 때 쓴다.
  window.RGDTest = {
    freeze: function (p) { setAnimating(false); animPhase = p || 0; render(); }
  };

  // θ만 바뀔 때는 (E)(F)만 다시 그리면 된다 — 장 데이터도 (F) 곡선도 λ에만
  // 의존하므로 재계산이 필요 없다. 그래서 change가 아니라 input에서 즉시
  // 반응해도 무리가 없다.
  function renderThetaOnly() {
    renderTab2();
  }

  // ---------------------------------------------------------------------
  // 배선
  // ---------------------------------------------------------------------

  var aSlider = document.getElementById('aSlider');
  var lamSlider = document.getElementById('lambdaSlider');
  var lockBox = document.getElementById('lockBox');

  function toLog(v, min, max) { return min * Math.pow(max / min, v / 1000); }

  function fromLog(x, min, max) {
    var f = Math.log(x / min) / Math.log(max / min);
    return Math.round(Math.min(Math.max(f, 0), 1) * 1000);
  }

  function syncSliders() {
    aSlider.value = fromLog(aUm, A_MIN, A_MAX);
    lamSlider.value = fromLog(lambdaUm, LAM_MIN, LAM_MAX);
  }

  // 고정 상태에서 한쪽이 범위 끝에 닿으면 두 슬라이더를 그 자리에서 함께
  // 멈춘다 — 비율 유지가 이 체크박스의 약속이므로, 조용히 어긋나게 두지 않는다.
  // 두 범위를 λ/a 로 환산해 겹치는 구간으로 자르면 한 번에 처리된다.
  function setPair(a, lam) {
    atLimit = false;
    if (lockRatio) {
      var loA = Math.max(A_MIN, LAM_MIN / lockedRatio);
      var hiA = Math.min(A_MAX, LAM_MAX / lockedRatio);
      var ca = Math.min(Math.max(a, loA), hiA);
      // "이번에 잘렸는가"가 아니라 "결과가 경계에 앉아 있는가"로 판정한다.
      // input 다음에 오는 change 가 이미 잘린 값을 다시 넣으므로, 잘림
      // 여부로 판정하면 표시가 곧바로 지워진다.
      // a 자신의 끝은 슬라이더 위치로 이미 보이므로, λ 쪽 범위가 먼저
      // 걸려서 a 가 제 끝에 닿지 못한 경우만 알린다.
      atLimit =
        (LAM_MIN / lockedRatio > A_MIN && ca <= loA * (1 + 1e-9)) ||
        (LAM_MAX / lockedRatio < A_MAX && ca >= hiA * (1 - 1e-9));
      a = ca;
      lam = ca * lockedRatio;
    } else {
      a = Math.min(Math.max(a, A_MIN), A_MAX);
      lam = Math.min(Math.max(lam, LAM_MIN), LAM_MAX);
    }
    aUm = a;
    lambdaUm = lam;
    rebuildGrids(); // 격자는 a 에 의존한다
  }

  function driveA(a) { setPair(a, lockRatio ? a * lockedRatio : lambdaUm); }
  function driveLam(l) { setPair(lockRatio ? l / lockedRatio : aUm, l); }

  // 슬라이더는 input이 아니라 change 이벤트에서 재계산한다(기존 관례) —
  // 드래그 중 매 프레임 장 전체를 다시 도는 것을 피한다. 판독창과 짝
  // 슬라이더는 input 에서 즉시 따라간다.
  function wireSlider(el, drive, min, max) {
    el.addEventListener('input', function () {
      drive(toLog(parseFloat(el.value), min, max));
      syncSliders();
      updateReadout();
      clearPresetHighlight();
    });
    // change 에서는 슬라이더 값을 다시 읽지 않는다. 위치가 1/1000 로
    // 양자화되어 있어 되읽으면 input 이 정한 값이 미세하게 어긋나고,
    // 그 오차 때문에 "경계에 앉아 있다"는 판정이 곧바로 풀린다.
    el.addEventListener('change', function () { recompute(); });
  }
  wireSlider(aSlider, driveA, A_MIN, A_MAX);
  wireSlider(lamSlider, driveLam, LAM_MIN, LAM_MAX);

  lockBox.addEventListener('change', function () {
    lockRatio = this.checked;
    lockedRatio = lambdaUm / aUm; // 누른 시점의 비를 고정한다
    atLimit = false;
    updateNotices();
  });

  function clearPresetHighlight() {
    var btns = document.querySelectorAll('.presets button');
    for (var i = 0; i < btns.length; i++) btns[i].classList.remove('active');
  }

  var presetButtons = document.querySelectorAll('.presets button');
  for (var b = 0; b < presetButtons.length; b++) {
    presetButtons[b].addEventListener('click', function () {
      clearPresetHighlight();
      this.classList.add('active');
      var pa = parseFloat(this.getAttribute('data-a'));
      var pl = parseFloat(this.getAttribute('data-lambda'));
      // 프리셋은 두 값을 한꺼번에 정한다. 고정이 켜져 있으면 비를 새 값으로
      // 갱신한다 — 체크는 그대로 두되 프리셋과 어긋나지 않게 한다.
      lockRatio = false;
      setPair(pa, pl);
      lockRatio = lockBox.checked;
      lockedRatio = pl / pa;
      syncSliders();
      recompute();
    });
  }

  document.getElementById('termSelect').addEventListener('change', function () {
    term = this.value;
    render(); // 세 항이 캐싱되어 있으므로 재계산 없이 표시만 바꾼다
  });

  document.getElementById('scaleSelect').addEventListener('change', function () {
    scaleMode = this.value;
    render();
  });

  // 표시 범위 — λ 기준(기본) ↔ a 기준. 범위가 바뀌면 표본 위치가 달라지므로
  // 다시 그리는 것으로는 안 되고 장을 재계산해야 한다.
  document.getElementById('rangeSelect').addEventListener('change', function () {
    rangeMode = this.value;
    recompute();
  });

  document.getElementById('arrowBox').addEventListener('change', function () {
    showArrows = this.checked;
    render();
  });

  document.getElementById('playBtn').addEventListener('click', function () {
    setAnimating(!animating);
  });

  // 속도 0.25× ~ 2× (로그). 가운데(500)가 1×.
  var speedSlider = document.getElementById('speedSlider');
  speedSlider.addEventListener('input', function () {
    animSpeed = 0.25 * Math.pow(8, parseFloat(this.value) / 1000);
    document.getElementById('speedReadout').textContent = animSpeed.toFixed(2) + '×';
    if (!animating) render();
  });

  // 탭 전환 — λ/a 는 두 탭이 공유하므로 재계산 없이 다시 그리기만 한다.
  var tabButtons = document.querySelectorAll('.tab-btn');
  function setTab(n) {
    activeTab = n;
    for (var i = 0; i < tabButtons.length; i++) {
      tabButtons[i].classList.toggle('on', parseInt(tabButtons[i].getAttribute('data-tab'), 10) === n);
    }
    document.getElementById('tab1').hidden = (n !== 1);
    document.getElementById('tab2').hidden = (n !== 2);
    var groups = document.querySelectorAll('.ctl-group[data-tab]');
    for (var j = 0; j < groups.length; j++) {
      groups[j].hidden = parseInt(groups[j].getAttribute('data-tab'), 10) !== n;
    }
    // 열려 있던 ⓘ는 닫는다(다른 탭의 것이 남아 있으면 혼란스럽다)
    var overlays = document.querySelectorAll('.info-overlay');
    for (var q = 0; q < overlays.length; q++) overlays[q].hidden = true;
    var ibtns = document.querySelectorAll('.info-btn');
    for (var r = 0; r < ibtns.length; r++) ibtns[r].classList.remove('on');
    render();
  }
  for (var t = 0; t < tabButtons.length; t++) {
    tabButtons[t].addEventListener('click', function () {
      setTab(parseInt(this.getAttribute('data-tab'), 10));
    });
  }

  // 방향(θ)은 슬라이더가 아니라 산란 패턴 위를 클릭·드래그해 고른다 —
  // "왼쪽 그래프의 한 점을 고르면 오른쪽이 그 방향의 위상자 합을 보여준다"가
  // 이 탭의 핵심이므로, 고르는 동작이 그래프 위에서 일어나야 한다.
  var scatterCanvas = document.getElementById('canvasScatter');
  var dragging = false;

  function pickTheta(ev) {
    if (!scatterGeom) return;
    var r = scatterCanvas.getBoundingClientRect();
    var px = ev.clientX - r.left;
    var py = ev.clientY - r.top;
    var deg;
    if (scatterGeom.mode === 'cart') {
      deg = ((px - scatterGeom.mL) / scatterGeom.pw) * 180;
    } else {
      // 극좌표: 중심에서의 각도. 위/아래 대칭이므로 절댓값을 쓴다.
      deg = (Math.atan2(Math.abs(scatterGeom.cy - py), px - scatterGeom.cx) * 180) / Math.PI;
    }
    deg = Math.min(180, Math.max(0, deg));
    theta = (deg * Math.PI) / 180;
    renderThetaOnly();
  }

  scatterCanvas.addEventListener('mousedown', function (ev) {
    // ⓘ가 열려 있으면 그쪽이 우선
    if (!document.getElementById('infoScatter').hidden) return;
    dragging = true;
    pickTheta(ev);
    ev.preventDefault();
  });
  window.addEventListener('mousemove', function (ev) { if (dragging) pickTheta(ev); });
  window.addEventListener('mouseup', function () { dragging = false; });

  // "방향 분포로 보기" — 직교 좌표 ↔ 극좌표
  var polarBtn = document.getElementById('polarBtn');
  polarBtn.addEventListener('click', function () {
    polarMode = !polarMode;
    polarBtn.classList.toggle('on', polarMode);
    polarBtn.textContent = polarMode ? '그래프로 보기' : '방향 분포로 보기';
    renderScatterPattern();
  });

  // "수치" 토글 — 상태는 세션 동안 유지된다(λ를 바꿔도 꺼지지 않는다).
  var metricsBtn = document.getElementById('metricsBtn');
  metricsBtn.addEventListener('click', function () {
    showMetrics = !showMetrics;
    metricsBtn.classList.toggle('on', showMetrics);
    render();
  });

  // ⓘ 설명 오버레이 — 한 번에 하나만 열린다
  var infoButtons = document.querySelectorAll('.info-btn');
  for (var n = 0; n < infoButtons.length; n++) {
    infoButtons[n].addEventListener('click', function () {
      var target = document.getElementById(this.getAttribute('data-info'));
      var wasOpen = !target.hidden;
      var panels = document.querySelectorAll('.info-overlay');
      for (var q = 0; q < panels.length; q++) panels[q].hidden = true;
      for (var r = 0; r < infoButtons.length; r++) infoButtons[r].classList.remove('on');
      if (!wasOpen) {
        target.hidden = false;
        this.classList.add('on');
      }
    });
  }

  // 창 크기가 바뀌면 패널 가로세로비가 달라져 화면 격자 칸 수가 바뀔 수 있다.
  // 칸 수가 바뀌면 표본 위치가 달라지므로 다시 그리는 것만으로는 안 되고
  // 장을 재계산해야 한다.
  var resizeTimer = null;
  window.addEventListener('resize', function () {
    clearTimeout(resizeTimer);
    resizeTimer = setTimeout(function () {
      if (updateGeometry()) recompute();
      else render();
    }, 120);
  });

  // ---------------------------------------------------------------------

  syncSliders();
  setTab(1); // 탭 1이 기본. 컨트롤 그룹 표시 상태를 초기화한다
  recompute();

  // 최초 recompute 시점에 레이아웃이 아직 안 잡혀 clientWidth가 0이면 격자
  // 칸 수를 잡을 수 없다. load 후 한 번 더 재계산한다.
  window.addEventListener('load', function () {
    if (updateGeometry()) recompute(); else render();
    rafId = requestAnimationFrame(frame);
  });
})();
