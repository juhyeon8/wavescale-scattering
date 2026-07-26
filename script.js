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

  var Z_HALF = 5; // 표시 영역 z ∈ [−5a, +5a]
  var X_HALF = 3; // 표시 영역 x ∈ [−3a, +3a]
  var NZ = 120;   // 화면 격자 (z)
  var NX = 72;    // 화면 격자 (x) — 10a×6a 이므로 120:72 = 1.667 로 픽셀 정사각형

  // 마스크 반경은 1.5a로 고정한다 — 축소하지 않는다.
  // nd_display=14, 셀 중심 배치이므로 최외곽 쌍극자는 r ≈ 0.93~1.00a 에 있다.
  // 1.2a로 줄이면 최근접 쌍극자까지 0.21a 밖에 안 되어, 1/R³ 발산 때문에
  // 단일 쌍극자 하나가 전체 장의 최대 13%를 차지한다 → 격자점 위치가
  // 알갱이 무늬로 드러나 논문 그림으로 쓸 수 없다. 1.5a에서는 1~3%다.
  var MASK_R = 1.5;

  var ALPHA_TOTAL = 1;
  var LAMBDA_MIN = Math.PI; // x = 2.0000 정확히 (x ≤ 2 제한의 상한)
  var LAMBDA_MAX = 1000;

  // 파장이 표시 범위(10a)의 2배를 넘으면 입사파가 거의 단색으로 보인다.
  // 물리적으로 정확한 결과(준정전기 영역)이지만 화면만 보면 고장처럼 보이므로
  // 캡션으로 설명한다.
  var UNIFORM_CAPTION_LAMBDA = 20;

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

  var grid = RGD.buildDipoleGrid(1, ND_DISPLAY);
  var lambda = LAMBDA_MIN;
  var term = 'full';
  var scaleMode = 'auto';
  // 지표는 기본 숨김. 켠 상태는 세션 동안 유지된다(λ를 바꿔도 꺼지지 않는다).
  var showMetrics = false;

  // 화면 점마다 세 항을 캐싱한다 — (B)의 [전체]/[복사]/[유도]/[정전] 토글이
  // 재계산 없이 표시만 바꾸는 일이 되도록. 재계산은 λ가 바뀔 때만 한다.
  // 표시량이 Re(E_x) 순간 스냅샷(t=0)이므로 E_x의 실수부만 저장하면 된다.
  var fieldRe = {
    radiation: new Float64Array(NZ * NX),
    induction: new Float64Array(NZ * NX),
    static: new Float64Array(NZ * NX)
  };
  var incidentRe = new Float64Array(NZ * NX);
  var maskFlag = new Uint8Array(NZ * NX);
  var metrics = { radiation: 0, induction: 0, static: 0 };

  // ---------------------------------------------------------------------
  // 좌표
  // ---------------------------------------------------------------------

  var hz = (2 * Z_HALF) / NZ;
  var hx = (2 * X_HALF) / NX;

  // 셀 중심 배치 — x 표본이 x=0 에 대해 대칭이어야 좌우 반전이 정확해진다
  // (1단계 물리 격자가 셀 중심으로 Im(F)/N=0 을 보장했던 것과 같은 논리).
  // x=0 표본은 존재하지 않고, 행 j 와 NX−1−j 가 정확한 거울쌍이 된다.
  function zAt(i) { return -Z_HALF + (i + 0.5) * hz; }
  function xAt(j) { return X_HALF - (j + 0.5) * hx; } // 행 0 이 화면 위쪽(+x)

  // ---------------------------------------------------------------------
  // 장 계산
  // ---------------------------------------------------------------------

  function compute() {
    var k = (2 * Math.PI) / lambda;
    var halfRows = NX / 2; // 행 0..35 가 x>0, 나머지는 거울상

    for (var j = 0; j < halfRows; j++) {
      var x = xAt(j);
      var jm = NX - 1 - j; // 거울 행
      for (var i = 0; i < NZ; i++) {
        var z = zAt(i);
        var idx = j * NZ + i;
        var idxm = jm * NZ + i;

        var inc = Math.cos(k * z);
        incidentRe[idx] = inc;
        incidentRe[idxm] = inc; // 입사파는 x에 무관

        // 마스크 안쪽은 계산 자체를 생략한다(Born 근사 적용 범위 밖 + 성능).
        if (x * x + z * z < MASK_R * MASK_R) {
          maskFlag[idx] = 1;
          maskFlag[idxm] = 1;
          continue;
        }
        maskFlag[idx] = 0;
        maskFlag[idxm] = 0;

        var f = RGD.scatteredFieldXZPlane(x, z, grid, k, ALPHA_TOTAL);
        for (var t = 0; t < 3; t++) {
          var tn = TERM_NAMES[t];
          var v = f[tn].ex.re;
          fieldRe[tn][idx] = v;
          // x 거울 관계: E_x(−x_o, z_o) = +E_x(x_o, z_o) (짝함수).
          // 커널이 dx에 대해 E_x는 짝·E_z는 홀이라는 성질과, 격자가 x→−x
          // 대칭이라는 조건이 함께 있어야 성립한다. 표시량이 Re(E_x)이므로
          // 부호 그대로 복사한다(E_z를 쓰게 되면 부호를 뒤집어야 한다).
          fieldRe[tn][idxm] = v;
        }
      }
    }

    if (DEBUG_SYMMETRY) assertMirror(k);
    computeMetrics(k);
  }

  // 미러링으로 채운 행 하나를 직접 계산해 대조한다.
  function assertMirror(k) {
    var j = NX - 1 - 8; // 거울로 채워진 행
    var x = xAt(j);
    var worst = 0;
    for (var i = 0; i < NZ; i += 7) {
      var idx = j * NZ + i;
      if (maskFlag[idx]) continue;
      var f = RGD.scatteredFieldXZPlane(x, zAt(i), grid, k, ALPHA_TOTAL);
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
      var f = RGD.scatteredFieldXZPlane(0, 2, grid, k, ALPHA_TOTAL, { term: tn });
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

  // (A) 파장 눈금 막대 — 논문 표 6과 같은 형식.
  // λ가 창(10a)보다 크면 화살표가 패널 밖으로 잘려 나가게 그대로 둔다.
  // 그 잘림 자체가 "파장이 이 화면보다 크다"는 메시지다.
  function drawWavelengthBar(v) {
    var ctx = v.ctx;
    var y = v.xToPx(-2.3);
    var x0 = v.zToPx(-lambda / 2);
    var x1 = v.zToPx(lambda / 2);

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

    var label = 'λ = ' + fmtLambda(lambda) + ' a';
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
  var CORNER_R = Math.sqrt(Z_HALF * Z_HALF + X_HALF * X_HALF); // 5.83a

  function drawBoundaryCircle(v, k) {
    var R = 1 / k;
    var badge = document.getElementById('zoneBadge');

    if (R > CORNER_R) {
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

  // 마스크 안내
  function drawMaskLabel(v) {
    var ctx = v.ctx;
    var px = v.zToPx(0);
    var py = v.xToPx(0);
    ctx.save();
    ctx.font = '10px "Malgun Gothic", sans-serif';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    ctx.fillStyle = '#fff';
    ctx.fillText('모형 적용', px, py - 6);
    ctx.fillText('범위 밖', px, py + 6);
    ctx.restore();
  }

  // ---------------------------------------------------------------------
  // 표시 배율 (B 전용)
  // ---------------------------------------------------------------------

  function currentScale() {
    if (scaleMode !== 'auto') return parseFloat(scaleMode);
    var max = 0;
    for (var p = 0; p < NZ * NX; p++) {
      if (maskFlag[p]) continue;
      var v = Math.abs(termValue(p));
      if (v > max) max = v;
    }
    if (max <= 0) return 1;
    // 최댓값이 색 스케일 상한(1)의 80% 정도가 되게 하되, 1,10,100… 눈금으로 맞춘다
    var raw = 0.8 / max;
    if (raw <= 1) return 1;
    return Math.pow(10, Math.floor(Math.log(raw) / Math.LN10));
  }

  function termValue(p) {
    if (term === 'full') {
      return fieldRe.radiation[p] + fieldRe.induction[p] + fieldRe.static[p];
    }
    return fieldRe[term][p];
  }

  // ---------------------------------------------------------------------
  // 그리기
  // ---------------------------------------------------------------------

  function render() {
    var k = (2 * Math.PI) / lambda;

    var srcA = paint('A', function (p) { return incidentRe[p]; });
    var vA = blit(document.getElementById('canvasA'), srcA);
    if (vA) { drawMaskLabel(vA); drawWavelengthBar(vA); }

    var mag = currentScale();
    var srcB = paint('B', function (p) { return termValue(p) * mag; });
    var vB = blit(document.getElementById('canvasB'), srcB);
    if (vB) {
      drawMaskLabel(vB);
      drawBoundaryCircle(vB, k);
      if (showMetrics) drawReferencePoint(vB);
    }

    // (C)는 항상 세 항을 모두 합친 원본값(배율 1)으로 계산한다.
    // 배율이 걸린 산란파를 더하면 (C)가 물리적으로 무의미해진다.
    // (C)에는 이미지 + 캡션만 둔다 — 수치와 기준점 마커는 (B)로 옮겼다.
    var srcC = paint('C', function (p) {
      return incidentRe[p] + fieldRe.radiation[p] + fieldRe.induction[p] + fieldRe.static[p];
    });
    var vC = blit(document.getElementById('canvasC'), srcC);
    if (vC) drawMaskLabel(vC);

    document.getElementById('scaleBadge').textContent = '×' + fmtInt(mag);
    document.getElementById('termLabel').textContent = TERM_LABELS[term];

    // (A) 캡션 — 파장이 표시 범위보다 크면 균일하게 보이는 이유를 설명한다
    var capA = document.getElementById('captionA');
    if (lambda >= UNIFORM_CAPTION_LAMBDA) {
      capA.textContent = '파장이 표시 범위보다 커서 입사파가 균일하게 보입니다';
      capA.className = 'caption warn';
    } else {
      capA.textContent = '표시 영역 z ∈ [−5a, +5a], x ∈ [−3a, +3a]';
      capA.className = 'caption';
    }

    renderMetrics();
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

  function updateReadout() {
    var x = (2 * Math.PI) / lambda;
    document.getElementById('lambdaReadout').innerHTML =
      'λ/a = ' + fmtLambda(lambda) + '<span class="sep">·</span>x = ' +
      (x >= 0.001 ? x.toFixed(4) : x.toExponential(3));
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
      compute();
      render();
      status.textContent = '';
    }, 0);
  }

  // ---------------------------------------------------------------------
  // 배선
  // ---------------------------------------------------------------------

  var slider = document.getElementById('lambdaSlider');

  function sliderToLambda(v) {
    var f = v / 1000;
    return LAMBDA_MIN * Math.pow(LAMBDA_MAX / LAMBDA_MIN, f);
  }

  function lambdaToSlider(l) {
    var f = Math.log(l / LAMBDA_MIN) / Math.log(LAMBDA_MAX / LAMBDA_MIN);
    return Math.round(Math.min(Math.max(f, 0), 1) * 1000);
  }

  // 슬라이더는 input이 아니라 change 이벤트에서 재계산한다(기존 시뮬레이션과
  // 동일한 관례) — 드래그 중 매 프레임 302만 쌍을 다시 도는 것을 피한다.
  slider.addEventListener('input', function () {
    lambda = sliderToLambda(parseFloat(slider.value));
    updateReadout();
    clearPresetHighlight();
  });
  slider.addEventListener('change', function () {
    lambda = sliderToLambda(parseFloat(slider.value));
    recompute();
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
      lambda = parseFloat(this.getAttribute('data-lambda'));
      // 21cm 프리셋(λ/a = 2.1e6)은 슬라이더 범위 밖이므로 슬라이더는 끝에 붙는다
      slider.value = lambdaToSlider(lambda);
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

  var resizeTimer = null;
  window.addEventListener('resize', function () {
    clearTimeout(resizeTimer);
    resizeTimer = setTimeout(render, 120); // 장 데이터는 그대로, 다시 그리기만
  });

  // ---------------------------------------------------------------------

  slider.value = lambdaToSlider(lambda);
  recompute();

  // 최초 recompute 시점에 레이아웃이 아직 안 잡혀 clientWidth가 0이면 blit이
  // 그리지 않고 돌아간다. load 후 한 번 더 그려 둔다(장 데이터는 캐시되어 있으므로
  // 다시 계산하지 않는다).
  window.addEventListener('load', render);
})();
