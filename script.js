(function () {
  "use strict";

  // =====================================================================
  // 0. 상수
  // =====================================================================
  const C_LIGHT = 2.99792458e8; // m/s
  const TWO_PI = Math.PI * 2;
  const VMAX = 1.5;             // 색 포화 기준 전기장 [V/m] (세 패널 공유 고정 스케일)
  const A_RATIO_MAX = 0.30;     // 반지름 상한 = 간격의 30% (단극자 모델 유효 범위)

  // 상태 (화면 표시 단위:  a,d = mm,  λ = cm,  A = V/m)
  const state = {
    a_mm: 0.5, d_mm: 10, N: 30, lam_cm: 12.2, amp: 1.0,
    polParallel: true,   // true=E∥wire(차폐), false=E⊥wire(통과)
    playing: true,
    phase: 0,            // 애니메이션 위상 ωt [rad]
  };

  // =====================================================================
  // 1. 베셀/한켈 함수  (Abramowitz & Stegun 9.4 다항 근사)
  //    H0(x) = J0(x) + i Y0(x)   (x>0 실수)
  // =====================================================================
  function besselJ0(x) {
    const ax = Math.abs(x);
    if (ax < 3) {
      const y = (x / 3) * (x / 3);
      return 1 + y * (-2.2499997 + y * (1.2656208 + y * (-0.3163866 +
        y * (0.0444479 + y * (-0.0039444 + y * 0.0002100)))));
    }
    const z = 3 / ax;
    const f = 0.79788456 + z * (-0.00000077 + z * (-0.00552740 + z * (-0.00009512 +
      z * (0.00137237 + z * (-0.00072805 + z * 0.00014476)))));
    const t = ax - 0.78539816 + z * (-0.04166397 + z * (-0.00003954 + z * (0.00262573 +
      z * (-0.00054125 + z * (-0.00029333 + z * 0.00013558)))));
    return f / Math.sqrt(ax) * Math.cos(t);
  }
  function besselY0(x) {
    // x>0 만 사용
    if (x < 3) {
      const y = (x / 3) * (x / 3);
      const poly = 0.36746691 + y * (0.60559366 + y * (-0.74350384 + y * (0.25300117 +
        y * (-0.04261214 + y * (0.00427916 + y * (-0.00024846))))));
      return (2 / Math.PI) * Math.log(x / 2) * besselJ0(x) + poly;
    }
    const z = 3 / x;
    const f = 0.79788456 + z * (-0.00000077 + z * (-0.00552740 + z * (-0.00009512 +
      z * (0.00137237 + z * (-0.00072805 + z * 0.00014476)))));
    const t = x - 0.78539816 + z * (-0.04166397 + z * (-0.00003954 + z * (0.00262573 +
      z * (-0.00054125 + z * (-0.00029333 + z * 0.00013558)))));
    return f / Math.sqrt(x) * Math.sin(t);
  }
  // H0(x) -> {re, im}
  function hankel0(x) { return { re: besselJ0(x), im: besselY0(x) }; }

  // =====================================================================
  // 2. 복소 선형계 풀이  Z c = b   (N x N, 부분 피벗 가우스 소거)
  //    Z: Float64Array 길이 N*N*2 (행우선, [re,im] 교차),  b/x: 길이 N*2
  // =====================================================================
  function solveComplex(N, Z, b) {
    // 증강행렬 복사
    const M = Z.slice();
    const x = b.slice();
    for (let col = 0; col < N; col++) {
      // 피벗 선택 (최대 절댓값)
      let piv = col, best = -1;
      for (let r = col; r < N; r++) {
        const re = M[(r * N + col) * 2], im = M[(r * N + col) * 2 + 1];
        const mag = re * re + im * im;
        if (mag > best) { best = mag; piv = r; }
      }
      if (piv !== col) {
        for (let k = 0; k < N; k++) {
          const i1 = (col * N + k) * 2, i2 = (piv * N + k) * 2;
          let t = M[i1]; M[i1] = M[i2]; M[i2] = t;
          t = M[i1 + 1]; M[i1 + 1] = M[i2 + 1]; M[i2 + 1] = t;
        }
        let t = x[col * 2]; x[col * 2] = x[piv * 2]; x[piv * 2] = t;
        t = x[col * 2 + 1]; x[col * 2 + 1] = x[piv * 2 + 1]; x[piv * 2 + 1] = t;
      }
      // 정규화 인자
      const pr = M[(col * N + col) * 2], pi = M[(col * N + col) * 2 + 1];
      const pden = pr * pr + pi * pi;
      // 아래 행 소거
      for (let r = 0; r < N; r++) {
        if (r === col) continue;
        const fr0 = M[(r * N + col) * 2], fi0 = M[(r * N + col) * 2 + 1];
        // factor = M[r,col] / M[col,col]
        const fr = (fr0 * pr + fi0 * pi) / pden;
        const fi = (fi0 * pr - fr0 * pi) / pden;
        if (fr === 0 && fi === 0) continue;
        for (let k = col; k < N; k++) {
          const ar = M[(col * N + k) * 2], ai = M[(col * N + k) * 2 + 1];
          // M[r,k] -= factor * M[col,k]
          M[(r * N + k) * 2]     -= fr * ar - fi * ai;
          M[(r * N + k) * 2 + 1] -= fr * ai + fi * ar;
        }
        const br = x[col * 2], bi = x[col * 2 + 1];
        x[r * 2]     -= fr * br - fi * bi;
        x[r * 2 + 1] -= fr * bi + fi * br;
      }
    }
    // 대각 나누기
    for (let i = 0; i < N; i++) {
      const dr = M[(i * N + i) * 2], di = M[(i * N + i) * 2 + 1];
      const den = dr * dr + di * di;
      const xr = x[i * 2], xi = x[i * 2 + 1];
      x[i * 2]     = (xr * dr + xi * di) / den;
      x[i * 2 + 1] = (xi * dr - xr * di) / den;
    }
    return x; // 길이 N*2
  }

  // =====================================================================
  // 3. 물리 계산: 도선 전류 + 격자 위 복소장 (입사/산란, A=1 기준)
  // =====================================================================
  const solver = {
    k: 0, aEff_m: 0, wiresY: [],   // 도선 y좌표 [m] (x=0)
    cRe: null, cIm: null,          // 도선별 산란 계수 (A=1)
    gridW: 0, gridH: 0, Xw: 0, Yw: 0,
    incRe: null, incIm: null, scRe: null, scIm: null, // 격자 복소장 (A=1)
    tau: 1,
  };

  function recompute() {
    const lam_m = state.lam_cm / 100;
    const d_m = state.d_mm / 1000;
    const aEff_m = Math.min(state.a_mm, A_RATIO_MAX * state.d_mm) / 1000;
    const k = TWO_PI / lam_m;
    const N = state.N;

    // 도선 위치 (y중심 0)
    const wiresY = new Float64Array(N);
    for (let n = 0; n < N; n++) wiresY[n] = (n - (N - 1) / 2) * d_m;

    // --- MoM 행렬: Z_mn = H0(k*ρ_mn),  Z_mm = H0(k*aEff),  b_m = -1 (입사=+1)
    const Z = new Float64Array(N * N * 2);
    const b = new Float64Array(N * 2);
    const Hself = hankel0(k * aEff_m);
    for (let m = 0; m < N; m++) {
      b[m * 2] = -1; b[m * 2 + 1] = 0;
      for (let n = 0; n < N; n++) {
        let h;
        if (m === n) h = Hself;
        else h = hankel0(k * Math.abs(wiresY[m] - wiresY[n]));
        Z[(m * N + n) * 2] = h.re;
        Z[(m * N + n) * 2 + 1] = h.im;
      }
    }
    const c = solveComplex(N, Z, b);
    const cRe = new Float64Array(N), cIm = new Float64Array(N);
    for (let n = 0; n < N; n++) { cRe[n] = c[n * 2]; cIm[n] = c[n * 2 + 1]; }

    // --- 화면 창(window): 배열 전체가 항상 보이도록 자동 확대
    //     Yw = max(3λ × 종횡비,  배열전체높이 × 1.25) 로 잡고,
    //     Xw = Yw / 종횡비 로 유지 → 원통파 원형 왜곡 없음
    const arrHalf = (N - 1) * d_m / 2;
    const aspect = layout.bandH / layout.bandW;
    const YwMin = 3 * lam_m * aspect;          // 파장 기준 최소 세로
    const YwArr = arrHalf * 1.25;              // 배열 전체 보기 기준
    const Yw = Math.max(YwMin, YwArr);
    const Xw = Yw / aspect;
    const gridW = layout.gridW;
    const gridH = Math.max(40, Math.round(gridW * aspect));

    const incRe = new Float32Array(gridW * gridH);
    const incIm = new Float32Array(gridW * gridH);
    const scRe = new Float32Array(gridW * gridH);
    const scIm = new Float32Array(gridW * gridH);
    const aEff2 = aEff_m;

    for (let gj = 0; gj < gridH; gj++) {
      const wy = Yw - (gj + 0.5) / gridH * 2 * Yw;
      for (let gi = 0; gi < gridW; gi++) {
        const wx = -Xw + (gi + 0.5) / gridW * 2 * Xw;
        const idx = gj * gridW + gi;
        // 입사파 Ψ = e^{+ikx}  (왼→오른 진행, 시간인자 e^{-iωt}와 결합하면 cos(kx-ωt))
        const ph = k * wx;
        incRe[idx] = Math.cos(ph);
        incIm[idx] = Math.sin(ph);
        // 산란파 Σ c_n H0(k ρ_n)
        let sr = 0, si = 0;
        for (let n = 0; n < N; n++) {
          const dy = wy - wiresY[n];
          let r = Math.sqrt(wx * wx + dy * dy);
          if (r < aEff2) r = aEff2;       // 도선 내부 특이점 클램프
          const x = k * r;
          const jr = besselJ0(x), yi = besselY0(x);
          // c_n * (jr + i yi)
          sr += cRe[n] * jr - cIm[n] * yi;
          si += cRe[n] * yi + cIm[n] * jr;
        }
        scRe[idx] = sr; scIm[idx] = si;
      }
    }

    Object.assign(solver, {
      k, aEff_m, wiresY, cRe, cIm, gridW, gridH, Xw, Yw,
      incRe, incIm, scRe, scIm,
    });

    // E⊥wire 산란 약화 계수 τ (정성적): 얇은 도선은 거의 투명, ka 클수록 약간 산란
    const ka = k * aEff_m;
    solver.tau = state.polParallel ? 1 : Math.min(0.35, 2.0 * ka * ka);

    computeTransmittance();
    updateInfo();
  }

  // 전력 투과율 T = <|E_total|^2> / |E_inc|^2  (왼쪽 투과영역, A=1 기준이라 A 무관)
  //  측정점: 배열 바로 뒤 고정 거리(파장에 무관) · 배열 중앙부.
  //  (멀리서 재면 유한배열 가장자리 회절이 섞여 차폐경향이 흐트러지므로 바로 뒤에서 잰다.)
  let transmittance = 0;
  function computeTransmittance() {
    const k = solver.k, N = state.N, wiresY = solver.wiresY;
    const d_m = state.d_mm / 1000;
    const arrHalf = (N - 1) * d_m / 2;
    const xMeas = +Math.max(0.02, 1.5 * d_m);           // 배열 오른쪽(투과영역) 고정 거리 [m]
    const yHalf = Math.min(0.3 * arrHalf, 0.06);       // 중앙부만 (가장자리 회피)
    const samples = 41;
    let sum = 0;
    for (let s = 0; s < samples; s++) {
      const wy = -yHalf + (s / (samples - 1)) * 2 * yHalf;
      // 입사
      let tr = Math.cos(k * xMeas), ti = -Math.sin(k * xMeas);
      // 산란 (τ 반영)
      let sr = 0, si = 0;
      for (let n = 0; n < N; n++) {
        const dy = wy - wiresY[n];
        let r = Math.sqrt(xMeas * xMeas + dy * dy);
        if (r < solver.aEff_m) r = solver.aEff_m;
        const x = k * r;
        const jr = besselJ0(x), yi = besselY0(x);
        sr += solver.cRe[n] * jr - solver.cIm[n] * yi;
        si += solver.cRe[n] * yi + solver.cIm[n] * jr;
      }
      tr += solver.tau * sr; ti += solver.tau * si;
      sum += tr * tr + ti * ti;
    }
    transmittance = sum / samples; // |E_inc|^2 = 1
  }

  // =====================================================================
  // 4. 레이아웃 / 캔버스
  // =====================================================================
  const canvas = document.getElementById("canvas");
  const ctx = canvas.getContext("2d");
  const offscreen = document.createElement("canvas");
  const offctx = offscreen.getContext("2d");

  const layout = {
    cssW: 0, cssH: 0, marginL: 12, marginR: 12, marginT: 10, marginB: 10,
    gap: 14, bandX: 0, bandW: 0, bandH: 0, bandY: [0, 0, 0],
    gridW: 360,
  };

  const BAND_TITLES = ["① 입사파", "② 산란파", "③ 중첩 (입사 + 산란)"];

  function resize() {
    const rect = canvas.parentElement.getBoundingClientRect();
    layout.cssW = rect.width; layout.cssH = rect.height;
    const dpr = window.devicePixelRatio || 1;
    canvas.width = Math.round(layout.cssW * dpr);
    canvas.height = Math.round(layout.cssH * dpr);
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);

    layout.bandX = layout.marginL;
    layout.bandW = layout.cssW - layout.marginL - layout.marginR;
    const totalH = layout.cssH - layout.marginT - layout.marginB - 2 * layout.gap;
    layout.bandH = totalH / 3;
    for (let i = 0; i < 3; i++)
      layout.bandY[i] = layout.marginT + i * (layout.bandH + layout.gap);

    offscreen.width = layout.gridW;
    offscreen.height = Math.max(40, Math.round(layout.gridW * layout.bandH / layout.bandW));

    recompute();
    drawFrame();
  }

  // =====================================================================
  // 5. 색 매핑 + 프레임 렌더
  // =====================================================================
  function colorFor(v, out, o) {
    // v: 정규화된 [-1,1] 범위 기대값 (클램프)
    let t = v; if (t > 1) t = 1; else if (t < -1) t = -1;
    let r, g, bl;
    if (t >= 0) { r = 255; g = 255 - t * 205; bl = 255 - t * 215; }   // 흰→빨강
    else { const u = -t; r = 255 - u * 215; g = 255 - u * 165; bl = 255 - u * 35; } // 흰→파랑
    out[o] = r; out[o + 1] = g; out[o + 2] = bl; out[o + 3] = 255;
  }

  function drawFrame() {
    ctx.clearRect(0, 0, layout.cssW, layout.cssH);
    const gw = solver.gridW, gh = solver.gridH;
    const A = state.amp, tau = solver.tau;
    const cosP = Math.cos(state.phase), sinP = Math.sin(state.phase);

    if (offscreen.width !== gw || offscreen.height !== gh) {
      offscreen.width = gw; offscreen.height = gh;
    }
    const img = offctx.createImageData(gw, gh);
    const data = img.data;

    for (let band = 0; band < 3; band++) {
      // 밴드별 복소장 선택
      for (let p = 0; p < gw * gh; p++) {
        let fr, fi;
        if (band === 0) { fr = solver.incRe[p]; fi = solver.incIm[p]; }
        else if (band === 1) { fr = solver.scRe[p] * tau; fi = solver.scIm[p] * tau; }
        else { fr = solver.incRe[p] + solver.scRe[p] * tau; fi = solver.incIm[p] + solver.scIm[p] * tau; }
        const val = (fr * cosP + fi * sinP) * A;   // 실수 전기장 [V/m]
        colorFor(val / VMAX, data, p * 4);
      }
      offctx.putImageData(img, 0, 0);
      // 밴드 위치에 확대 그리기
      const by = layout.bandY[band];
      ctx.imageSmoothingEnabled = true;
      ctx.drawImage(offscreen, layout.bandX, by, layout.bandW, layout.bandH);
      drawOverlay(band, by);
    }
  }

  function worldToBand(wx, wy, by) {
    const sx = layout.bandX + (wx + solver.Xw) / (2 * solver.Xw) * layout.bandW;
    const sy = by + (solver.Yw - wy) / (2 * solver.Yw) * layout.bandH;
    return { x: sx, y: sy };
  }

  function drawOverlay(band, by) {
    const bx = layout.bandX, bw = layout.bandW, bh = layout.bandH;
    // 테두리
    ctx.strokeStyle = "#c8c8ce"; ctx.lineWidth = 1;
    ctx.strokeRect(bx + 0.5, by + 0.5, bw - 1, bh - 1);

    // 도선 (산란/중첩 칸에 표시; 입사 칸엔 위치만 옅게)
    const sPx = layout.bandW / (2 * solver.Xw);     // px per meter (가로=세로 동일)
    const dPx = state.d_mm / 1000 * sPx;            // 도선 간격 [px]
    // 반지름: 물리값(a)은 λ 대비 매우 작아 픽셀로 보이지 않으므로 a/d 비율 × 5배 과장 표시
    //         (물리 계산은 aEff_m으로 정확히 수행, 표시만 과장)
    const aRatio = solver.aEff_m / (state.d_mm / 1000);  // a/d
    const rPx = Math.min(dPx * 0.48, Math.max(2.5, dPx * aRatio * 5));
    const showWires = band !== 0;
    // 중심선
    const top = worldToBand(0, solver.Yw, by), bot = worldToBand(0, -solver.Yw, by);
    ctx.save();
    ctx.strokeStyle = band === 0 ? "#e4e4e8" : "#9aa0aa";
    ctx.setLineDash([3, 4]); ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(top.x, top.y); ctx.lineTo(bot.x, bot.y); ctx.stroke();
    ctx.restore();
    for (let n = 0; n < state.N; n++) {
      const p = worldToBand(0, solver.wiresY[n], by);
      if (p.y < by - 4 || p.y > by + bh + 4) continue;
      ctx.beginPath(); ctx.arc(p.x, p.y, rPx, 0, TWO_PI);
      if (showWires) { ctx.fillStyle = "#3a3a40"; ctx.fill(); ctx.lineWidth = 1; ctx.strokeStyle = "#1c1c1f"; ctx.stroke(); }
      else { ctx.fillStyle = "rgba(120,120,128,0.45)"; ctx.fill(); }
    }

    // 진행방향 화살표 (입사 칸, 왼→오른)
    if (band === 0) {
      const ay = by + 16, ax = bx + 16;
      ctx.fillStyle = "#444"; ctx.strokeStyle = "#444"; ctx.lineWidth = 1.5;
      ctx.beginPath(); ctx.moveTo(ax, ay); ctx.lineTo(ax + 34, ay); ctx.stroke();
      ctx.beginPath(); ctx.moveTo(ax + 34, ay); ctx.lineTo(ax + 27, ay - 4); ctx.lineTo(ax + 27, ay + 4); ctx.closePath(); ctx.fill();
      ctx.font = "11px sans-serif"; ctx.textAlign = "left";
      ctx.fillText("입사파 진행 →", ax, ay - 8);
      // 편광 표시
      drawPolIndicator(bx + 14, by + 36);
    }

    // 밴드 제목
    ctx.font = "bold 13px sans-serif"; ctx.textAlign = "left";
    ctx.fillStyle = "#10193a";
    ctx.fillText(BAND_TITLES[band], bx + 10, by + bh - 10);
    if (band === 2) {
      ctx.font = "11px sans-serif"; ctx.fillStyle = "#5a5a62";
      ctx.fillText("오른쪽=차폐영역 · 왼쪽=반사 간섭무늬(차폐 강할수록 정재파에 근접)", bx + 120, by + bh - 10);
    }
  }

  function drawPolIndicator(x, y) {
    ctx.save();
    ctx.textAlign = "left"; ctx.font = "11px sans-serif";
    if (state.polParallel) {
      // 도선과 평행 = 화면 안↔밖 ⊙
      ctx.strokeStyle = "#7a1fa0"; ctx.fillStyle = "#7a1fa0"; ctx.lineWidth = 1.5;
      ctx.beginPath(); ctx.arc(x + 6, y + 4, 6, 0, TWO_PI); ctx.stroke();
      ctx.beginPath(); ctx.arc(x + 6, y + 4, 1.6, 0, TWO_PI); ctx.fill();
      ctx.fillText("E ∥ 도선 (화면 안↔밖)", x + 18, y + 8);
    } else {
      // 도선과 수직 = 화면 위↔아래 ↕
      ctx.strokeStyle = "#1a8a4a"; ctx.fillStyle = "#1a8a4a"; ctx.lineWidth = 1.6;
      ctx.beginPath(); ctx.moveTo(x + 6, y - 3); ctx.lineTo(x + 6, y + 11); ctx.stroke();
      ctx.beginPath(); ctx.moveTo(x + 6, y - 4); ctx.lineTo(x + 3, y); ctx.lineTo(x + 9, y); ctx.closePath(); ctx.fill();
      ctx.beginPath(); ctx.moveTo(x + 6, y + 12); ctx.lineTo(x + 3, y + 8); ctx.lineTo(x + 9, y + 8); ctx.closePath(); ctx.fill();
      ctx.fillText("E ⊥ 도선 (화면 위↔아래)", x + 18, y + 8);
    }
    ctx.restore();
  }

  // =====================================================================
  // 6. 정보 표시
  // =====================================================================
  function updateInfo() {
    const lam_m = state.lam_cm / 100, d_m = state.d_mm / 1000;
    const dlam = d_m / lam_m;
    const f_GHz = C_LIGHT / lam_m / 1e9;
    const aEff = Math.min(state.a_mm, A_RATIO_MAX * state.d_mm);
    document.getElementById("infoBox").innerHTML =
      `파장 λ = <b>${state.lam_cm.toFixed(1)} cm</b> &nbsp;(f ≈ <b>${f_GHz.toFixed(2)} GHz</b>)<br>` +
      `간격 d = <b>${state.d_mm.toFixed(1)} mm</b> · 반지름 a = <b>${aEff.toFixed(2)} mm</b><br>` +
      `도선 수 N = <b>${state.N}</b> · 편광 <b>${state.polParallel ? "E∥wire" : "E⊥wire"}</b><br>` +
      `<b>d/λ = ${dlam.toFixed(3)}</b><br>` +
      `전력 투과율 <b>T = ${(transmittance * 100).toFixed(1)} %</b>`;
  }

  // =====================================================================
  // 7. UI 바인딩
  // =====================================================================
  function syncLabels() {
    const aEff = Math.min(state.a_mm, A_RATIO_MAX * state.d_mm);
    const aMax = (A_RATIO_MAX * state.d_mm);
    document.getElementById("aVal").textContent =
      state.a_mm.toFixed(2) + " mm" + (state.a_mm > aMax + 1e-9 ? " →" + aEff.toFixed(2) : "");
    document.getElementById("dVal").textContent = state.d_mm.toFixed(1) + " mm";
    document.getElementById("nVal").textContent = state.N + " 개";
    document.getElementById("lamVal").textContent = state.lam_cm.toFixed(1) + " cm";
    document.getElementById("ampVal").textContent = state.amp.toFixed(2) + " V/m";
    document.getElementById("aSlider").max = Math.max(0.05, aMax).toFixed(2);
    document.getElementById("polLegend").innerHTML = state.polParallel
      ? "전기장이 도선과 <b>평행</b> → 전류가 잘 유도되어 <b>차폐</b>."
      : "전기장이 도선과 <b>수직</b> → 도선이 거의 투명, <b>통과</b> (정성적 표현).";
  }

  let recomputeTimer = null;
  function scheduleRecompute() {
    if (recomputeTimer) clearTimeout(recomputeTimer);
    recomputeTimer = setTimeout(() => { recompute(); drawFrame(); }, 60);
  }

  function bindSlider(id, key, parse, heavy) {
    const el = document.getElementById(id);
    el.addEventListener("input", () => {
      state[key] = parse(el.value);
      syncLabels();
      if (heavy) scheduleRecompute();   // 격자 재계산 필요 (a,d,N,λ)
      else drawFrame();                 // A는 재계산 불필요
    });
  }
  bindSlider("aSlider", "a_mm", parseFloat, true);
  bindSlider("dSlider", "d_mm", parseFloat, true);
  bindSlider("nSlider", "N", v => parseInt(v, 10), true);
  bindSlider("lamSlider", "lam_cm", parseFloat, true);
  bindSlider("ampSlider", "amp", parseFloat, false);

  document.getElementById("polPar").addEventListener("click", () => setPol(true));
  document.getElementById("polPerp").addEventListener("click", () => setPol(false));
  function setPol(par) {
    state.polParallel = par;
    document.getElementById("polPar").classList.toggle("active", par);
    document.getElementById("polPerp").classList.toggle("active", !par);
    syncLabels();
    // τ 갱신 + 투과율 재계산 (격자 복소장은 그대로, τ만 변경)
    const ka = solver.k * solver.aEff_m;
    solver.tau = par ? 1 : Math.min(0.35, 2.0 * ka * ka);
    computeTransmittance(); updateInfo(); drawFrame();
  }

  const playBtn = document.getElementById("playBtn");
  const phaseWrap = document.getElementById("phaseWrap");
  const phaseSlider = document.getElementById("phaseSlider");
  playBtn.addEventListener("click", () => {
    state.playing = !state.playing;
    playBtn.textContent = state.playing ? "‖ 일시정지" : "▶ 재생";
    phaseWrap.classList.toggle("on", !state.playing);
    document.getElementById("phaseHint").style.display = state.playing ? "none" : "block";
  });
  phaseSlider.addEventListener("input", () => {
    state.phase = parseFloat(phaseSlider.value) * Math.PI / 180;
    document.getElementById("phaseVal").textContent = phaseSlider.value + "°";
    if (!state.playing) drawFrame();
  });

  // =====================================================================
  // 8. 애니메이션 루프
  // =====================================================================
  function loop() {
    if (state.playing) {
      state.phase = (state.phase + 0.06) % TWO_PI;   // 일정한 화면상 속도
      drawFrame();
    }
    requestAnimationFrame(loop);
  }

  // =====================================================================
  // 9. 콘솔 자가검증
  // =====================================================================
  function selfCheck() {
    // 한켈 근사 점검 (알려진 값과 비교)
    console.log("[검증] J0(1)=", besselJ0(1).toFixed(6), "(기대 0.765198)");
    console.log("[검증] Y0(1)=", besselY0(1).toFixed(6), "(기대 0.088257)");
    console.log("[검증] J0(5)=", besselJ0(5).toFixed(6), "(기대 -0.177597)");
    console.log("[검증] Y0(0.01)=", besselY0(0.01).toFixed(4), "(유한, 발산 아님)");

    // 차폐 경향: 같은 배열에서 파장만 바꿔 투과율 T 비교 (짧은 λ vs 긴 λ)
    const save = { ...state };
    function Tfor(lam_cm) {
      state.lam_cm = lam_cm; state.polParallel = true;
      recompute();
      return transmittance;
    }
    const tShort = Tfor(2);    // d/λ 큼 → 누설
    const tLong = Tfor(25);    // d/λ 작음 → 차폐
    console.log(`[검증] 차폐 경향: T(λ=2cm)=${(tShort * 100).toFixed(1)}%  vs  T(λ=25cm)=${(tLong * 100).toFixed(1)}%  → 긴 파장에서 작아야 함:`, tLong < tShort ? "OK" : "확인필요");
    Object.assign(state, save);
    recompute();
  }

  // =====================================================================
  // 시작
  // =====================================================================
  syncLabels();
  window.addEventListener("resize", resize);
  resize();
  selfCheck();
  // selfCheck가 state를 복원하므로 라벨/화면 재동기화
  syncLabels(); recompute(); drawFrame();
  requestAnimationFrame(loop);
})();
