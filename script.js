(function () {
  "use strict";

  // =====================================================================
  // 0. 상수
  // =====================================================================
  const C_LIGHT = 2.99792458e8;
  const TWO_PI = Math.PI * 2;
  const VMAX = 1.5;           // 색 포화 기준 [V/m]
  const A_RATIO_MAX = 0.30;   // 반지름 상한 = 간격의 30%
  const N_MAX_INF = 80;       // 무한 배열 탭 N 상한
  const D_MIN_INF = 3;        // 무한 배열 탭 d 최솟값 [mm] (HTML min 속성과 일치)

  // =====================================================================
  // 1. 상태
  //    shared  : 두 탭 공유 (λ, A, 편광, 애니메이션)
  //    tabState: 탭별 독립 (d, a, N)
  //    activeTab: 0=무한 배열(기본), 1=유한 배열
  // =====================================================================
  const shared = {
    lam_cm: 12.2, amp: 1.0, polParallel: true,
    playing: true, phase: 0,
  };
  const tabState = [
    { d_mm: 10, a_mm: 0.5, N: 0 },   // Tab 0: 무한 배열 (N 자동)
    { d_mm: 10, a_mm: 0.5, N: 30 },  // Tab 1: 유한 배열
  ];
  let activeTab = 0;

  // =====================================================================
  // 2. 베셀/한켈 함수  (Abramowitz & Stegun 9.4 다항 근사)
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
  function hankel0(x) { return { re: besselJ0(x), im: besselY0(x) }; }

  // =====================================================================
  // 3. 복소 선형계 Z c = b  (N×N, 부분 피벗 가우스 소거)
  // =====================================================================
  function solveComplex(N, Z, b) {
    const M = Z.slice();
    const x = b.slice();
    for (let col = 0; col < N; col++) {
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
      const pr = M[(col * N + col) * 2], pi = M[(col * N + col) * 2 + 1];
      const pden = pr * pr + pi * pi;
      for (let r = 0; r < N; r++) {
        if (r === col) continue;
        const fr0 = M[(r * N + col) * 2], fi0 = M[(r * N + col) * 2 + 1];
        const fr = (fr0 * pr + fi0 * pi) / pden;
        const fi = (fi0 * pr - fr0 * pi) / pden;
        if (fr === 0 && fi === 0) continue;
        for (let k = col; k < N; k++) {
          const ar = M[(col * N + k) * 2], ai = M[(col * N + k) * 2 + 1];
          M[(r * N + k) * 2]     -= fr * ar - fi * ai;
          M[(r * N + k) * 2 + 1] -= fr * ai + fi * ar;
        }
        const br = x[col * 2], bi = x[col * 2 + 1];
        x[r * 2]     -= fr * br - fi * bi;
        x[r * 2 + 1] -= fr * bi + fi * br;
      }
    }
    for (let i = 0; i < N; i++) {
      const dr = M[(i * N + i) * 2], di = M[(i * N + i) * 2 + 1];
      const den = dr * dr + di * di;
      const xr = x[i * 2], xi = x[i * 2 + 1];
      x[i * 2]     = (xr * dr + xi * di) / den;
      x[i * 2 + 1] = (xi * dr - xr * di) / den;
    }
    return x;
  }

  // =====================================================================
  // 4. 솔버 출력 (recompute 결과 저장)
  // =====================================================================
  const solver = {
    k: 0, aEff_m: 0, wiresY: [],
    cRe: null, cIm: null,
    gridW: 0, gridH: 0, Xw: 0, Yw: 0,
    incRe: null, incIm: null, scRe: null, scIm: null,
    tau: 1,
  };

  // =====================================================================
  // 5. 물리 계산
  // =====================================================================
  function recompute() {
    if (!layout.bandW || !layout.bandH) return;

    const ts = tabState[activeTab];
    const lam_m = shared.lam_cm / 100;
    const d_m = ts.d_mm / 1000;
    const aEff_m = Math.min(ts.a_mm, A_RATIO_MAX * ts.d_mm) / 1000;
    const k = TWO_PI / lam_m;
    const aspect = layout.bandH / layout.bandW;

    let N, Xw, Yw;

    if (activeTab === 0) {
      // 무한 배열: 창을 파장 기준으로 고정, N을 창 크기에 맞게 자동 결정
      Yw = 3 * lam_m * aspect;
      Xw = 3 * lam_m;
      N = Math.max(4, Math.min(N_MAX_INF, Math.floor(2 * Yw / d_m)));
      ts.N = N;
      const el = document.getElementById("autoNVal");
      if (el) el.textContent = N;
    } else {
      // 유한 배열: 배열 전체가 보이도록 창 자동 확대
      N = ts.N;
      const arrHalf = (N - 1) * d_m / 2;
      const YwMin = 3 * lam_m * aspect;
      const YwArr = arrHalf * 1.25;
      Yw = Math.max(YwMin, YwArr);
      Xw = Yw / aspect;
    }

    // 도선 위치 (y=0 중심)
    const wiresY = new Float64Array(N);
    for (let n = 0; n < N; n++) wiresY[n] = (n - (N - 1) / 2) * d_m;

    // MoM 행렬: Z_mn=H0(k·ρ_mn), Z_mm=H0(k·a), b_m=-1
    const Z = new Float64Array(N * N * 2);
    const b = new Float64Array(N * 2);
    const Hself = hankel0(k * aEff_m);
    for (let m = 0; m < N; m++) {
      b[m * 2] = -1; b[m * 2 + 1] = 0;
      for (let n = 0; n < N; n++) {
        const h = (m === n) ? Hself : hankel0(k * Math.abs(wiresY[m] - wiresY[n]));
        Z[(m * N + n) * 2] = h.re;
        Z[(m * N + n) * 2 + 1] = h.im;
      }
    }
    const c = solveComplex(N, Z, b);
    const cRe = new Float64Array(N), cIm = new Float64Array(N);
    for (let n = 0; n < N; n++) { cRe[n] = c[n * 2]; cIm[n] = c[n * 2 + 1]; }

    // 격자 복소장 사전 계산 (A=1 기준)
    const gridW = layout.gridW;
    const gridH = Math.max(40, Math.round(gridW * aspect));
    const incRe = new Float32Array(gridW * gridH);
    const incIm = new Float32Array(gridW * gridH);
    const scRe  = new Float32Array(gridW * gridH);
    const scIm  = new Float32Array(gridW * gridH);

    for (let gj = 0; gj < gridH; gj++) {
      const wy = Yw - (gj + 0.5) / gridH * 2 * Yw;
      for (let gi = 0; gi < gridW; gi++) {
        const wx = -Xw + (gi + 0.5) / gridW * 2 * Xw;
        const idx = gj * gridW + gi;
        // 입사파 e^{ikx}  (왼→오른: Re[e^{ikx}·e^{-iωt}] = cos(kx−ωt))
        const ph = k * wx;
        incRe[idx] = Math.cos(ph);
        incIm[idx] = Math.sin(ph);
        // 산란파 Σ c_n H0(k ρ_n)
        let sr = 0, si = 0;
        for (let n = 0; n < N; n++) {
          const dy = wy - wiresY[n];
          let r = Math.sqrt(wx * wx + dy * dy);
          if (r < aEff_m) r = aEff_m;
          const x = k * r;
          const jr = besselJ0(x), yi = besselY0(x);
          sr += cRe[n] * jr - cIm[n] * yi;
          si += cRe[n] * yi + cIm[n] * jr;
        }
        scRe[idx] = sr; scIm[idx] = si;
      }
    }

    Object.assign(solver, { k, aEff_m, wiresY, cRe, cIm, gridW, gridH, Xw, Yw, incRe, incIm, scRe, scIm });

    // E⊥wire 약화 계수 τ (정성적)
    const ka = k * aEff_m;
    solver.tau = shared.polParallel ? 1 : Math.min(0.35, 2.0 * ka * ka);

    computeTransmittance();
    updateInfo();
  }

  // 전력 투과율 T = <|E_total|²> / |E_inc|²
  //   측정점: 배열 오른쪽(투과영역) 고정 거리 · 중앙부 샘플
  let transmittance = 0;
  function computeTransmittance() {
    const ts = tabState[activeTab];
    const k = solver.k;
    const wiresY = solver.wiresY;
    const N = wiresY.length;
    const d_m = ts.d_mm / 1000;
    const xMeas = +Math.max(0.02, 1.5 * d_m);         // 배열 오른쪽 고정 거리
    const yHalf = Math.min(0.25 * solver.Yw, 0.06);   // 중앙부 샘플 범위
    const samples = 41;
    let sum = 0;
    for (let s = 0; s < samples; s++) {
      const wy = -yHalf + (s / (samples - 1)) * 2 * yHalf;
      // 입사파 복소 진폭 at xMeas: e^{ikx}
      let tr = Math.cos(k * xMeas), ti = Math.sin(k * xMeas);
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
    transmittance = sum / samples;
  }

  // =====================================================================
  // 6. 레이아웃 / 캔버스
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
  // 7. 색 매핑 + 프레임 렌더
  // =====================================================================
  function colorFor(v, out, o) {
    let t = v; if (t > 1) t = 1; else if (t < -1) t = -1;
    let r, g, bl;
    if (t >= 0) { r = 255; g = 255 - t * 205; bl = 255 - t * 215; }
    else { const u = -t; r = 255 - u * 215; g = 255 - u * 165; bl = 255 - u * 35; }
    out[o] = r; out[o + 1] = g; out[o + 2] = bl; out[o + 3] = 255;
  }

  function drawFrame() {
    ctx.clearRect(0, 0, layout.cssW, layout.cssH);
    const gw = solver.gridW, gh = solver.gridH;
    if (!gw || !gh) return;
    const A = shared.amp, tau = solver.tau;
    const cosP = Math.cos(shared.phase), sinP = Math.sin(shared.phase);

    if (offscreen.width !== gw || offscreen.height !== gh) {
      offscreen.width = gw; offscreen.height = gh;
    }
    const img = offctx.createImageData(gw, gh);
    const data = img.data;

    for (let band = 0; band < 3; band++) {
      for (let p = 0; p < gw * gh; p++) {
        let fr, fi;
        if (band === 0) { fr = solver.incRe[p]; fi = solver.incIm[p]; }
        else if (band === 1) { fr = solver.scRe[p] * tau; fi = solver.scIm[p] * tau; }
        else { fr = solver.incRe[p] + solver.scRe[p] * tau; fi = solver.incIm[p] + solver.scIm[p] * tau; }
        const val = (fr * cosP + fi * sinP) * A;
        colorFor(val / VMAX, data, p * 4);
      }
      offctx.putImageData(img, 0, 0);
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
    ctx.strokeStyle = "#c8c8ce"; ctx.lineWidth = 1;
    ctx.strokeRect(bx + 0.5, by + 0.5, bw - 1, bh - 1);

    const ts = tabState[activeTab];
    const sPx = layout.bandW / (2 * solver.Xw);
    const dPx = ts.d_mm / 1000 * sPx;
    const aRatio = solver.aEff_m / (ts.d_mm / 1000);
    // 반지름: a/d 비율×5배 과장 표시 (물리 계산은 aEff_m으로 정확히 수행)
    const rPx = Math.min(dPx * 0.48, Math.max(2.5, dPx * aRatio * 5));

    // 중심선 (도선 배열 위치)
    const top = worldToBand(0, solver.Yw, by), bot = worldToBand(0, -solver.Yw, by);
    ctx.save();
    ctx.strokeStyle = band === 0 ? "#e4e4e8" : "#9aa0aa";
    ctx.setLineDash([3, 4]); ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(top.x, top.y); ctx.lineTo(bot.x, bot.y); ctx.stroke();
    ctx.restore();

    // 도선
    const N = ts.N;
    for (let n = 0; n < N; n++) {
      const p = worldToBand(0, solver.wiresY[n], by);
      if (p.y < by - 4 || p.y > by + bh + 4) continue;
      ctx.beginPath(); ctx.arc(p.x, p.y, rPx, 0, TWO_PI);
      if (band !== 0) {
        ctx.fillStyle = "#3a3a40"; ctx.fill();
        ctx.lineWidth = 1; ctx.strokeStyle = "#1c1c1f"; ctx.stroke();
      } else {
        ctx.fillStyle = "rgba(120,120,128,0.45)"; ctx.fill();
      }
    }

    // 진행방향 화살표 (입사 칸, 왼→오른)
    if (band === 0) {
      const ay = by + 16, ax = bx + 16;
      ctx.fillStyle = "#444"; ctx.strokeStyle = "#444"; ctx.lineWidth = 1.5;
      ctx.beginPath(); ctx.moveTo(ax, ay); ctx.lineTo(ax + 34, ay); ctx.stroke();
      ctx.beginPath(); ctx.moveTo(ax + 34, ay); ctx.lineTo(ax + 27, ay - 4); ctx.lineTo(ax + 27, ay + 4); ctx.closePath(); ctx.fill();
      ctx.font = "11px sans-serif"; ctx.textAlign = "left"; ctx.fillStyle = "#444";
      ctx.fillText("입사파 진행 →", ax, ay - 8);
      drawPolIndicator(bx + 14, by + 36);
    }

    // 밴드 제목
    ctx.font = "bold 13px sans-serif"; ctx.textAlign = "left";
    ctx.fillStyle = "#10193a";
    ctx.fillText(BAND_TITLES[band], bx + 10, by + bh - 10);

    if (band === 2) {
      ctx.font = "11px sans-serif"; ctx.fillStyle = "#5a5a62";
      const caption = activeTab === 0
        ? "오른쪽=차폐영역 · 왼쪽=반사 간섭무늬 · 도선이 창에 가득 참(무한 배열 근사)"
        : "오른쪽=차폐영역 · 왼쪽=반사 간섭무늬(차폐 강할수록 정재파에 근접)";
      ctx.fillText(caption, bx + 120, by + bh - 10);
    }
  }

  function drawPolIndicator(x, y) {
    ctx.save();
    ctx.textAlign = "left"; ctx.font = "11px sans-serif";
    if (shared.polParallel) {
      ctx.strokeStyle = "#7a1fa0"; ctx.fillStyle = "#7a1fa0"; ctx.lineWidth = 1.5;
      ctx.beginPath(); ctx.arc(x + 6, y + 4, 6, 0, TWO_PI); ctx.stroke();
      ctx.beginPath(); ctx.arc(x + 6, y + 4, 1.6, 0, TWO_PI); ctx.fill();
      ctx.fillText("E ∥ 도선 (화면 안↔밖)", x + 18, y + 8);
    } else {
      ctx.strokeStyle = "#1a8a4a"; ctx.fillStyle = "#1a8a4a"; ctx.lineWidth = 1.6;
      ctx.beginPath(); ctx.moveTo(x + 6, y - 3); ctx.lineTo(x + 6, y + 11); ctx.stroke();
      ctx.beginPath(); ctx.moveTo(x + 6, y - 4); ctx.lineTo(x + 3, y); ctx.lineTo(x + 9, y); ctx.closePath(); ctx.fill();
      ctx.beginPath(); ctx.moveTo(x + 6, y + 12); ctx.lineTo(x + 3, y + 8); ctx.lineTo(x + 9, y + 8); ctx.closePath(); ctx.fill();
      ctx.fillText("E ⊥ 도선 (화면 위↔아래)", x + 18, y + 8);
    }
    ctx.restore();
  }

  // =====================================================================
  // 8. 정보 표시
  // =====================================================================
  function updateInfo() {
    const ts = tabState[activeTab];
    const lam_m = shared.lam_cm / 100;
    const d_m = ts.d_mm / 1000;
    const dlam = d_m / lam_m;
    const f_GHz = C_LIGHT / lam_m / 1e9;
    const aEff = Math.min(ts.a_mm, A_RATIO_MAX * ts.d_mm);
    const modeStr = activeTab === 0 ? "무한 배열" : "유한 배열";
    document.getElementById("infoBox").innerHTML =
      `[${modeStr}] &nbsp;N = <b>${ts.N}${activeTab === 0 ? " (자동)" : ""}</b><br>` +
      `파장 λ = <b>${shared.lam_cm.toFixed(1)} cm</b> &nbsp;(f ≈ <b>${f_GHz.toFixed(2)} GHz</b>)<br>` +
      `간격 d = <b>${ts.d_mm.toFixed(1)} mm</b> · 반지름 a = <b>${aEff.toFixed(2)} mm</b><br>` +
      `편광 <b>${shared.polParallel ? "E∥wire" : "E⊥wire"}</b> &nbsp;· &nbsp;<b>d/λ = ${dlam.toFixed(3)}</b><br>` +
      `전력 투과율 <b>T = ${(transmittance * 100).toFixed(1)} %</b>`;
  }

  function syncLabels() {
    // Tab 0 슬라이더 라벨
    const ts0 = tabState[0];
    const aMax0 = A_RATIO_MAX * ts0.d_mm;
    const aEff0 = Math.min(ts0.a_mm, aMax0);
    document.getElementById("a0Val").textContent =
      ts0.a_mm.toFixed(2) + " mm" + (ts0.a_mm > aMax0 + 1e-9 ? " →" + aEff0.toFixed(2) : "");
    document.getElementById("d0Val").textContent = ts0.d_mm.toFixed(1) + " mm";
    document.getElementById("a0Slider").max = Math.max(0.05, aMax0).toFixed(2);

    // Tab 1 슬라이더 라벨
    const ts1 = tabState[1];
    const aMax1 = A_RATIO_MAX * ts1.d_mm;
    const aEff1 = Math.min(ts1.a_mm, aMax1);
    document.getElementById("a1Val").textContent =
      ts1.a_mm.toFixed(2) + " mm" + (ts1.a_mm > aMax1 + 1e-9 ? " →" + aEff1.toFixed(2) : "");
    document.getElementById("d1Val").textContent = ts1.d_mm.toFixed(1) + " mm";
    document.getElementById("n1Val").textContent = ts1.N + " 개";
    document.getElementById("a1Slider").max = Math.max(0.05, aMax1).toFixed(2);

    // 공유 슬라이더 라벨
    document.getElementById("lamVal").textContent = shared.lam_cm.toFixed(1) + " cm";
    document.getElementById("ampVal").textContent = shared.amp.toFixed(2) + " V/m";

    // 편광 범례
    document.getElementById("polLegend").innerHTML = shared.polParallel
      ? "전기장이 도선과 <b>평행</b> → 전류가 잘 유도되어 <b>차폐</b>."
      : "전기장이 도선과 <b>수직</b> → 도선이 거의 투명, <b>통과</b> (정성적 표현).";
  }

  // =====================================================================
  // 9. UI 바인딩
  // =====================================================================
  let recomputeTimer = null;
  function scheduleRecompute() {
    if (recomputeTimer) clearTimeout(recomputeTimer);
    recomputeTimer = setTimeout(() => { recompute(); drawFrame(); }, 60);
  }

  // 탭별 슬라이더: 해당 탭 상태 업데이트, 현재 탭일 때만 재계산
  function bindTabSlider(id, tabIdx, key, parse) {
    document.getElementById(id).addEventListener("input", function () {
      tabState[tabIdx][key] = parse(this.value);
      syncLabels();
      if (activeTab === tabIdx) scheduleRecompute();
    });
  }
  bindTabSlider("a0Slider", 0, "a_mm", parseFloat);
  bindTabSlider("d0Slider", 0, "d_mm", parseFloat);
  bindTabSlider("a1Slider", 1, "a_mm", parseFloat);
  bindTabSlider("d1Slider", 1, "d_mm", parseFloat);
  bindTabSlider("n1Slider", 1, "N", v => parseInt(v, 10));

  // 공유 슬라이더
  document.getElementById("lamSlider").addEventListener("input", function () {
    shared.lam_cm = parseFloat(this.value);
    syncLabels(); scheduleRecompute();
  });
  document.getElementById("ampSlider").addEventListener("input", function () {
    shared.amp = parseFloat(this.value);
    syncLabels(); drawFrame();
  });

  // 탭 전환
  document.querySelectorAll(".tabBtn").forEach(btn => {
    btn.addEventListener("click", function () {
      const newTab = parseInt(this.dataset.tab);
      if (newTab === activeTab) return;
      activeTab = newTab;
      document.querySelectorAll(".tabBtn").forEach((b, i) => b.classList.toggle("active", i === newTab));
      document.querySelectorAll(".tabPane").forEach((p, i) => p.classList.toggle("active", i === newTab));
      recompute(); drawFrame();
    });
  });

  // 편광 토글
  document.getElementById("polPar").addEventListener("click", () => setPol(true));
  document.getElementById("polPerp").addEventListener("click", () => setPol(false));
  function setPol(par) {
    shared.polParallel = par;
    document.getElementById("polPar").classList.toggle("active", par);
    document.getElementById("polPerp").classList.toggle("active", !par);
    syncLabels();
    const ka = solver.k * solver.aEff_m;
    solver.tau = par ? 1 : Math.min(0.35, 2.0 * ka * ka);
    computeTransmittance(); updateInfo(); drawFrame();
  }

  // 재생/일시정지
  const playBtn = document.getElementById("playBtn");
  const phaseWrap = document.getElementById("phaseWrap");
  const phaseSlider = document.getElementById("phaseSlider");
  playBtn.addEventListener("click", () => {
    shared.playing = !shared.playing;
    playBtn.textContent = shared.playing ? "‖ 일시정지" : "▶ 재생";
    phaseWrap.classList.toggle("on", !shared.playing);
    document.getElementById("phaseHint").style.display = shared.playing ? "none" : "block";
  });
  phaseSlider.addEventListener("input", () => {
    shared.phase = parseFloat(phaseSlider.value) * Math.PI / 180;
    document.getElementById("phaseVal").textContent = phaseSlider.value + "°";
    if (!shared.playing) drawFrame();
  });

  // =====================================================================
  // 10. 애니메이션 루프
  // =====================================================================
  function loop() {
    if (shared.playing) {
      shared.phase = (shared.phase + 0.06) % TWO_PI;
      drawFrame();
    }
    requestAnimationFrame(loop);
  }

  // =====================================================================
  // 11. 콘솔 자가검증
  // =====================================================================
  function selfCheck() {
    console.log("[검증] J0(1)=", besselJ0(1).toFixed(6), "(기대 0.765198)");
    console.log("[검증] Y0(1)=", besselY0(1).toFixed(6), "(기대 0.088257)");
    console.log("[검증] J0(5)=", besselJ0(5).toFixed(6), "(기대 -0.177597)");
    console.log("[검증] Y0(0.01)=", besselY0(0.01).toFixed(4), "(유한, 발산 아님)");

    const savedLam = shared.lam_cm, savedPol = shared.polParallel, savedTab = activeTab;
    activeTab = 0; shared.polParallel = true;
    shared.lam_cm = 2;  recompute(); const tShort = transmittance;
    shared.lam_cm = 25; recompute(); const tLong = transmittance;
    console.log(
      `[검증] T(λ=2cm)=${(tShort*100).toFixed(1)}%  T(λ=25cm)=${(tLong*100).toFixed(1)}%` +
      `  → 긴 파장이 더 작아야:`, tLong < tShort ? "OK" : "확인필요"
    );
    shared.lam_cm = savedLam; shared.polParallel = savedPol; activeTab = savedTab;
  }

  // =====================================================================
  // 시작
  // =====================================================================
  syncLabels();
  window.addEventListener("resize", resize);
  resize();       // layout 확정 + recompute + drawFrame
  selfCheck();    // 베셀·차폐 경향 검증 (내부에서 recompute 호출)
  recompute();    // selfCheck 후 원래 상태로 복원
  drawFrame();
  requestAnimationFrame(loop);
})();
