"use strict";
(function (global) {
  function clamp(v, lo, hi) { return v < lo ? lo : (v > hi ? hi : v); }

  var GAMMA = 0.45;
  function divergingColor(value, scale) {
    var t = clamp(value / scale, -1, 1);
    var mag = Math.pow(Math.abs(t), GAMMA);
    if (t >= 0) return [255, Math.round(255 * (1 - mag)), Math.round(255 * (1 - mag))];
    return [Math.round(255 * (1 - mag)), Math.round(255 * (1 - mag)), 255];
  }

  // field: Float32Array(Nx*Ny), i=x열 j=y행. 화면 y는 위가 0이므로 j를 뒤집음.
  function paintField(ctx, field, Nx, Ny, scale, geom) {
    const img = ctx.createImageData(Nx, Ny);
    const d = img.data;
    for (let j = 0; j < Ny; j++) {
      for (let i = 0; i < Nx; i++) {
        const v = field[i * Ny + j];
        const c = divergingColor(v, scale);
        const px = ((Ny - 1 - j) * Nx + i) * 4;
        d[px] = c[0]; d[px + 1] = c[1]; d[px + 2] = c[2]; d[px + 3] = 255;
      }
    }
    ctx.imageSmoothingEnabled = false;
    ctx.putImageData(img, 0, 0);
    if (geom && geom.zoom && geom.zoom !== 1) {
      ctx.drawImage(ctx.canvas, 0, 0, Nx, Ny, 0, 0, geom.zoom * Nx, geom.zoom * Ny);
    }
  }

  function drawGuide(ctx, geom) {
    const z = geom.zoom || 1, Ny = geom.Ny;
    const wt = geom.wallThick || 1; // 벽 두께(셀 수)
    const x0 = geom.plateColStart * z;
    const xw = (geom.plateColEnd - geom.plateColStart) * z;
    // 화면 y: screen_y = (Ny-1-j)*z (j=0이 화면 아래)
    // 상단 벽: j = jTop(내면) … jTop+wt-1(외면) → screen_y 작은 쪽부터 아래로
    const yTopInner = (Ny - 1 - geom.jTop) * z;
    const yTopOuter = (Ny - 1 - (geom.jTop + wt - 1)) * z; // yTopOuter < yTopInner
    // 하단 벽: j = jBot(내면) … jBot-wt+1(외면)
    const yBotInner = (Ny - 1 - geom.jBot) * z;
    const yBotOuter = (Ny - 1 - (geom.jBot - wt + 1)) * z; // yBotOuter > yBotInner
    ctx.fillStyle = "rgba(40,40,40,0.85)";
    ctx.fillRect(x0, yTopOuter, xw, yTopInner - yTopOuter + 1); // 상단 벽
    ctx.fillRect(x0, yBotInner, xw, yBotOuter - yBotInner + 1); // 하단 벽
    const sx = geom.sourceI * z, sy = (Ny - 1 - geom.sourceJ) * z;
    ctx.fillStyle = "#ff6600"; ctx.beginPath();
    ctx.arc(sx, sy, 4, 0, Math.PI * 2); ctx.fill();
  }

  function drawGraph(ctx, row, geom, opts) {
    const W = ctx.canvas.width, H = ctx.canvas.height;
    ctx.clearRect(0, 0, W, H);
    const pml = (geom && geom.pml) ? geom.pml : 0;
    const nx  = (geom && geom.Nx)  ? geom.Nx  : row.length;
    const i0 = pml, i1 = nx - 1 - pml;
    let mx = 1e-9;
    for (let i = i0; i <= i1; i++) mx = Math.max(mx, row[i]);
    ctx.strokeStyle = opts && opts.evanescent ? "#c0392b" : "#2f6feb";
    ctx.lineWidth = 2; ctx.beginPath();
    for (let i = i0; i <= i1; i++) {
      const x = (i - i0) / (i1 - i0) * W;
      const y = H - (row[i] / mx) * (H - 8) - 4;
      if (i === i0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    }
    ctx.stroke();
    if (geom && geom.mouthI != null) {
      const xm = (geom.mouthI - i0) / (i1 - i0) * W;
      ctx.strokeStyle = "#999"; ctx.setLineDash([4, 4]);
      ctx.beginPath(); ctx.moveTo(xm, 0); ctx.lineTo(xm, H); ctx.stroke();
      ctx.setLineDash([]);
    }
  }

  function drawMarkers(ctx, markerA, markerB, geom) {
    var z = geom.zoom || 1, Ny = geom.Ny;
    var sy = (Ny - 1 - geom.jMid) * z;
    ctx.lineWidth = 2;
    ctx.strokeStyle = "#e67e22";
    ctx.beginPath(); ctx.arc(markerA * z, sy, 6, 0, Math.PI * 2); ctx.stroke();
    ctx.strokeStyle = "#2980b9";
    ctx.beginPath(); ctx.arc(markerB * z, sy, 6, 0, Math.PI * 2); ctx.stroke();
  }

  function drawTimeGraph(ctx, ezA, ezB, head, count, period) {
    var W = ctx.canvas.width, H = ctx.canvas.height;
    ctx.clearRect(0, 0, W, H);
    var BUF = ezA.length;
    var show = Math.min(count, Math.round(3 * period), BUF);
    if (show < 2) return;
    // 자동 스케일: 표시 구간 내 최대값 기준
    var mx = 1e-9;
    for (var m = 0; m < show; m++) {
      var mi = (head - show + m + BUF) % BUF;
      var av = Math.abs(ezA[mi]); if (av > mx) mx = av;
      var bv = Math.abs(ezB[mi]); if (bv > mx) mx = bv;
    }
    ctx.strokeStyle = "#e0e0e0"; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(0, H / 2); ctx.lineTo(W, H / 2); ctx.stroke();
    ctx.strokeStyle = "#e67e22"; ctx.lineWidth = 1.5; ctx.beginPath();
    for (var k = 0; k < show; k++) {
      var idx = (head - show + k + BUF) % BUF;
      var x = k / (show - 1) * W;
      var y = H / 2 - (ezA[idx] / mx) * (H / 2 - 3);
      if (k === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    }
    ctx.stroke();
    ctx.strokeStyle = "#2980b9"; ctx.lineWidth = 1.5; ctx.beginPath();
    for (var k2 = 0; k2 < show; k2++) {
      var idx2 = (head - show + k2 + BUF) % BUF;
      var x2 = k2 / (show - 1) * W;
      var y2 = H / 2 - (ezB[idx2] / mx) * (H / 2 - 3);
      if (k2 === 0) ctx.moveTo(x2, y2); else ctx.lineTo(x2, y2);
    }
    ctx.stroke();
  }

  const API = { divergingColor, paintField, drawGuide, drawGraph, drawMarkers, drawTimeGraph };
  if (typeof module !== "undefined" && module.exports) module.exports = API;
  else { global.WaveSim = global.WaveSim || {}; global.WaveSim.render = API; }
})(typeof globalThis !== "undefined" ? globalThis : this);
