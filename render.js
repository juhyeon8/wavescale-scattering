"use strict";
(function (global) {
  function clamp(v, lo, hi) { return v < lo ? lo : (v > hi ? hi : v); }

  function divergingColor(value, scale) {
    const t = clamp(value / scale, -1, 1);
    if (t >= 0) return [255, Math.round(255 * (1 - t)), Math.round(255 * (1 - t))];
    const a = -t;
    return [Math.round(255 * (1 - a)), Math.round(255 * (1 - a)), 255];
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
    const yTop = (Ny - 1 - geom.jTop) * z, yBot = (Ny - 1 - geom.jBot) * z;
    ctx.strokeStyle = "#444"; ctx.lineWidth = 3;
    ctx.beginPath();
    ctx.moveTo(geom.plateColStart * z, yTop); ctx.lineTo(geom.plateColEnd * z, yTop);
    ctx.moveTo(geom.plateColStart * z, yBot); ctx.lineTo(geom.plateColEnd * z, yBot);
    ctx.stroke();
    const sx = geom.sourceI * z, sy = (Ny - 1 - geom.sourceJ) * z;
    ctx.fillStyle = "#111"; ctx.beginPath();
    ctx.arc(sx, sy, 5, 0, Math.PI * 2); ctx.fill();
  }

  function drawGraph(ctx, row, geom, opts) {
    const W = ctx.canvas.width, H = ctx.canvas.height;
    ctx.clearRect(0, 0, W, H);
    let mx = 1e-9; for (let i = 0; i < row.length; i++) mx = Math.max(mx, row[i]);
    ctx.strokeStyle = opts && opts.evanescent ? "#c0392b" : "#2f6feb";
    ctx.lineWidth = 2; ctx.beginPath();
    for (let i = 0; i < row.length; i++) {
      const x = i / (row.length - 1) * W;
      const y = H - (row[i] / mx) * (H - 8) - 4;
      if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    }
    ctx.stroke();
    if (geom && geom.mouthI != null) {
      const xm = geom.mouthI / (row.length - 1) * W;
      ctx.strokeStyle = "#999"; ctx.setLineDash([4, 4]);
      ctx.beginPath(); ctx.moveTo(xm, 0); ctx.lineTo(xm, H); ctx.stroke();
      ctx.setLineDash([]);
    }
  }

  const API = { divergingColor, paintField, drawGuide, drawGraph };
  if (typeof module !== "undefined" && module.exports) module.exports = API;
  else { global.WaveSim = global.WaveSim || {}; global.WaveSim.render = API; }
})(typeof globalThis !== "undefined" ? globalThis : this);
