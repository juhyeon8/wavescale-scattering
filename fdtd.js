"use strict";
(function (global) {
  function FDTD2D(Nx, Ny, courant) {
    this.Nx = Nx; this.Ny = Ny;
    this.S = (courant === undefined) ? 0.5 : courant;
    const sz = Nx * Ny;
    this.Ezx = new Float32Array(sz);    // 분리장 x성분
    this.Ezy = new Float32Array(sz);    // 분리장 y성분
    this.Ez  = new Float32Array(sz);    // = Ezx + Ezy
    this.Hx = new Float32Array(sz);
    this.Hy = new Float32Array(sz);
    this.mask = null;
    this.sources = [];
    this.n = 0;
    this.accRe = new Float32Array(sz);
    this.accIm = new Float32Array(sz);
    this.accCount = 0;
    this.phasorLambda = 0;
    // 셀별 업데이트 계수(기본=무손실)
    this.caEx = new Float32Array(sz); this.cbEx = new Float32Array(sz); // Ezx (x손실)
    this.caEy = new Float32Array(sz); this.cbEy = new Float32Array(sz); // Ezy (y손실)
    this.daHx = new Float32Array(sz); this.dbHx = new Float32Array(sz); // Hx (y손실)
    this.daHy = new Float32Array(sz); this.dbHy = new Float32Array(sz); // Hy (x손실)
    this._setLossless();
  }
  FDTD2D.prototype.idx = function (i, j) { return i * this.Ny + j; };

  FDTD2D.prototype._setLossless = function () {
    const S = this.S, sz = this.Nx * this.Ny;
    for (let k = 0; k < sz; k++) {
      this.caEx[k] = 1; this.cbEx[k] = S;
      this.caEy[k] = 1; this.cbEy[k] = S;
      this.daHx[k] = 1; this.dbHx[k] = S;
      this.daHy[k] = 1; this.dbHy[k] = S;
    }
  };

  // 네 경계에 두께 D의 분리장(Berenger) PML. PEC로 종단(외곽 링은 갱신 안 함).
  FDTD2D.prototype.setupPML = function (D, sigmaMax, m) {
    D = D || 10; sigmaMax = (sigmaMax === undefined) ? 0.8 : sigmaMax; m = m || 3;
    this._setLossless();
    const Nx = this.Nx, Ny = this.Ny, S = this.S;
    function eCoef(sig) { return { a: (1 - sig / 2) / (1 + sig / 2), b: S / (1 + sig / 2) }; }
    // x방향 PML(좌·우) → Ezx, Hy
    for (let i = 0; i < Nx; i++) {
      let depth = 0;
      if (i < D) depth = (D - i) / D;
      else if (i >= Nx - D) depth = (i - (Nx - D - 1)) / D;
      if (depth <= 0) continue;
      const c = eCoef(sigmaMax * Math.pow(depth, m));
      for (let j = 0; j < Ny; j++) {
        const k = i * Ny + j;
        this.caEx[k] = c.a; this.cbEx[k] = c.b;
        this.daHy[k] = c.a; this.dbHy[k] = c.b;
      }
    }
    // y방향 PML(상·하) → Ezy, Hx
    for (let j = 0; j < Ny; j++) {
      let depth = 0;
      if (j < D) depth = (D - j) / D;
      else if (j >= Ny - D) depth = (j - (Ny - D - 1)) / D;
      if (depth <= 0) continue;
      const c = eCoef(sigmaMax * Math.pow(depth, m));
      for (let i = 0; i < Nx; i++) {
        const k = i * Ny + j;
        this.caEy[k] = c.a; this.cbEy[k] = c.b;
        this.daHx[k] = c.a; this.dbHx[k] = c.b;
      }
    }
  };

  FDTD2D.prototype.reset = function () {
    this.Ezx.fill(0); this.Ezy.fill(0); this.Ez.fill(0);
    this.Hx.fill(0); this.Hy.fill(0);
    this.accRe.fill(0); this.accIm.fill(0); this.accCount = 0; this.n = 0;
  };
  FDTD2D.prototype.setPECMask = function (mask) { this.mask = mask; };
  FDTD2D.prototype.setSources = function (list) { this.sources = list || []; };

  FDTD2D.prototype.step = function (p) {
    const Nx = this.Nx, Ny = this.Ny;
    const Ez = this.Ez, Ezx = this.Ezx, Ezy = this.Ezy, Hx = this.Hx, Hy = this.Hy;
    // --- H 업데이트 ---
    for (let i = 0; i < Nx; i++) {
      for (let j = 0; j < Ny - 1; j++) {
        const k = i * Ny + j;
        Hx[k] = this.daHx[k] * Hx[k] - this.dbHx[k] * (Ez[k + 1] - Ez[k]); // ∂Ez/∂y
      }
    }
    for (let i = 0; i < Nx - 1; i++) {
      for (let j = 0; j < Ny; j++) {
        const k = i * Ny + j;
        Hy[k] = this.daHy[k] * Hy[k] + this.dbHy[k] * (Ez[k + Ny] - Ez[k]); // ∂Ez/∂x
      }
    }
    // --- E 업데이트(외곽 링 제외 → PEC 종단) ---
    for (let i = 1; i < Nx - 1; i++) {
      for (let j = 1; j < Ny - 1; j++) {
        const k = i * Ny + j;
        Ezx[k] = this.caEx[k] * Ezx[k] + this.cbEx[k] * (Hy[k] - Hy[k - Ny]);
        Ezy[k] = this.caEy[k] * Ezy[k] - this.cbEy[k] * (Hx[k] - Hx[k - 1]);
        Ez[k] = Ezx[k] + Ezy[k];
      }
    }
    // --- 소스(소프트, Ezx/Ezy 절반씩) ---
    const ramp = Math.min(1, this.n / Math.max(1, p.rampSteps));
    const phase = 2 * Math.PI * this.S * this.n / p.lambdaCells;
    const drive = p.amp * ramp * Math.sin(phase);
    for (let s = 0; s < this.sources.length; s++) {
      const src = this.sources[s], k = src.i * Ny + src.j;
      Ezx[k] += 0.5 * drive * src.weight;
      Ezy[k] += 0.5 * drive * src.weight;
      Ez[k] = Ezx[k] + Ezy[k];
    }
    // --- PEC(도체판) ---
    if (this.mask) {
      const m = this.mask;
      for (let k = 0; k < m.length; k++) if (m[k]) { Ezx[k] = 0; Ezy[k] = 0; Ez[k] = 0; }
    }
    this.n++;
  };

  FDTD2D.prototype.beginPhasor = function (lambdaCells) {
    this.accRe.fill(0); this.accIm.fill(0); this.accCount = 0;
    this.phasorLambda = lambdaCells;
  };
  FDTD2D.prototype.accumulate = function () {
    const ph = 2 * Math.PI * this.S * this.n / this.phasorLambda;
    const c = Math.cos(ph), s = Math.sin(ph), Ez = this.Ez;
    const re = this.accRe, im = this.accIm;
    for (let k = 0; k < Ez.length; k++) { re[k] += Ez[k] * c; im[k] += Ez[k] * s; }
    this.accCount++;
  };
  FDTD2D.prototype.phasorComplex = function () {
    const f = this.accCount ? 2 / this.accCount : 0;
    return { re: this.accRe, im: this.accIm, factor: f };
  };
  FDTD2D.prototype.phasorAmp = function (out) {
    const f = this.accCount ? 2 / this.accCount : 0;
    const re = this.accRe, im = this.accIm;
    for (let k = 0; k < out.length; k++) out[k] = f * Math.hypot(re[k], im[k]);
    return out;
  };

  const API = { FDTD2D };
  if (typeof module !== "undefined" && module.exports) module.exports = API;
  else { global.WaveSim = global.WaveSim || {}; global.WaveSim.FDTD2D = FDTD2D; }
})(typeof globalThis !== "undefined" ? globalThis : this);
