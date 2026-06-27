"use strict";
(function (global) {
  const isNode = (typeof module !== "undefined" && module.exports);
  const FDTD2D = isNode ? require("./fdtd.js").FDTD2D : global.WaveSim.FDTD2D;
  const P = isNode ? require("./physics.js") : global.WaveSim.physics;

  function WaveSim(cfg) {
    this.cfg = cfg;
    this.Nx = cfg.Nx; this.Ny = cfg.Ny; this.aCells = cfg.aCells;
    this.mouthI = cfg.mouthI;
    this.lambdaCells = 14; this.amp = 1; this.rampSteps = 20;
    // 판 위치: 중심행 기준 위/아래로 a/2
    this.jMid = Math.round(cfg.Ny / 2);
    this.jBot = this.jMid - Math.round(cfg.aCells / 2);
    this.jTop = this.jMid + Math.round(cfg.aCells / 2);
    this.sourceI = Math.round(cfg.Nx / 2); this.sourceJ = this.jMid;

    this.incident = new FDTD2D(cfg.Nx, cfg.Ny, cfg.courant);
    this.total = new FDTD2D(cfg.Nx, cfg.Ny, cfg.courant);
    const pml = cfg.pml || 10;
    this.incident.setupPML(pml); this.total.setupPML(pml); // 두 격자 동일 PML
    // 전체 격자에만 PEC 평행판 마스크
    const mask = new Uint8Array(cfg.Nx * cfg.Ny);
    if (cfg.plateColEnd >= cfg.plateColStart) {
      for (let i = cfg.plateColStart; i <= cfg.plateColEnd; i++) {
        mask[this.total.idx(i, this.jBot)] = 1;
        mask[this.total.idx(i, this.jTop)] = 1;
      }
    }
    this.total.setPECMask(mask);
    this._tmp = new Float32Array(cfg.Nx * cfg.Ny);
  }

  WaveSim.prototype.setSourceCell = function (i, j) {
    i = Math.max(1, Math.min(this.Nx - 2, Math.round(i)));
    // 판 사이 구간(도파관 안)일 때만 j 클램프; 자유공간 바깥은 자유
    if (i >= this.cfg.plateColStart && i <= this.cfg.plateColEnd &&
        this.cfg.plateColEnd >= this.cfg.plateColStart) {
      j = Math.max(this.jBot + 2, Math.min(this.jTop - 2, Math.round(j)));
    } else {
      j = Math.max(1, Math.min(this.Ny - 2, Math.round(j)));
    }
    this.sourceI = i; this.sourceJ = j;
    const src = [{ i: i, j: j, weight: 1 }];
    this.incident.setSources(src);
    this.total.setSources(src.map(function (s) { return { i: s.i, j: s.j, weight: s.weight }; }));
    return { i: i, j: j };
  };
  WaveSim.prototype.setLambda = function (l) { this.lambdaCells = l; };
  WaveSim.prototype.setAmp = function (a) { this.amp = a; };

  WaveSim.prototype.step = function () {
    const p = { amp: this.amp, lambdaCells: this.lambdaCells, rampSteps: this.rampSteps };
    this.incident.step(p); this.total.step(p);
  };
  // 테스트 전용 배치 측정: 런타임 애니메이션 루프에서는 호출 금지
  WaveSim.prototype.measurePhasors = function (periods) {
    periods = periods || 1;
    const period = Math.round(this.lambdaCells / this.incident.S);
    this.incident.beginPhasor(this.lambdaCells);
    this.total.beginPhasor(this.lambdaCells);
    for (let n = 0; n < periods * period; n++) {
      this.step(); this.incident.accumulate(); this.total.accumulate();
    }
  };
  // 런타임용: 애니메이션 루프가 매 스텝 accumulateMeasure()로 누적
  WaveSim.prototype.periodSteps = function () {
    return Math.round(this.lambdaCells / this.incident.S);
  };
  WaveSim.prototype.beginMeasure = function () {
    this.incident.beginPhasor(this.lambdaCells);
    this.total.beginPhasor(this.lambdaCells);
  };
  WaveSim.prototype.accumulateMeasure = function () {
    this.incident.accumulate(); this.total.accumulate();
  };

  WaveSim.prototype.incidentAmp = function (out) { return this.incident.phasorAmp(out); };
  WaveSim.prototype.totalAmp = function (out) { return this.total.phasorAmp(out); };
  WaveSim.prototype.scatteredRe = function (out) {
    var ci = this.incident.phasorComplex(), ct = this.total.phasorComplex();
    for (var k = 0; k < out.length; k++) {
      out[k] = ct.factor * ct.re[k] - ci.factor * ci.re[k];
    }
    return out;
  };

  WaveSim.prototype.scatteredAmp = function (out) {
    const ci = this.incident.phasorComplex(), ct = this.total.phasorComplex();
    for (let k = 0; k < out.length; k++) {
      const re = ct.factor * ct.re[k] - ci.factor * ci.re[k];
      const im = ct.factor * ct.im[k] - ci.factor * ci.im[k];
      out[k] = Math.hypot(re, im);
    }
    return out;
  };
  WaveSim.prototype.incidentFrame = function () { return this.incident.Ez; };
  WaveSim.prototype.totalFrame = function () { return this.total.Ez; };
  WaveSim.prototype.scatteredFrame = function () {
    const t = this.total.Ez, inc = this.incident.Ez, out = this._tmp;
    for (let k = 0; k < out.length; k++) out[k] = t[k] - inc[k];
    return out;
  };
  WaveSim.prototype.centerlineAmp = function () {
    this.total.phasorAmp(this._tmp);
    const row = new Float32Array(this.Nx);
    for (let i = 0; i < this.Nx; i++) row[i] = this._tmp[this.total.idx(i, this.jMid)];
    return row;
  };
  WaveSim.prototype.regime = function () {
    return P.regimeOf(this.sourceI, this.sourceJ,
      { mouthI: this.mouthI, jBot: this.jBot, jTop: this.jTop });
  };
  WaveSim.prototype.reflectionPercent = function (centerline) {
    var row = centerline || this.centerlineAmp();
    var mn = Infinity, mx = 0;
    for (var i = this.mouthI + 5; i < this.Nx - 5; i++) {
      if (row[i] <= 0) continue;
      if (row[i] < mn) mn = row[i];
      if (row[i] > mx) mx = row[i];
    }
    if (!isFinite(mn) || mn <= 0) return null;
    var swr = mx / mn;
    return P.swrToReflection(swr) * 100;
  };

  WaveSim.prototype.cutoffInfo = function () {
    const kappa = P.attenuationConstant(this.aCells, this.lambdaCells);
    return {
      lambdaCells: this.lambdaCells,
      lambdaCCells: P.cutoffWavelength(this.aCells),
      modeCount: P.propagatingModeCount(this.aCells, this.lambdaCells),
      evanescent: kappa !== null,
      kappa: kappa,
    };
  };

  const API = { WaveSim };
  if (isNode) module.exports = API;
  else {
    global.WaveSim = global.WaveSim || {};
    global.WaveSim.WaveSim = WaveSim;
    global.WaveSim.Sim = WaveSim; // main.js에서 WaveSim.Sim으로 참조
  }
})(typeof globalThis !== "undefined" ? globalThis : this);
