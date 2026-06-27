"use strict";
const assert = require("node:assert");
const { test, approx, done } = require("./_harness.js");
const { WaveSim } = require("../sim.js");

function baseCfg() {
  return { Nx: 120, Ny: 60, aCells: 24, mouthI: 40,
           plateColStart: 40, plateColEnd: 119, courant: 0.5 };
}

test("소스 클램프: 판과 최소 2셀", () => {
  const s = new WaveSim(baseCfg());
  // 판은 중심행 ±a/2. 판 바로 위에 두려 해도 2셀 안으로 못 들어감
  const jBot = s.jBot, jTop = s.jTop;
  const c = s.setSourceCell(80, jBot + 1);
  assert.ok(c.j >= jBot + 2 && c.j <= jTop - 2, "클램프됨: " + c.j);
});

test("두 격자 소스 인덱스 동일", () => {
  const s = new WaveSim(baseCfg());
  s.setSourceCell(80, 30);
  assert.deepStrictEqual(s.incident.sources, s.total.sources);
});

test("판이 없으면 산란≈0", () => {
  const cfg = baseCfg(); cfg.plateColStart = 0; cfg.plateColEnd = -1; // 판 없음
  const s = new WaveSim(cfg);
  s.setSourceCell(60, 30); s.setLambda(14);
  for (let n = 0; n < 80; n++) s.step();
  const out = new Float32Array(cfg.Nx * cfg.Ny);
  s.measurePhasors(2);
  s.scatteredAmp(out);
  let mx = 0; for (let k = 0; k < out.length; k++) mx = Math.max(mx, out[k]);
  assert.ok(mx < 1e-3, "산란 거의 0: " + mx);
});

test("scatteredRe: 판 없으면 ≈ 0", () => {
  const cfg = baseCfg(); cfg.plateColStart = 0; cfg.plateColEnd = -1;
  const s = new WaveSim(cfg);
  s.setSourceCell(60, 30); s.setLambda(14);
  for (let n = 0; n < 80; n++) s.step();
  s.measurePhasors(2);
  const out = new Float32Array(cfg.Nx * cfg.Ny);
  s.scatteredRe(out);
  let mx = 0; for (let k = 0; k < out.length; k++) mx = Math.max(mx, Math.abs(out[k]));
  assert.ok(mx < 1e-3, "산란 Re 거의 0: " + mx);
});

test("산란 순간장 = 전체 − 입사", () => {
  const s = new WaveSim(baseCfg());
  s.setSourceCell(80, 30); s.setLambda(14);
  for (let n = 0; n < 30; n++) s.step();
  const sc = s.scatteredFrame();
  const k = s.total.idx(50, 30);
  approx(sc[k], s.total.Ez[k] - s.incident.Ez[k], 1e-6);
});

test("cutoffInfo: λ>2a 면 소멸·모드0", () => {
  const s = new WaveSim(baseCfg());
  s.setLambda(60); // λ=60 > 2a=48
  const ci = s.cutoffInfo();
  assert.strictEqual(ci.modeCount, 0);
  assert.strictEqual(ci.evanescent, true);
});

test("reflectionPercent: 정재파비에서 % 산출", () => {
  const s = new WaveSim(baseCfg());
  // 합성 중심선: 내부 구간에 정재파(min/max) 주입해 SWR 검증
  const row = new Float32Array(s.Nx);
  for (let i = 0; i < s.Nx; i++) {
    row[i] = (i > s.mouthI) ? (2 + Math.cos(i * 0.5)) : 0; // max3 min1 → SWR3
  }
  const pct = s.reflectionPercent(row);
  // SWR=3 → |R|=0.5 → 50%
  approx(pct, 50, 1, "반사율%");
});
done();
