"use strict";
const assert = require("node:assert");
const { test, approx, done } = require("./_harness.js");
const { WaveSim } = require("../sim.js");
const P = require("../physics.js");

test("내부여기 λ>2a: 중심선 |E| 지수감쇠가 이론 κ와 일치", () => {
  const Nx = 200, Ny = 60, a = 20;
  const cfg = { Nx, Ny, aCells: a, mouthI: 0,
                plateColStart: 0, plateColEnd: Nx - 1, courant: 0.5 };
  const s = new WaveSim(cfg);
  const lambda = 56; // λ>2a=40 → 소멸
  s.setLambda(lambda);
  // 순수 n=1 모드 소스를 i=30 열 전체(판 사이)에 sin(πy/a) 가중으로 주입
  const i0 = 30, src = [];
  for (let j = s.jBot + 1; j < s.jTop; j++) {
    const y = j - s.jBot;
    src.push({ i: i0, j: j, weight: Math.sin(Math.PI * y / a) });
  }
  s.total.setSources(src); s.incident.setSources([]);
  // 정상상태까지(소멸파는 수렴에 충분한 스텝 필요)
  for (let n = 0; n < 2000; n++) s.step();
  s.total.beginPhasor(lambda);
  const period = Math.round(lambda / 0.5);
  for (let n = 0; n < 3 * period; n++) { s.step(); s.total.accumulate(); }
  const amp = new Float32Array(Nx * Ny);
  s.total.phasorAmp(amp);
  // 소스에서 충분히 떨어진 구간(근접장 제외)에서 중심선 샘플
  const xs = [], ys = [], jMid = s.jMid;
  for (let i = i0 + 10; i <= i0 + 40; i++) { xs.push(i - i0); ys.push(amp[s.total.idx(i, jMid)]); }
  const fit = P.fitExponential(xs, ys);
  const kappaTheory = P.attenuationConstant(a, lambda);
  // FDTD 수치분산 감안, 상대오차 20% 이내
  approx(fit.kappa / kappaTheory, 1, 0.2, "κ 상대오차");
});

test("내부여기 λ<2a: 전파(감쇠 거의 없음)", () => {
  const Nx = 200, Ny = 60, a = 20;
  const cfg = { Nx, Ny, aCells: a, mouthI: 0,
                plateColStart: 0, plateColEnd: Nx - 1, courant: 0.5 };
  const s = new WaveSim(cfg);
  const lambda = 28; // λ<2a=40 → 전파
  s.setLambda(lambda);
  const i0 = 30, src = [];
  for (let j = s.jBot + 1; j < s.jTop; j++) {
    const y = j - s.jBot;
    src.push({ i: i0, j: j, weight: Math.sin(Math.PI * y / a) });
  }
  s.total.setSources(src); s.incident.setSources([]);
  for (let n = 0; n < 2000; n++) s.step();
  s.total.beginPhasor(lambda);
  const period = Math.round(lambda / 0.5);
  for (let n = 0; n < 3 * period; n++) { s.step(); s.total.accumulate(); }
  const amp = new Float32Array(Nx * Ny);
  s.total.phasorAmp(amp);
  const a1 = amp[s.total.idx(i0 + 20, s.jMid)];
  const a2 = amp[s.total.idx(i0 + 50, s.jMid)];
  assert.ok(a2 > 0.4 * a1, "전파면 멀리서도 진폭 유지: " + a1 + " -> " + a2);
});

// ── 마커 A-B 위상차: total 위상자 + 소스 같은 쪽 마커 ───────────────
function phaseTestSetup(lambda, a) {
  const Nx = 200, Ny = 60;
  const cfg = { Nx, Ny, aCells: a, mouthI: 0,
                plateColStart: 0, plateColEnd: Nx - 1, courant: 0.5 };
  const s = new WaveSim(cfg);
  s.setLambda(lambda);
  const i0 = 60, src = [];
  for (let j = s.jBot + 1; j < s.jTop; j++) {
    src.push({ i: i0, j: j, weight: Math.sin(Math.PI * (j - s.jBot) / a) });
  }
  s.total.setSources(src); s.incident.setSources([]);
  for (let n = 0; n < 3000; n++) s.step();
  s.total.beginPhasor(lambda);
  const period = Math.round(lambda / 0.5);
  for (let n = 0; n < 4 * period; n++) { s.step(); s.total.accumulate(); }
  const pc = s.total.phasorComplex();
  // 소스 오른쪽 같은 쪽 마커
  const mA = i0 + Math.max(10, Math.round(0.5 * a));
  const mB = i0 + Math.max(20, Math.round(1.2 * a));
  const kA = mA * Ny + s.jMid, kB = mB * Ny + s.jMid;
  const reA = pc.factor * pc.re[kA], imA = pc.factor * pc.im[kA];
  const reB = pc.factor * pc.re[kB], imB = pc.factor * pc.im[kB];
  let dph = (Math.atan2(imB, reB) - Math.atan2(imA, reA)) * 180 / Math.PI;
  while (dph > 180) dph -= 360;
  while (dph < -180) dph += 360;
  return dph;
}

test("total 위상자 A-B: 내부 차단(λ>2a) → 위상차 ≤ 20°(제자리 진동)", () => {
  // λ=56, a=20: κ≈0.110/cell, mA=i0+10, mB=i0+24
  // FDTD 수치분산으로 ~1°/셀 위상 축적 → 14셀 간격에서 최대 ~16° 허용
  const dph = phaseTestSetup(56, 20);
  assert.ok(Math.abs(dph) < 20, `차단 위상차 ${dph.toFixed(1)}° ≥ 20°`);
});

test("total 위상자 A-B: 내부 전파(λ<2a) → 위상차 > 15°(진행파)", () => {
  // λ=28, a=20: k_guide≈0.161/cell, spacing=14 → dph≈129°
  const dph = phaseTestSetup(28, 20);
  assert.ok(Math.abs(dph) > 15, `전파 위상차 ${dph.toFixed(1)}° ≤ 15°`);
});

done();
