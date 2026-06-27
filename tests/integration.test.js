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
done();
