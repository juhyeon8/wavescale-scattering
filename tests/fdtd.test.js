"use strict";
const assert = require("node:assert");
const { test, approx, done } = require("./_harness.js");
const { FDTD2D } = require("../fdtd.js");

test("초기화: 필드 0, 인덱스 규칙", () => {
  const g = new FDTD2D(8, 6);
  assert.strictEqual(g.Ez.length, 48);
  assert.strictEqual(g.idx(2, 3), 2 * 6 + 3);
  assert.strictEqual(g.Ez[0], 0);
});

test("점 소스 → 좌우/상하 대칭장(자유공간)", () => {
  const N = 41, g = new FDTD2D(N, N);
  const c = (N - 1) / 2;
  g.setSources([{ i: c, j: c, weight: 1 }]);
  for (let n = 0; n < 20; n++) g.step({ amp: 1, lambdaCells: 12, rampSteps: 5 });
  // 중심 기준 대칭: (c+5,c) 와 (c-5,c) 거의 동일
  approx(g.Ez[g.idx(c + 5, c)], g.Ez[g.idx(c - 5, c)], 1e-4, "좌우대칭");
  approx(g.Ez[g.idx(c, c + 5)], g.Ez[g.idx(c, c - 5)], 1e-4, "상하대칭");
});

test("안정성: 200스텝 후 NaN/발산 없음", () => {
  const N = 60, g = new FDTD2D(N, N);
  g.setupPML(10, 0.8, 3); // 흡수경계로 공진 누적 없이 유계
  g.setSources([{ i: 30, j: 30, weight: 1 }]);
  for (let n = 0; n < 200; n++) g.step({ amp: 1, lambdaCells: 16, rampSteps: 10 });
  let maxAbs = 0;
  for (let k = 0; k < g.Ez.length; k++) {
    assert.ok(Number.isFinite(g.Ez[k]), "유한값");
    maxAbs = Math.max(maxAbs, Math.abs(g.Ez[k]));
  }
  assert.ok(maxAbs < 100, "발산하지 않음: " + maxAbs);
});

test("PEC 마스크 셀은 Ez=0 유지", () => {
  const N = 40, g = new FDTD2D(N, N);
  const mask = new Uint8Array(N * N);
  mask[g.idx(20, 10)] = 1;
  g.setPECMask(mask);
  g.setSources([{ i: 20, j: 20, weight: 1 }]);
  for (let n = 0; n < 30; n++) g.step({ amp: 1, lambdaCells: 12, rampSteps: 5 });
  assert.strictEqual(g.Ez[g.idx(20, 10)], 0);
});

test("위상자 누적: 정상상태 진폭 > 0", () => {
  const N = 50, g = new FDTD2D(N, N);
  g.setSources([{ i: 25, j: 25, weight: 1 }]);
  const lambda = 14;
  for (let n = 0; n < 60; n++) g.step({ amp: 1, lambdaCells: lambda, rampSteps: 10 });
  g.beginPhasor(lambda);
  const period = Math.round(lambda / 0.5); // T_steps = λ/S
  for (let n = 0; n < period; n++) { g.step({ amp: 1, lambdaCells: lambda, rampSteps: 10 }); g.accumulate(); }
  const amp = new Float32Array(N * N);
  g.phasorAmp(amp);
  assert.ok(amp[g.idx(27, 25)] > 0, "소스 근처 진폭 양수");
});

test("PML: 경계 반사 감소 — PML이 PEC 박스보다 후기 잔향 작음", () => {
  // 2D 점원은 하이겐스 후행파(wake)가 남아 절대 잔향은 0이 안 된다.
  // 따라서 '절대 잔향'이 아니라, 반사하는 PEC 박스 대비 흡수하는 PML의
  // 후기 잔향 감소로 흡수를 검증한다(공통 wake는 상쇄됨).
  function run(usePML) {
    const N = 80, g = new FDTD2D(N, N);
    if (usePML) g.setupPML(10, 0.8, 3);
    const c = 40; g.setSources([{ i: c, j: c, weight: 1 }]);
    const probe = g.idx(c + 20, c); // 경계 쪽 탐침
    let early = 0, late = 0;
    for (let n = 0; n < 250; n++) {
      const amp = n < 40 ? 1 : 0;     // 짧은 톤버스트 후 소스 끔
      g.step({ amp: amp, lambdaCells: 12, rampSteps: 5 });
      const v = Math.abs(g.Ez[probe]);
      if (n < 120) early = Math.max(early, v);
      else late = Math.max(late, v);
    }
    return { early: early, late: late };
  }
  const pec = run(false), pml = run(true);
  assert.ok(pec.early > 1e-3 && pml.early > 1e-3, "초기 통과 신호 존재");
  assert.ok(pml.late < 0.5 * pec.late,
    "PML 후기 잔향이 PEC 박스보다 작음: pml=" + pml.late + " pec=" + pec.late);
});
done();
