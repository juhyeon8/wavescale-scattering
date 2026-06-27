"use strict";
const assert = require("node:assert");
const { test, approx, done } = require("./_harness.js");
const P = require("../physics.js");

test("cutoffWavelength = 2a", () => {
  assert.strictEqual(P.cutoffWavelength(24), 48);
});
test("regimeOf: x안쪽+y판사이만 internal, 그 외 aperture", () => {
  const geom = { mouthI: 100, jBot: 20, jTop: 40 };
  assert.strictEqual(P.regimeOf(40, 30, geom), "aperture");   // x가 입구 왼쪽
  assert.strictEqual(P.regimeOf(120, 30, geom), "internal");  // x 안쪽 + y 판 사이
  assert.strictEqual(P.regimeOf(120, 10, geom), "aperture");  // x 안쪽이나 y가 판 밖
});
test("attenuationConstant: λ>2a 양수, λ=2a 임계=0(소멸), λ<2a null", () => {
  // a=10 → kc=π/10. λ=42 → k<kc → 소멸
  approx(P.attenuationConstant(10, 42), Math.sqrt((Math.PI/10)**2 - (2*Math.PI/42)**2), 1e-9);
  approx(P.attenuationConstant(10, 20), 0, 1e-9);             // λ=2a 임계 → κ=0(소멸)
  assert.strictEqual(P.attenuationConstant(10, 15), null);    // λ<2a → 전파
});
test("propagatingModeCount", () => {
  assert.strictEqual(P.propagatingModeCount(10, 42), 0); // λ>2a → 0
  assert.strictEqual(P.propagatingModeCount(10, 15), 1); // 2a/λ=1.33 → n=1
  assert.strictEqual(P.propagatingModeCount(10, 8), 2);  // 2a/λ=2.5 → n=1,2
});
test("swrToReflection", () => {
  approx(P.swrToReflection(1), 0, 1e-12);
  approx(P.swrToReflection(3), 0.5, 1e-12);
});
test("modeExcitationStrength: 중심선 n=1 배=1, 판 근처≈0", () => {
  approx(P.modeExcitationStrength(5, 10, 1), 1, 1e-9);     // y=a/2
  approx(P.modeExcitationStrength(0.01, 10, 1), 0, 1e-2);  // 판 근처
});
test("fitExponential: 합성 지수데이터 복원", () => {
  const trueAmp = 2.5, trueKappa = 0.3;
  const xs = [], ys = [];
  for (let x = 0; x <= 20; x++) { xs.push(x); ys.push(trueAmp * Math.exp(-trueKappa * x)); }
  const r = P.fitExponential(xs, ys);
  approx(r.kappa, trueKappa, 1e-6, "kappa");
  approx(r.amp, trueAmp, 1e-5, "amp");
});
done();
