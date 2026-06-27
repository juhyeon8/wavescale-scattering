"use strict";
const assert = require("node:assert");
const { test, done } = require("./_harness.js");
const R = require("../render.js");

test("divergingColor: 0=백, 양수=적쪽, 음수=청쪽", () => {
  assert.deepStrictEqual(R.divergingColor(0, 1), [255, 255, 255]);
  const pos = R.divergingColor(1, 1);   // 최대 양수
  assert.ok(pos[0] > pos[2], "양수는 적색 우세");
  const neg = R.divergingColor(-1, 1);  // 최대 음수
  assert.ok(neg[2] > neg[0], "음수는 청색 우세");
});
test("divergingColor: 스케일 밖 값은 클램프", () => {
  assert.deepStrictEqual(R.divergingColor(5, 1), R.divergingColor(1, 1));
});
done();
