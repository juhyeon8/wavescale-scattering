"use strict";
const assert = require("node:assert");
const { test, done } = require("./_harness.js");
const U = require("../ui.js");

test("panel1Label: 영역별 라벨", () => {
  assert.strictEqual(U.panel1Label("aperture"), "① 입사파 (원통파)");
  assert.strictEqual(U.panel1Label("internal"), "① 자유공간 기준장 (판 없을 때)");
});
test("regimeBadge", () => {
  assert.ok(U.regimeBadge("aperture").includes("자유공간"));
  assert.ok(U.regimeBadge("internal").includes("도파관 내부"));
});
test("cutoffBadge: 항상 '도파관 내부' 한정 문구 포함", () => {
  const info = { modeCount: 0, evanescent: true };
  assert.ok(U.cutoffBadge(info, "aperture").includes("도파관 내부"));
  assert.ok(U.cutoffBadge({ modeCount: 1, evanescent: false }, "internal").includes("전파"));
});
done();
