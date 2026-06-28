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

// ── computeMarkers ──────────────────────────────────────────────────
const PLATE = { plateColStart: 40, plateColEnd: 159, mouthI: 40 };

test("computeMarkers 내부: 소스 중앙 → 오른쪽 같은 쪽 배치", () => {
  const { markerA, markerB } = U.computeMarkers("internal", 100, 24,
    PLATE.plateColStart, PLATE.plateColEnd, PLATE.mouthI);
  // markerA = 100+max(10,12)=112, markerB = 100+max(20,29)=129
  assert.ok(markerA > 100, "A가 소스 오른쪽: " + markerA);
  assert.ok(markerB > 100, "B가 소스 오른쪽: " + markerB);
  assert.ok(markerB - markerA >= 6, "간격≥6: " + (markerB - markerA));
});

test("computeMarkers 내부: 소스 오른쪽 끝 → 왼쪽 폴백", () => {
  const { markerA, markerB } = U.computeMarkers("internal", 150, 24,
    PLATE.plateColStart, PLATE.plateColEnd, PLATE.mouthI);
  // 오른쪽: 162,179 → clamp 154,154 → gap=0 < 6 → 왼쪽 폴백
  // markerA=150-12=138, markerB=150-29=121 → 모두 < 150
  assert.ok(markerA < 150, "A가 소스 왼쪽: " + markerA);
  assert.ok(markerB < 150, "B가 소스 왼쪽: " + markerB);
});

test("computeMarkers 외부(aperture): 마커가 입구 안쪽", () => {
  const { markerA, markerB } = U.computeMarkers("aperture", 10, 24,
    PLATE.plateColStart, PLATE.plateColEnd, PLATE.mouthI);
  assert.ok(markerA >= PLATE.mouthI, "A가 입구 이후: " + markerA);
  assert.ok(markerB > markerA, "B가 A보다 오른쪽: " + markerB);
  assert.ok(markerB <= PLATE.plateColEnd - 5, "B가 끝 이전: " + markerB);
});

test("computeMarkers: 두 레짐 모두 마커가 도파관 내부 범위", () => {
  [["internal", 100], ["aperture", 10]].forEach(([regime, src]) => {
    const { markerA, markerB } = U.computeMarkers(regime, src, 24,
      PLATE.plateColStart, PLATE.plateColEnd, PLATE.mouthI);
    assert.ok(markerA >= PLATE.plateColStart + 5,
      `${regime} A(${markerA}) < plateColStart+5`);
    assert.ok(markerB <= PLATE.plateColEnd - 5,
      `${regime} B(${markerB}) > plateColEnd-5`);
  });
});

done();
