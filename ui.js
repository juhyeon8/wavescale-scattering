"use strict";
(function (global) {
  function panel1Label(regime) {
    return regime === "internal"
      ? "① 자유공간 기준장 (판 없을 때)"
      : "① 입사파 (원통파)";
  }
  function regimeBadge(regime) {
    return regime === "internal"
      ? "선원: 도파관 내부 (내부 여기)"
      : "선원: 자유공간 (개구 결합)";
  }
  function cutoffBadge(info, regime) {
    const head = "도파관 내부: ";
    if (info.evanescent) {
      const tail = regime === "internal"
        ? "전파 모드 없음(소멸) — 선원 양쪽으로 소멸파만"
        : "전파 모드 없음(소멸)";
      return head + tail;
    }
    return head + "전파 모드 " + info.modeCount + "개";
  }
  function clampDragToCell(px, py, geom) {
    const i = Math.round(px / (geom.zoom || 1));
    const j = Math.round(geom.Ny - 1 - py / (geom.zoom || 1));
    return { i: i, j: j };
  }

  // 마커 A·B 배치: internal → 소스 같은 쪽, aperture → 입구 기준 내부 고정
  function computeMarkers(regime, sourceI, aCells, plateColStart, plateColEnd, mouthI) {
    var markerA, markerB;
    if (regime === "internal") {
      // 오른쪽 우선 배치
      markerA = sourceI + Math.max(10, Math.round(0.5 * aCells));
      markerB = sourceI + Math.max(20, Math.round(1.2 * aCells));
      markerA = Math.min(markerA, plateColEnd - 5);
      markerB = Math.min(markerB, plateColEnd - 5);
      // 오른쪽 공간 부족 → 왼쪽 폴백(같은 쪽 유지)
      if (markerB - markerA < 6) {
        markerA = sourceI - Math.max(10, Math.round(0.5 * aCells));
        markerB = sourceI - Math.max(20, Math.round(1.2 * aCells));
        markerA = Math.max(markerA, plateColStart + 5);
        markerB = Math.max(markerB, plateColStart + 5);
      }
    } else {
      // aperture: 입구로부터 내부 고정 위치
      var base = mouthI + Math.round(aCells * 0.5);
      markerA = Math.min(base, plateColEnd - 20);
      markerB = Math.min(markerA + Math.round(aCells * 0.7), plateColEnd - 5);
    }
    return { markerA: markerA, markerB: markerB };
  }

  const API = { panel1Label, regimeBadge, cutoffBadge, clampDragToCell, computeMarkers };
  if (typeof module !== "undefined" && module.exports) module.exports = API;
  else { global.WaveSim = global.WaveSim || {}; global.WaveSim.ui = API; }
})(typeof globalThis !== "undefined" ? globalThis : this);
