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

  const API = { panel1Label, regimeBadge, cutoffBadge, clampDragToCell };
  if (typeof module !== "undefined" && module.exports) module.exports = API;
  else { global.WaveSim = global.WaveSim || {}; global.WaveSim.ui = API; }
})(typeof globalThis !== "undefined" ? globalThis : this);
