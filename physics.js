"use strict";
(function (global) {
  function cutoffWavelength(aCells) { return 2 * aCells; }

  function regimeOf(sourceI, sourceJ, geom) {
    const inside = sourceI >= geom.mouthI &&
                   sourceJ > geom.jBot && sourceJ < geom.jTop;
    return inside ? "internal" : "aperture";
  }

  function attenuationConstant(aCells, lambdaCells) {
    const k = 2 * Math.PI / lambdaCells;
    const kc = Math.PI / aCells;
    if (k > kc) return null;                          // 전파(λ<2a)
    return Math.sqrt(Math.max(0, kc * kc - k * k));   // 소멸(λ≥2a), 임계 λ=2a→0
  }

  function propagatingModeCount(aCells, lambdaCells) {
    // n < 2a/λ 인 정수 n≥1 의 개수
    return Math.max(0, Math.floor(2 * aCells / lambdaCells - 1e-9));
  }

  function swrToReflection(swr) {
    return (swr - 1) / (swr + 1);
  }

  function modeExcitationStrength(yCells, aCells, n) {
    return Math.abs(Math.sin(n * Math.PI * yCells / aCells));
  }

  const API = {
    cutoffWavelength, regimeOf, attenuationConstant,
    propagatingModeCount, swrToReflection, modeExcitationStrength,
  };
  if (typeof module !== "undefined" && module.exports) module.exports = API;
  else { global.WaveSim = global.WaveSim || {}; global.WaveSim.physics = API; }
})(typeof globalThis !== "undefined" ? globalThis : this);
