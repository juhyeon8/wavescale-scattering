"use strict";
(function () {
  var NS = window.WaveSim;
  var DX_CM = 0.5; // cm per cell
  var Nx = 480, Ny = 240;
  var aCells = 24, mouthI = Math.round(Nx * 0.25);
  var pml = 10;
  // plateColEnd: PML 바로 앞(469)까지만 → PML이 자동으로 오른쪽 끝 역할
  // wallThick: 벽 두께(셀) — 수치 누설 억제
  var cfg = { Nx: Nx, Ny: Ny, aCells: aCells, mouthI: mouthI,
              plateColStart: mouthI, plateColEnd: Nx - 1 - pml,
              courant: 0.5, pml: pml, wallThick: 2 };
  var sim = new NS.Sim(cfg);
  sim.setLambda(28); sim.setSourceCell(mouthI - 15, sim.jMid);

  var cv1 = document.getElementById("cv1"), cv2 = document.getElementById("cv2"),
      cv3 = document.getElementById("cv3"), cvG = document.getElementById("cvGraph"),
      cvTime = document.getElementById("cvTime");
  cv1.width = Nx; cv1.height = Ny;
  cv2.width = Nx; cv2.height = Ny;
  cv3.width = Nx; cv3.height = Ny;
  cvG.width = Nx; cvG.height = 120;
  cvTime.width = Nx; cvTime.height = 80;
  var g1 = cv1.getContext("2d"), g2 = cv2.getContext("2d"),
      g3 = cv3.getContext("2d"), gG = cvG.getContext("2d"),
      gT = cvTime.getContext("2d");
  var geom = { Nx: Nx, Ny: Ny, zoom: 1, jBot: sim.jBot, jTop: sim.jTop,
               plateColStart: cfg.plateColStart, plateColEnd: sim.plateColEnd,
               sourceI: sim.sourceI, sourceJ: sim.sourceJ, mouthI: mouthI,
               pml: pml, wallThick: cfg.wallThick || 2, jMid: sim.jMid };

  var playing = true, speed = 2, scale = 0.3;
  var measCount = 0, period = sim.periodSteps();
  var graphRow = new Float32Array(Nx);

  var markerA, markerB;
  function updateMarkers() {
    var m = NS.ui.computeMarkers(sim.regime(), sim.sourceI, sim.aCells,
      cfg.plateColStart, sim.plateColEnd, mouthI);
    markerA = m.markerA; markerB = m.markerB;
  }
  updateMarkers();

  var BUF = 1024;
  var ezAbuf = new Float32Array(BUF);
  var ezBbuf = new Float32Array(BUF);
  var bufHead = 0;

  // 위상차: 정확히 주기 경계에서만 저장(부분 주기 DFT 앨리어싱 방지)
  var lastDph = null;       // null = 아직 유효한 측정 없음
  var phasePeriods = 0;     // 마지막 파라미터 변경 이후 완료된 주기 수

  sim.beginMeasure();

  function onParamChange() {
    updateMarkers();
    lastDph = null; phasePeriods = 0;
    sim.beginMeasure(); measCount = 0; period = sim.periodSteps();
  }

  function syncLabels() {
    var r = sim.regime(), info = sim.cutoffInfo();
    document.getElementById("label1").textContent = NS.ui.panel1Label(r);
    document.getElementById("regimeBadge").textContent = NS.ui.regimeBadge(r);
    document.getElementById("cutoffBadge").textContent = NS.ui.cutoffBadge(info, r);
    document.getElementById("lambdaVal").textContent =
      (sim.lambdaCells * DX_CM).toFixed(0) + " cm";
    document.getElementById("aVal").textContent =
      (sim.aCells * DX_CM).toFixed(0) + " cm (λc=" + (2 * sim.aCells * DX_CM).toFixed(0) + " cm)";
    document.getElementById("ampVal").textContent = sim.amp.toFixed(1);

    var lbl2 = document.getElementById("label2");
    if (lbl2) lbl2.textContent = "② 산란파 (순간장, E_total − E_incident)";

    // κ 검증: 내부여기 & 소멸일 때 양방향 fit 평균
    var pml = cfg.pml || 10;
    var kInfo = document.getElementById("kappaInfo");
    if (r === "internal" && info.evanescent && info.kappa) {
      var si = sim.sourceI;
      var xsR = [], ysR = [], xsL = [], ysL = [];
      for (var i = si + 8; i <= Math.min(Nx - pml - 2, si + 40); i++) {
        xsR.push(i - si); ysR.push(graphRow[i] + 1e-9);
      }
      // 왼쪽 fit: 입구(mouthI)에서 충분히 떨어진 내부 구간만 (입구 방사·결합 영향 배제)
      for (var i2 = Math.max(mouthI + 10, si - 40); i2 <= si - 8; i2++) {
        xsL.push(si - i2); ysL.push(graphRow[i2] + 1e-9);
      }
      var kR = xsR.length > 1 ? NS.physics.fitExponential(xsR, ysR).kappa : 0;
      var kL = xsL.length > 1 ? NS.physics.fitExponential(xsL, ysL).kappa : kR;
      var kMeas = (kR + kL) / 2;
      kInfo.textContent = "감쇠상수 κ — 측정(양방향) " + (kMeas / DX_CM).toFixed(3) +
        " /cm vs 이론 " + (info.kappa / DX_CM).toFixed(3) + " /cm";
    } else if (r === "aperture" && info.evanescent && info.kappa) {
      var xsFit = [], ysFit = [];
      for (var ia = mouthI + 15; ia <= Nx - pml - 5; ia++) {
        xsFit.push(ia - mouthI); ysFit.push(graphRow[ia] + 1e-9);
      }
      if (xsFit.length > 2) {
        var kFit = NS.physics.fitExponential(xsFit, ysFit).kappa;
        kInfo.textContent = "감쇠상수 κ — 측정(입구→내부) " +
          (kFit / DX_CM).toFixed(3) + " /cm vs 이론 " +
          (info.kappa / DX_CM).toFixed(3) + " /cm";
      }
    } else {
      kInfo.textContent = "";
    }

    // 반사율: 개구 결합 + 전파일 때 SWR로 추정
    var reflEl = document.getElementById("reflInfo");
    if (r === "aperture" && !info.evanescent) {
      var pct = reflectionPercent(graphRow);
      reflEl.textContent = pct == null ? "" : "입구 반사율 ≈ " + pct.toFixed(0) + "%";
    } else if (info.evanescent) {
      reflEl.textContent = "입구 반사 ≈ 거의 전반사 (전파 모드 없음)";
    } else {
      reflEl.textContent = "";
    }

    // 마커 A-B 위상차: 정확한 주기 경계에서 저장된 값 표시
    var phaseEl = document.getElementById("phaseInfo");
    if (phaseEl && lastDph !== null) {
      phaseEl.textContent = "A-B 위상차: " + lastDph.toFixed(1) +
        "° — " + (Math.abs(lastDph) < 20 ? "제자리 진동(차단)" : "진행파");
    } else if (phaseEl) {
      phaseEl.textContent = "";
    }

    // 모드 여기 캡션
    var y = sim.sourceJ - sim.jBot;
    var strength = NS.physics.modeExcitationStrength(y, sim.aCells, 1);
    document.getElementById("modeCaption").textContent =
      "n=1 모드 여기 세기(선원 y위치): " + (strength * 100).toFixed(0) +
      "% — 중심선=최대, 판 근처=약함";
  }

  function reflectionPercent(row) {
    var pml = cfg.pml || 10;
    var mn = Infinity, mx = 0;
    for (var i = mouthI + 5; i < Nx - pml; i++) {
      if (row[i] <= 0) continue;
      if (row[i] < mn) mn = row[i];
      if (row[i] > mx) mx = row[i];
    }
    if (!isFinite(mn) || mn <= 0 || mx <= 0) return null;
    var swr = mx / mn;
    return NS.physics.swrToReflection(swr) * 100;
  }

  function geomLive() {
    geom.sourceI = sim.sourceI; geom.sourceJ = sim.sourceJ; return geom;
  }

  function frame() {
    if (playing) {
      for (var s = 0; s < speed; s++) {
        sim.step();
        sim.accumulateMeasure();
        var ezTot = sim.totalFrame();
        ezAbuf[bufHead % BUF] = ezTot[markerA * Ny + sim.jMid];
        ezBbuf[bufHead % BUF] = ezTot[markerB * Ny + sim.jMid];
        bufHead++;
        if (++measCount >= period) {
          graphRow = sim.centerlineAmp();
          phasePeriods++;
          window._dbgBufHead = bufHead;
          window._dbgPhasePeriods = phasePeriods;
          // 위상: 정확히 완료된 주기에서만 측정(앨리어싱 제로)
          // 과도기 감쇠 대기(λ=56 기준 4 주기 ≈ 3.7 sec)
          if (phasePeriods >= 2) {
            var cpx = sim.total.phasorComplex();
            var kA = markerA * Ny + sim.jMid, kB = markerB * Ny + sim.jMid;
            var reA = cpx.factor * cpx.re[kA], imA = cpx.factor * cpx.im[kA];
            var reB = cpx.factor * cpx.re[kB], imB = cpx.factor * cpx.im[kB];
            var d = (Math.atan2(imB, reB) - Math.atan2(imA, reA)) * 180 / Math.PI;
            while (d > 180) d -= 360; while (d < -180) d += 360;
            lastDph = d;
            window._dbgLastDph = d;
          }
          sim.beginMeasure(); measCount = 0;
        }
      }
    }
    NS.render.paintField(g1, sim.incidentFrame(), Nx, Ny, scale, geom);
    NS.render.paintField(g2, sim.scatteredFrame(), Nx, Ny, scale, geom);
    NS.render.paintField(g3, sim.totalFrame(), Nx, Ny, scale, geom);
    var gl = geomLive();
    NS.render.drawGuide(g1, gl);
    NS.render.drawGuide(g2, gl);
    NS.render.drawGuide(g3, gl);
    NS.render.drawMarkers(g3, markerA, markerB, geom);
    NS.render.drawGraph(gG, graphRow, geom, { evanescent: sim.cutoffInfo().evanescent });
    NS.render.drawTimeGraph(gT, ezAbuf, ezBbuf, bufHead, bufHead, period);
    syncLabels();
    requestAnimationFrame(frame);
  }

  // 컨트롤 배선
  document.getElementById("lambda").addEventListener("input", function (e) {
    sim.setLambda(+e.target.value); onParamChange();
  });
  document.getElementById("amp").addEventListener("input", function (e) {
    sim.setAmp(+e.target.value);
  });
  document.getElementById("speed").addEventListener("input", function (e) {
    speed = +e.target.value;
    document.getElementById("speedVal").textContent = speed;
  });
  document.getElementById("playBtn").addEventListener("click", function (e) {
    playing = !playing;
    e.target.textContent = playing ? "⏸ 일시정지" : "▶ 재생";
  });
  document.getElementById("resetBtn").addEventListener("click", function () {
    sim.incident.reset(); sim.total.reset(); sim.beginMeasure(); measCount = 0;
    bufHead = 0; ezAbuf.fill(0); ezBbuf.fill(0);
  });
  document.getElementById("aGap").addEventListener("change", function (e) {
    rebuildSim(+e.target.value);
  });

  function rebuildSim(aCellsNew) {
    var si = sim.sourceI, lam = sim.lambdaCells, a = sim.amp;
    cfg.aCells = aCellsNew;
    var fresh = new NS.Sim(cfg);
    fresh.setLambda(lam); fresh.setAmp(a); fresh.setSourceCell(si, fresh.jMid);
    sim = fresh;
    geom.jBot = sim.jBot; geom.jTop = sim.jTop;
    geom.plateColEnd = sim.plateColEnd;
    updateMarkers();
    bufHead = 0; ezAbuf.fill(0); ezBbuf.fill(0);
    onParamChange();
  }

  // 마우스 드래그(③ 패널에서 선원 이동)
  var dragging = false;
  cv3.addEventListener("mousedown", function () { dragging = true; });
  window.addEventListener("mouseup", function () { dragging = false; });
  cv3.addEventListener("mousemove", function (e) {
    if (!dragging) return;
    var r = cv3.getBoundingClientRect();
    var px = (e.clientX - r.left) / r.width * Nx;
    var py = (e.clientY - r.top) / r.height * Ny;
    var cell = NS.ui.clampDragToCell(px, py, { zoom: 1, Ny: Ny });
    sim.setSourceCell(cell.i, cell.j);
    geom.sourceI = sim.sourceI; geom.sourceJ = sim.sourceJ;
    onParamChange();
  });

  document.getElementById("speedVal").textContent = speed;
  // 디버그용 전역 노출
  window._dbgBufHead = 0;
  window._dbgPhasePeriods = 0;
  window._dbgLastDph = null;
  requestAnimationFrame(frame);
})();
