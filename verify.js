// verify.js — RGD 물리 코어 검증 (console/fs/process는 여기서만 사용한다)
//
// 실행 순서(계산 순서 ≠ 파일 기록 조건):
//   1) sweep 계산 (60점 + 프리셋 2행) — 메모리에만 보관
//   2) 검증 8항목(1, 2, 3, 3′, 4, 5a, 5b, 6) 실행 — 3, 3′은 1의 sweep 데이터를 재사용
//   3) 전부 PASS일 때만 verify-results.json / sweep.csv를 기록. 하나라도 FAIL이면
//      파일을 쓰지 않고 process.exit(1)

'use strict';

var RGD = require('./physics.js');
var fs = require('fs');

var A = 1; // 반지름 고정
var ALPHA_TOTAL = 1;

function leastSquaresSlope(xs, ys) {
  var n = xs.length;
  var sx = 0, sy = 0, sxx = 0, sxy = 0;
  for (var i = 0; i < n; i++) {
    sx += xs[i];
    sy += ys[i];
    sxx += xs[i] * xs[i];
    sxy += xs[i] * ys[i];
  }
  return (n * sxy - sx * sy) / (n * sxx - sx * sx);
}

function fmtExp(v, digits) {
  return v.toExponential(digits === undefined ? 6 : digits);
}

// ---------------------------------------------------------------------
// STEP 1: sweep (메모리에만, 아직 파일에 쓰지 않는다)
// ---------------------------------------------------------------------

var SWEEP_N = 60;
var LAMBDA_A_MIN = Math.PI; // x = 2.0000 정확히 (λ/a=3.14로 두면 x=2.00064가 되어 x≤2 제한을 위반한다)
var LAMBDA_A_MAX = 1e3;

function computeRow(lambdaOverA) {
  var x = (2 * Math.PI) / lambdaOverA;
  var k = x; // a=1 이므로 k = 2π/λ = x
  var nd = RGD.chooseNd(x);
  var grid = RGD.buildDipoleGrid(A, nd);
  var r = RGD.sigma1D(grid, k, ALPHA_TOTAL, { n: 32 });
  return {
    lambda: lambdaOverA,
    x: x,
    nd: nd,
    N: grid.N,
    sigma: r.sigma,
    G: r.G,
    kh: r.kh,
    khOk: r.khOk
  };
}

var sweepRows = [];
(function buildSweep() {
  var logMin = Math.log(LAMBDA_A_MIN);
  var logMax = Math.log(LAMBDA_A_MAX);
  for (var i = 0; i < SWEEP_N; i++) {
    var t = i / (SWEEP_N - 1);
    var lambdaOverA = Math.exp(logMin + t * (logMax - logMin));
    sweepRows.push(computeRow(lambdaOverA));
  }
})();

var preset21cm = computeRow(2.1e6);
var presetVisible500nm = computeRow(5);

// ---------------------------------------------------------------------
// STEP 2: 검증 8항목
// ---------------------------------------------------------------------

var results = []; // {id, item, measured, measuredLabel, threshold, thresholdLabel, pass, group}

// --- 검증 1: 형상 인자 수렴 (절대 오차) ---
var TEST1_XS = [0.1, 0.5, 1.0, 2.0];
var test1PerX = TEST1_XS.map(function (x) {
  var k = x;
  var nd = RGD.chooseNd(x);
  var grid = RGD.buildDipoleGrid(A, nd);
  var N = grid.N;
  var maxErr = 0;
  for (var deg = 0; deg <= 180; deg++) {
    var theta = (deg * Math.PI) / 180;
    var q = RGD.qVector(k, theta, 0);
    var F = RGD.formFactor(grid.points, q);
    var Fmag = Math.sqrt(F.re * F.re + F.im * F.im) / N;
    var u = 2 * x * Math.sin(theta / 2);
    var f = RGD.analyticShapeFactor(u);
    var err = Math.abs(Fmag - f);
    if (err > maxErr) maxErr = err;
  }
  return { x: x, nd: nd, maxErr: maxErr };
});
var test1MaxErr = Math.max.apply(
  null,
  test1PerX.map(function (r) {
    return r.maxErr;
  })
);
var test1X2 = test1PerX[test1PerX.length - 1]; // x=2.0 — 이산화 신호 진단용
results.push({
  id: '1',
  item: '형상 인자 수렴 (절대오차, θ 181점)',
  measured: fmtExp(test1MaxErr, 4),
  threshold: '0.01 (절대)',
  pass: test1MaxErr <= 0.01,
  group: 'discretization'
});

// --- 검증 2: 전방 산란 불변 (입사파 위상 누락 검출) ---
var test2MaxDev = 0;
TEST1_XS.forEach(function (x) {
  var k = x;
  var nd = RGD.chooseNd(x);
  var grid = RGD.buildDipoleGrid(A, nd);
  var q = RGD.qVector(k, 0, 0);
  var F = RGD.formFactor(grid.points, q);
  var FN = Math.sqrt(F.re * F.re + F.im * F.im) / grid.N;
  var dev = Math.abs(FN - 1);
  if (dev > test2MaxDev) test2MaxDev = dev;
});
results.push({
  id: '2',
  item: '전방 산란 불변 (θ=0, F/N=1)',
  measured: fmtExp(test2MaxDev, 4),
  threshold: '1e-12',
  pass: test2MaxDev <= 1e-12,
  group: 'exact'
});

// --- 검증 3: 레일리 극한 — 기울기 ---
var test3Pts = sweepRows.filter(function (r) {
  return r.x < 0.05;
});
var test3Slope = leastSquaresSlope(
  test3Pts.map(function (r) {
    return Math.log(1 / r.lambda);
  }),
  test3Pts.map(function (r) {
    return Math.log(r.sigma);
  })
);
results.push({
  id: '3',
  item: '레일리 극한 — 기울기 (x<0.05, ' + test3Pts.length + '점 최소제곱)',
  measured: test3Slope.toFixed(6),
  threshold: '|기울기-4| < 0.01',
  pass: Math.abs(test3Slope - 4) < 0.01,
  group: 'discretization'
});

// --- 검증 3′: 레일리 극한 — 계수 (21cm 프리셋 재사용) ---
var test3primeGminus1 = preset21cm.G - 1;
results.push({
  id: "3'",
  item: '레일리 극한 — 계수 (21cm 프리셋, G→1.000)',
  measured: fmtExp(test3primeGminus1, 6) + '  [G=' + preset21cm.G.toFixed(9) + ']',
  threshold: '|G-1| < 1e-3',
  pass: Math.abs(test3primeGminus1) < 1e-3,
  group: 'exact'
});

// --- 검증 4: N 비의존성 (x=2.0, 최악 조건) ---
var TEST4_X = 2.0;
var TEST4_NDS = [20, 28, 40, 56];
function computeSigma1DAt(x, nd) {
  var k = x;
  var grid = RGD.buildDipoleGrid(A, nd);
  var r = RGD.sigma1D(grid, k, ALPHA_TOTAL, { n: 32 });
  return { nd: nd, N: grid.N, sigma: r.sigma };
}
var test4Rows = TEST4_NDS.map(function (nd) {
  return computeSigma1DAt(TEST4_X, nd);
});
var test4Mean =
  test4Rows.reduce(function (s, r) {
    return s + r.sigma;
  }, 0) / test4Rows.length;
var test4MaxDev = Math.max.apply(
  null,
  test4Rows.map(function (r) {
    return Math.abs(r.sigma - test4Mean) / test4Mean;
  })
);
results.push({
  id: '4',
  item: 'N 비의존성 (x=2.0, n_d=20/28/40/56)',
  measured: (test4MaxDev * 100).toFixed(4) + '%',
  threshold: '0.5%',
  pass: test4MaxDev <= 0.005,
  group: 'discretization'
});

// --- 검증 5a: 격자 이방성 (2D 대조, x=2.0, 최악 조건) ---
var TEST5_X = 2.0;
var test5Grid = RGD.buildDipoleGrid(A, RGD.chooseNd(TEST5_X));
var test5a1D = RGD.sigma1D(test5Grid, TEST5_X, ALPHA_TOTAL, { n: 32 });
var test5a2D = RGD.sigma2D(test5Grid, TEST5_X, ALPHA_TOTAL, { nTheta: 32, nPhi: 16 });
var test5aRel = Math.abs(test5a2D.sigma - test5a1D.sigma) / test5a1D.sigma;
results.push({
  id: '5a',
  item: '격자 이방성 — 2D 대조 (x=2.0)',
  measured: (test5aRel * 100).toFixed(4) + '%',
  threshold: '0.5%',
  pass: test5aRel <= 0.005,
  group: 'discretization'
});

// --- 검증 5b: 격자 이방성 (φ=0 vs φ=45°, 직접 진단, x=2.0) ---
var TEST5B_THETAS_DEG = [60, 90, 120];
var test5bPerTheta = TEST5B_THETAS_DEG.map(function (deg) {
  var theta = (deg * Math.PI) / 180;
  var q0 = RGD.qVector(TEST5_X, theta, 0);
  var F0 = RGD.formFactor(test5Grid.points, q0);
  var mag0 = Math.sqrt(F0.re * F0.re + F0.im * F0.im);
  var q45 = RGD.qVector(TEST5_X, theta, Math.PI / 4);
  var F45 = RGD.formFactor(test5Grid.points, q45);
  var mag45 = Math.sqrt(F45.re * F45.re + F45.im * F45.im);
  var dev = Math.abs(mag0 - mag45) / mag0;
  return { thetaDeg: deg, mag0: mag0, mag45: mag45, dev: dev };
});
var test5bMaxDev = Math.max.apply(
  null,
  test5bPerTheta.map(function (r) {
    return r.dev;
  })
);
results.push({
  id: '5b',
  item: '격자 이방성 — φ=0 vs φ=45° 직접 진단 (x=2.0, θ∈{60,90,120}°)',
  measured: (test5bMaxDev * 100).toFixed(4) + '%',
  threshold: '0.5%',
  pass: test5bMaxDev <= 0.005,
  group: 'discretization'
});

// --- 검증 6: 격자 중심 확인 (Im(F)/N ≈ 0) ---
var TEST6_ANGLES_DEG = [
  [60, 0],
  [90, 45],
  [120, 90],
  [30, 200]
];
var test6MaxIm = 0;
TEST6_ANGLES_DEG.forEach(function (pair) {
  var theta = (pair[0] * Math.PI) / 180;
  var phi = (pair[1] * Math.PI) / 180;
  var q = RGD.qVector(TEST5_X, theta, phi);
  var F = RGD.formFactor(test5Grid.points, q);
  var im = Math.abs(F.im) / test5Grid.N;
  if (im > test6MaxIm) test6MaxIm = im;
});
results.push({
  id: '6',
  item: '격자 중심 확인 (Im(F)/N ≈ 0)',
  measured: fmtExp(test6MaxIm, 4),
  threshold: '1e-12',
  pass: test6MaxIm <= 1e-12,
  group: 'exact'
});

// --- 검증 7: 원거리 한계 — 근접장 코드가 검증된 원거리장 코어와 일치하는가 ---
//
// scatteredFieldAt(완전한 쌍극자 장)으로 계산한 |E_sca|²의 각도 분포가
// 이미 검증된 dSigmaDOmega(θ,φ)와 일치하는지 본다. 이것으로 새 근접장
// 코드가 기존에 검증된 원거리장 코어에 묶인다.
//
// 목표값 = (k²·(alphaTotal/N)/r)² · dSigmaDOmega(points, k, θ, φ)
//   개별 분극률 alphaTotal/N 을 써야 한다 — alphaTotal을 그대로 쓰면
//   N²배(n_d=20에서 약 1.8e7배) 어긋난다.
//
// 관측 거리 r = max(30a, 100/k):
//   원거리 조건은 r ≫ a 뿐 아니라 kr ≫ 1 이다. λ/a=1000이면
//   k = 2π/1000 = 0.00628 이라 r=30a에서 kr = 0.19 로 오히려 준정전기
//   영역이 된다. 그래서 kr ≥ 100 을 강제하는 100/k 하한을 둔다.
//
// 평가 각도를 sin²Θ ≥ 0.25 로 제한하는 근거:
//   단일 쌍극자 장을 R̂/t̂ 정규직교 기저로 분해하면(u = kR, c = cosΘ, s = sinΘ)
//   횡성분 k²s[1 − 1/u² + i/u], 종성분 k²c[2/u² − 2i/u] 이므로
//     |E|² = k⁴[ sin²Θ·(1 − 1/(kR)²) + 4cos²Θ/(kR)² ] + O(1/(kR)⁴)
//     상대 오차 = (4cos²Θ/sin²Θ − 1) / (kR)²
//   φ=0 평면의 θ=90°는 n̂ = x̂ = 편광축이라 sin²Θ = 0 (쌍극자 널)이다.
//   목표값이 0인데 계산값은 0이 아니므로 상대 오차가 발산해 반드시 FAIL한다.
//   sin²Θ ≥ 0.25 로 제한하면 4cos²Θ/sin²Θ ≤ 12 이므로 최대 오차는
//   11/(kR)² = 0.11% 로, 임계값 1% 대비 9배 여유가 있다.
//     - φ=90° 평면: sin²Θ = 1 — 전 각도 안전
//     - φ=0  평면: sin²Θ = cos²θ — θ ≤ 60° 또는 θ ≥ 120° 만 사용
var TEST7_X_VALUES = [2.0, 0.0063]; // sweep 양 끝
var TEST7_ALPHA_TOTAL = 1;
var TEST7_ANGLES = [];
// φ=90° 평면: sin²Θ = 1 이므로 전 각도 안전
[20, 45, 60, 90, 120, 150].forEach(function (thetaDeg) {
  TEST7_ANGLES.push({ thetaDeg: thetaDeg, phiDeg: 90 });
});
// φ=0 평면: sin²Θ = cos²θ 이므로 θ ≤ 60° 또는 θ ≥ 120° 만
[20, 45, 60, 120, 150].forEach(function (thetaDeg) {
  TEST7_ANGLES.push({ thetaDeg: thetaDeg, phiDeg: 0 });
});

var test7Rows = [];
var test7MaxDev = 0;
TEST7_X_VALUES.forEach(function (x) {
  var nd = RGD.chooseNd(x);
  var grid = RGD.buildDipoleGrid(1, nd);
  var k = x; // a = 1 이므로 k = x/a = x
  var r = Math.max(30, 100 / k);
  var prefactor = (k * k * (TEST7_ALPHA_TOTAL / grid.N)) / r;

  TEST7_ANGLES.forEach(function (ang) {
    var theta = (ang.thetaDeg * Math.PI) / 180;
    var phi = (ang.phiDeg * Math.PI) / 180;
    var obs = {
      x: r * Math.sin(theta) * Math.cos(phi),
      y: r * Math.sin(theta) * Math.sin(phi),
      z: r * Math.cos(theta)
    };

    var f = RGD.scatteredFieldAt(obs, grid, k, TEST7_ALPHA_TOTAL);
    // 세 항을 모두 합친 전체 장으로 비교한다(원거리에서는 복사항이 지배하지만
    // 남은 근접장 성분까지 포함한 값이 위 오차식의 대상이다).
    var mag2 = 0;
    ['ex', 'ey', 'ez'].forEach(function (cn) {
      var re = f.radiation[cn].re + f.induction[cn].re + f.static[cn].re;
      var im = f.radiation[cn].im + f.induction[cn].im + f.static[cn].im;
      mag2 += re * re + im * im;
    });

    var target = prefactor * prefactor * RGD.dSigmaDOmega(grid.points, k, theta, phi);
    var dev = Math.abs(mag2 - target) / target;
    if (dev > test7MaxDev) test7MaxDev = dev;
    test7Rows.push({
      x: x,
      thetaDeg: ang.thetaDeg,
      phiDeg: ang.phiDeg,
      kr: k * r,
      computed: mag2,
      target: target,
      dev: dev
    });
  });
});
results.push({
  id: '7',
  item: '원거리 한계 (|E_sca|² vs dSigmaDOmega, r=max(30a,100/k))',
  measured: fmtExp(test7MaxDev, 4),
  threshold: '1e-2 (1%)',
  pass: test7MaxDev <= 1e-2,
  group: 'discretization'
});

// --- 검증 8: 정전기 극한 — 준정전기 두 항(정전·유도) 검증 ---
//
// 검증 7은 원거리이므로 복사항(k²/R)만 검증한다. 1/R²(유도), 1/R³(정전)
// 항은 전혀 검증되지 않는데, 패널 (A)(B)(C)의 존재 이유가 바로 이 두 항이다.
//
// 두 항목 모두 근사가 아니라 정확한 항등식이므로 임계값은 기계 정밀도(1e-12)다.
//
// (8a) static 항이 [3R̂(R̂·p̂) − p̂]·e^{ikR}/R³ 와 일치하는가.
//      기준식에 e^{ikR}를 반드시 포함한다 — static 항은 e^{ikR}를 이미 갖고
//      있으므로 기준식에서 빼면 |e^{ikR} − 1| ≈ kR = 1e-4 가 그대로 상대
//      오차로 남아, 구현이 아무리 정확해도 통과할 수 없는 검증이 된다.
// (8b) |induction|/|static| = kR. 두 항이 같은 벡터 인자 [3R̂(R̂·p̂) − p̂]와
//      같은 e^{ikR}를 공유하므로(계수만 (−ik/R) 대 (1/R²)) 이 비 역시
//      정확한 항등식이다.
var TEST8_KR = 1e-4;
// R̂·p̂ = 0, 0.5, 1 인 세 방향 × 여러 R
var TEST8_DIRECTIONS = [
  { name: 'R̂·p̂=0', ux: 0, uy: 0, uz: 1 },
  { name: 'R̂·p̂=0.5', ux: 0.5, uy: 0, uz: Math.sqrt(3) / 2 },
  { name: 'R̂·p̂=1', ux: 1, uy: 0, uz: 0 },
  { name: 'R̂·p̂=0.5(사선)', ux: 0.5, uy: Math.sqrt(3) / 2, uz: 0 }
];
var TEST8_RADII = [0.5, 2, 37];

var test8aMaxDev = 0;
var test8bMaxDev = 0;
var test8Rows = [];
TEST8_DIRECTIONS.forEach(function (dir) {
  TEST8_RADII.forEach(function (R) {
    var k = TEST8_KR / R;
    var dx = dir.ux * R;
    var dy = dir.uy * R;
    var dz = dir.uz * R;
    var kern = RGD.dipoleFieldKernel(dx, dy, dz, k);

    // 기준: [3R̂(R̂·p̂) − p̂]·e^{ikR}/R³   (p̂ = x̂ 이므로 R̂·p̂ = ux)
    var c = dir.ux;
    var v2 = [3 * dir.ux * c - 1, 3 * dir.uy * c, 3 * dir.uz * c];
    var kR = k * R;
    var cosKR = Math.cos(kR);
    var sinKR = Math.sin(kR);
    var R3 = R * R * R;

    var diff2 = 0;
    var ref2 = 0;
    var staticMag2 = 0;
    var indMag2 = 0;
    ['ex', 'ey', 'ez'].forEach(function (cn, m) {
      var refRe = (v2[m] * cosKR) / R3;
      var refIm = (v2[m] * sinKR) / R3;
      var s = kern.static[cn];
      diff2 += (s.re - refRe) * (s.re - refRe) + (s.im - refIm) * (s.im - refIm);
      // 규격화 분모는 e^{ikR} 없는 |[3R̂(R̂·p̂) − p̂]/R³| (크기는 동일)
      ref2 += (v2[m] / R3) * (v2[m] / R3);
      staticMag2 += s.re * s.re + s.im * s.im;
      var ind = kern.induction[cn];
      indMag2 += ind.re * ind.re + ind.im * ind.im;
    });

    var dev8a = Math.sqrt(diff2) / Math.sqrt(ref2);
    if (dev8a > test8aMaxDev) test8aMaxDev = dev8a;

    var ratio = Math.sqrt(indMag2) / Math.sqrt(staticMag2);
    var dev8b = Math.abs(ratio - kR) / kR;
    if (dev8b > test8bMaxDev) test8bMaxDev = dev8b;

    test8Rows.push({
      dir: dir.name,
      R: R,
      kR: kR,
      dev8a: dev8a,
      ratio: ratio,
      dev8b: dev8b
    });
  });
});
results.push({
  id: '8a',
  item: '정전기 극한 — static = [3R̂(R̂·p̂)−p̂]·e^{ikR}/R³ (kR=1e-4)',
  measured: fmtExp(test8aMaxDev, 4),
  threshold: '1e-12',
  pass: test8aMaxDev <= 1e-12,
  group: 'exact'
});
results.push({
  id: '8b',
  item: '유도항 k 스케일링 — |induction|/|static| = kR',
  measured: fmtExp(test8bMaxDev, 4),
  threshold: '1e-12',
  pass: test8bMaxDev <= 1e-12,
  group: 'exact'
});

// --- 검증 9: UI 최적화 경로(대칭)가 브루트포스와 일치하는가 ---
//
// 검증 7·8은 scatteredFieldAt(브루트포스)만 검증한다. 정작 화면을 그리는
// scatteredFieldXZPlane(y 대칭 2배)과 렌더 루프의 x 거울 관계는 아무것도
// 검증하지 않으므로, 나중에 성능 최적화를 건드릴 때 조용히 깨질 수 있다.
//
// 임계값은 1e-12가 아니라 1e-10으로 둔다: 두 경로는 덧셈 순서가 다르므로
// 약 1400항 합산에서 반올림 잡음이 2e-13까지 올라간다. 1e-12는 여유가
// 5배뿐이라 오검출 위험이 있다. 느슨하게 잡아도 잃는 것이 없다 —
// 대칭 버그는 미묘하지 않기 때문이다(×2 누락은 정확히 2배, 부호 오류는
// 200% 어긋남). 1e-13(반올림 잡음)과 실제 버그(O(1)) 사이는 텅 비어 있다.
//
// 9b의 한계: x 대칭은 script.js 렌더 루프에 있으므로 Node인 verify.js가
// 직접 검증할 수 없다. 9b가 검증하는 것은 관계식 자체이지, 렌더 루프가
// 그것을 올바로 적용했는지가 아니다. 렌더 루프의 부호 실수는 script.js의
// DEBUG_SYMMETRY 런타임 어서션으로 잡는다.
var TEST9_ND = 14; // 표시 전용 격자(1행 근접장 렌더링이 쓰는 값)
var TEST9_GRID = RGD.buildDipoleGrid(1, TEST9_ND);
var TEST9_X_VALUES = [2.0, 0.0063];
var TEST9_POINTS = [
  [0.5, 2.0],
  [2.5, -3.0],
  [-1.7, 0.3],
  [3.0, 5.0],
  [1.1, -0.8]
];
var TEST9_TERMS = ['radiation', 'induction', 'static'];

// (9a) y 대칭 최적화 검증 — scatteredFieldXZPlane vs scatteredFieldAt
var test9aMaxDev = 0;
TEST9_X_VALUES.forEach(function (k) {
  TEST9_POINTS.forEach(function (p) {
    var brute = RGD.scatteredFieldAt({ x: p[0], y: 0, z: p[1] }, TEST9_GRID, k, 1);
    var opt = RGD.scatteredFieldXZPlane(p[0], p[1], TEST9_GRID, k, 1);
    TEST9_TERMS.forEach(function (t) {
      // 항별 스케일로 규격화한다 — 성분 하나가 대칭 때문에 0이 되는 지점이
      // 있어(예: x=0에서 E_z=0) 성분별 상대 오차는 정의되지 않는다.
      var scale = 0;
      ['ex', 'ez'].forEach(function (c) {
        scale = Math.max(scale, Math.sqrt(brute[t][c].re * brute[t][c].re + brute[t][c].im * brute[t][c].im));
      });
      ['ex', 'ez'].forEach(function (c) {
        var dre = brute[t][c].re - opt[t][c].re;
        var dim = brute[t][c].im - opt[t][c].im;
        var dev = Math.sqrt(dre * dre + dim * dim) / scale;
        if (dev > test9aMaxDev) test9aMaxDev = dev;
      });
    });
  });
});
results.push({
  id: '9a',
  item: 'y 대칭 최적화 (scatteredFieldXZPlane vs 브루트포스)',
  measured: fmtExp(test9aMaxDev, 4),
  threshold: '1e-10',
  pass: test9aMaxDev <= 1e-10,
  group: 'exact'
});

// (9b) x 거울 관계 검증 — E_x(−x) = +E_x(x),  E_z(−x) = −E_z(x)
var test9bMaxDev = 0;
TEST9_X_VALUES.forEach(function (k) {
  TEST9_POINTS.forEach(function (p) {
    var pos = RGD.scatteredFieldXZPlane(p[0], p[1], TEST9_GRID, k, 1);
    var neg = RGD.scatteredFieldXZPlane(-p[0], p[1], TEST9_GRID, k, 1);
    TEST9_TERMS.forEach(function (t) {
      var scale = 0;
      ['ex', 'ez'].forEach(function (c) {
        scale = Math.max(scale, Math.sqrt(pos[t][c].re * pos[t][c].re + pos[t][c].im * pos[t][c].im));
      });
      // E_x는 짝함수: neg − pos = 0
      var dxr = neg[t].ex.re - pos[t].ex.re;
      var dxi = neg[t].ex.im - pos[t].ex.im;
      var devX = Math.sqrt(dxr * dxr + dxi * dxi) / scale;
      if (devX > test9bMaxDev) test9bMaxDev = devX;
      // E_z는 홀함수: neg + pos = 0
      var dzr = neg[t].ez.re + pos[t].ez.re;
      var dzi = neg[t].ez.im + pos[t].ez.im;
      var devZ = Math.sqrt(dzr * dzr + dzi * dzi) / scale;
      if (devZ > test9bMaxDev) test9bMaxDev = devZ;
    });
  });
});
results.push({
  id: '9b',
  item: 'x 거울 관계 (E_x 짝함수 / E_z 홀함수)',
  measured: fmtExp(test9bMaxDev, 4),
  threshold: '1e-10',
  pass: test9bMaxDev <= 1e-10,
  group: 'exact'
});

// --- 회귀 체크: sweep 전 구간에서 chooseNd(x)가 상수(20)로 유지되는가 ---
var ndConstant = sweepRows.every(function (r) {
  return r.nd === 20;
}) && preset21cm.nd === 20 && presetVisible500nm.nd === 20;
results.push({
  id: 'nd',
  item: 'sweep 전 구간 n_d 상수성 (계단 방지 회귀 체크)',
  measured: ndConstant ? '상수(20)' : '비상수 — 계단 발생 가능',
  threshold: '상수(20)',
  pass: ndConstant,
  group: 'exact'
});

// ---------------------------------------------------------------------
// 출력
// ---------------------------------------------------------------------

console.log('=== RGD 물리 코어 검증 ===\n');

console.table(
  results.map(function (r) {
    return {
      항목: r.id,
      설명: r.item,
      측정값: r.measured,
      임계값: r.threshold,
      결과: r.pass ? 'PASS' : 'FAIL'
    };
  })
);

console.log('\n검증 5b 세부 (θ별 |F(φ=0)| vs |F(φ=45°)| 상대 편차):');
console.table(
  test5bPerTheta.map(function (r) {
    return {
      'θ(deg)': r.thetaDeg,
      '|F(φ=0)|': r.mag0,
      '|F(φ=45°)|': r.mag45,
      '상대편차': fmtExp(r.dev, 6)
    };
  })
);
console.log('검증 5b 최댓값: ' + fmtExp(test5bMaxDev, 6));

console.log('\n검증 7 세부 (원거리 |E_sca|² vs dSigmaDOmega):');
console.table(
  test7Rows.map(function (r) {
    return {
      x: r.x,
      'θ(deg)': r.thetaDeg,
      'φ(deg)': r.phiDeg,
      kr: r.kr.toFixed(1),
      '계산값': fmtExp(r.computed, 6),
      '목표값': fmtExp(r.target, 6),
      '상대편차': fmtExp(r.dev, 4)
    };
  })
);
console.log('검증 7 최댓값: ' + fmtExp(test7MaxDev, 6));

console.log('\n검증 8 세부 (kR=1e-4 정전기 극한):');
console.table(
  test8Rows.map(function (r) {
    return {
      '방향': r.dir,
      R: r.R,
      kR: fmtExp(r.kR, 2),
      '8a 상대편차': fmtExp(r.dev8a, 4),
      '|ind|/|sta|': fmtExp(r.ratio, 6),
      '8b 상대편차': fmtExp(r.dev8b, 4)
    };
  })
);

var allPass = results.every(function (r) {
  return r.pass;
});

// --- 첫 실행 판정 지침: 실패 항목을 그룹별로 진단한다 ---
var failed = results.filter(function (r) {
  return !r.pass;
});

if (failed.length > 0) {
  console.log('\n=== 실패 진단 ===');

  failed.forEach(function (r) {
    if (r.group === 'exact') {
      console.log(
        '\n[' + r.id + '] ' + r.item + ' — 기계 정밀도 항목 실패. n_d를 건드리지 말고 코드를 볼 것 (확실한 버그).'
      );
      return;
    }

    // group === 'discretization' → n_d 스캔 자동 실행
    console.log('\n[' + r.id + '] ' + r.item + ' — 이산화 그룹 실패. n_d 스캔:');
    var scanRows = [20, 28, 40, 56].map(function (nd) {
      var value;
      if (r.id === '1') {
        var k = 2.0;
        var grid = RGD.buildDipoleGrid(A, nd);
        var N = grid.N;
        var maxErr = 0;
        for (var deg = 0; deg <= 180; deg++) {
          var theta = (deg * Math.PI) / 180;
          var q = RGD.qVector(k, theta, 0);
          var F = RGD.formFactor(grid.points, q);
          var Fmag = Math.sqrt(F.re * F.re + F.im * F.im) / N;
          var u = 2 * k * Math.sin(theta / 2);
          var f = RGD.analyticShapeFactor(u);
          var err = Math.abs(Fmag - f);
          if (err > maxErr) maxErr = err;
        }
        value = maxErr;
      } else if (r.id === '3') {
        // 회귀 검정은 sweep 전체 구조에 의존하므로 n_d 스캔 대상에서 제외
        value = null;
      } else if (r.id === '4') {
        var rows = TEST4_NDS; // n_d 자체가 스캔 변수이므로 그대로 재사용
        value = null;
      } else if (r.id === '5a') {
        var g = RGD.buildDipoleGrid(A, nd);
        var s1 = RGD.sigma1D(g, TEST5_X, ALPHA_TOTAL, { n: 32 });
        var s2 = RGD.sigma2D(g, TEST5_X, ALPHA_TOTAL, { nTheta: 32, nPhi: 16 });
        value = Math.abs(s2.sigma - s1.sigma) / s1.sigma;
      } else if (r.id === '5b') {
        var g2 = RGD.buildDipoleGrid(A, nd);
        var maxDev = 0;
        TEST5B_THETAS_DEG.forEach(function (deg) {
          var theta = (deg * Math.PI) / 180;
          var q0 = RGD.qVector(TEST5_X, theta, 0);
          var F0 = RGD.formFactor(g2.points, q0);
          var mag0 = Math.sqrt(F0.re * F0.re + F0.im * F0.im);
          var q45 = RGD.qVector(TEST5_X, theta, Math.PI / 4);
          var F45 = RGD.formFactor(g2.points, q45);
          var mag45 = Math.sqrt(F45.re * F45.re + F45.im * F45.im);
          var dev = Math.abs(mag0 - mag45) / mag0;
          if (dev > maxDev) maxDev = dev;
        });
        value = maxDev;
      }
      return { n_d: nd, 값: value === null ? '(해당없음)' : fmtExp(value, 4) };
    });
    console.table(scanRows);

    if (r.id === '5b') {
      console.log('[5b] 검증용 n_d를 늘려 0.5%를 만족하는 최소 n_d 탐색 중...');
      var found = null;
      for (var cand = 20; cand <= 400; cand += 8) {
        var gscan = RGD.buildDipoleGrid(A, cand);
        var maxDevScan = 0;
        TEST5B_THETAS_DEG.forEach(function (deg) {
          var theta = (deg * Math.PI) / 180;
          var q0 = RGD.qVector(TEST5_X, theta, 0);
          var F0 = RGD.formFactor(gscan.points, q0);
          var mag0 = Math.sqrt(F0.re * F0.re + F0.im * F0.im);
          var q45 = RGD.qVector(TEST5_X, theta, Math.PI / 4);
          var F45 = RGD.formFactor(gscan.points, q45);
          var mag45 = Math.sqrt(F45.re * F45.re + F45.im * F45.im);
          var dev = Math.abs(mag0 - mag45) / mag0;
          if (dev > maxDevScan) maxDevScan = dev;
        });
        if (maxDevScan <= 0.005) {
          found = cand;
          break;
        }
      }
      if (found) {
        console.log('[5b] 최소 n_d = ' + found + ' (상용 sweep의 n_d=20은 그대로 유지, 검증용 n_d와는 별개)');
      } else {
        console.log('[5b] n_d ≤ 400 범위에서 0.5%를 만족하는 n_d를 찾지 못함');
      }
    }
  });
}

console.log('\n' + (allPass ? 'PASS: 전체 13항목 + 회귀 체크 통과' : 'FAIL: 하나 이상의 항목이 임계값을 초과함'));

if (!allPass) {
  process.exit(1);
}

// ---------------------------------------------------------------------
// STEP 3: 전부 PASS일 때만 파일 기록
// ---------------------------------------------------------------------

var resultsJson = {
  params: {
    a: A,
    alphaTotal: ALPHA_TOTAL,
    E0: 1,
    timeConvention: 'e^{-i*omega*t}',
    xRange: [0.0063, 2],
    sweepRange: { lambdaOverA: [LAMBDA_A_MIN, LAMBDA_A_MAX], n: SWEEP_N },
    generatedAt: new Date().toISOString()
  },
  checks: results.map(function (r) {
    return { id: r.id, item: r.item, measured: r.measured, threshold: r.threshold, pass: r.pass };
  })
};
fs.writeFileSync('verify-results.json', JSON.stringify(resultsJson, null, 2));

var csvLines = ['lambda,x,nd,N,sigma,G'];
sweepRows.concat([preset21cm, presetVisible500nm]).forEach(function (r) {
  csvLines.push([r.lambda, r.x, r.nd, r.N, r.sigma, r.G].join(','));
});
fs.writeFileSync('sweep.csv', csvLines.join('\n') + '\n');

console.log('\nverify-results.json, sweep.csv 기록 완료.');
