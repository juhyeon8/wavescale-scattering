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

console.log('\n' + (allPass ? 'PASS: 전체 8항목 + 회귀 체크 통과' : 'FAIL: 하나 이상의 항목이 임계값을 초과함'));

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
