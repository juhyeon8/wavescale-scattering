# 도파관 입사·산란 중첩 FDTD 시뮬레이션 — 구현 계획

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 선원이 만드는 원통파(입사)와 평행판 도파관의 산란파를 FDTD로 계산해, 입사/산란/중첩 3분할 + 진폭 그래프로 차단·소멸파·입구 반사·감쇠상수를 보여주는 더블클릭 실행 웹 시뮬레이션을 만든다.

**Architecture:** 2D TM_z FDTD 격자 2개(판 없는 자유공간=입사, 판 있는=전체)를 동시에 돌리고 산란=전체−입사로 분해한다. 물리 계산은 순수 함수(`physics.js`)와 엔진(`fdtd.js`)·컨트롤러(`sim.js`)로 분리해 Node로 단위·통합 테스트하고, 렌더링/UI(`render.js`/`main.js`)는 브라우저에서 검증한다.

**Tech Stack:** Vanilla JS(ES5 호환 클래식 스크립트), HTML5 Canvas 2D, `Float32Array`. 테스트는 Node 내장 `node:assert`만 사용(외부 의존성/설치 없음). 브라우저는 ES 모듈·fetch 미사용 → 더블클릭 실행.

## Global Constraints

- 브라우저 전달은 모두 **클래식 `<script src>` / `<link>`** (ES 모듈·로컬 `fetch` 금지) → `index.html` 더블클릭으로 열림.
- 모든 공개 JS 파일은 **UMD 패턴**으로 작성: Node에서는 `module.exports`, 브라우저에서는 전역 `window.WaveSim` 네임스페이스에 부착. (테스트 가능 + 더블클릭 양립)
- 응답·주석·UI 문구는 **한국어**.
- 물리 편광: **E_z (TM_z)**, 도체판 경계조건 **E_z = 0**, 차단파장 **λ_c = 2a**.
- FDTD 내부 단위는 **셀 단위**(dx=dy=1). 표시값만 `DX_CM`로 cm 환산. 기본 `DX_CM = 0.5 cm/cell`.
- Courant 수 **S = 0.5** (2D 안정조건 S ≤ 1/√2 만족).
- **흡수경계: 분리장(Berenger) PML** (두께 ~10셀, 차수 m=3) — Mur 미사용. 두 격자에 **동일 PML**을 적용해 산란 패널(②)의 경계 인공물을 최소화한다.
- **차단 통일:** `k ≤ kc` ⇔ `λ ≥ 2a` 를 **소멸**로 처리한다(임계 `λ = 2a` 포함 → `κ = 0`). `k > kc` ⇔ `λ < 2a` 만 전파. `attenuationConstant` 는 `k > kc` 일 때만 `null` 을 반환.
- **위상자(정상상태 진폭) 측정은 애니메이션 루프와 분리:** 매 스텝 누적하되 **한 주기마다 스냅샷만** 읽는다. 루프 안에서 별도 스테핑 버스트(`measurePhasors`)를 돌려 순간장 애니메이션을 깨뜨리지 말 것.
- 위상 정의(소스·위상자 공통): `phase(n) = 2π · S · n / λ_cells`.
- 기본값: `a = 24 cells`(= 12 cm → λ_c = 24 cm), λ 슬라이더 `18–120 cells`(=9–60 cm, **저해상도 격자 이방성으로 원통파가 일그러지지 않도록 하한 상향**), 격자 `Nx=480, Ny=240`.
- 선원은 도체판과 **최소 2셀** 간격으로 클램프하고, 두 격자에서 **동일한 양자화 인덱스**를 사용.
- 커밋 메시지 접두사는 기존 저장소 관례(`기능:`, `수정:`, `테스트:`)를 따른다.

---

## 파일 구조

| 파일 | 책임 | 테스트 |
|------|------|--------|
| `index.html` | 마크업: 3개 캔버스 패널 + 그래프 캔버스 + 우측 조절판 | 브라우저 수동 |
| `style.css` | 레이아웃·색상(기존 자료 톤) | 브라우저 수동 |
| `physics.js` | 순수 함수: 차단/영역/감쇠상수/모드수/SWR→반사율/모드여기/지수fit | `tests/physics.test.js` |
| `fdtd.js` | `FDTD2D` 엔진: 필드 업데이트·Mur 경계·PEC·소스·위상자 누적 | `tests/fdtd.test.js` |
| `sim.js` | `WaveSim` 컨트롤러: 두 격자 운용·산란 분해·측정(차단/반사/κ) | `tests/sim.test.js` |
| `render.js` | 캔버스 드로잉: 컬러맵·패널·도파관·그래프 (컬러맵은 순수 함수로 테스트) | 일부 단위 + 브라우저 |
| `main.js` | UI 배선·애니메이션 루프·마우스·라벨/배지(순수 헬퍼는 테스트) | 일부 단위 + 브라우저 |
| `tests/_harness.js` | 미니 테스트 러너(`test`/`done`) | — |

브라우저 로드 순서: `physics.js → fdtd.js → sim.js → render.js → main.js`.

---

## Task 1: 테스트 하니스 + 물리 순수 함수(차단·영역·반사)

**Files:**
- Create: `tests/_harness.js`
- Create: `physics.js`
- Test: `tests/physics.test.js`

**Interfaces:**
- Produces (`window.WaveSim.physics` / `require('../physics.js')`):
  - `cutoffWavelength(aCells) -> number` (= `2*aCells`)
  - `regimeOf(sourceI, sourceJ, geom) -> 'aperture' | 'internal'` (geom=`{mouthI, jBot, jTop}`; x가 입구 안쪽 **그리고** y가 판 사이일 때만 `internal`)
  - `attenuationConstant(aCells, lambdaCells) -> number | null` (소멸 κ; `k > kc`(λ<2a)일 때만 `null`, 임계 λ=2a는 κ=0)
  - `propagatingModeCount(aCells, lambdaCells) -> number`
  - `swrToReflection(swr) -> number` (0..1)
  - `modeExcitationStrength(yCells, aCells, n) -> number` (0..1, `|sin(nπy/a)|`)

- [ ] **Step 1: 테스트 하니스 작성**

`tests/_harness.js`:
```js
"use strict";
let pass = 0, fail = 0;
function test(name, fn) {
  try { fn(); pass++; console.log("  ok   -", name); }
  catch (e) { fail++; console.error("  FAIL -", name, "\n   ", e.message); }
}
function approx(actual, expected, tol, msg) {
  if (Math.abs(actual - expected) > tol)
    throw new Error((msg || "approx") + `: got ${actual}, expected ${expected} ±${tol}`);
}
function done() {
  console.log(`\n${pass} passed, ${fail} failed`);
  if (fail) process.exit(1);
}
module.exports = { test, approx, done };
```

- [ ] **Step 2: 실패하는 테스트 작성**

`tests/physics.test.js`:
```js
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
done();
```

- [ ] **Step 3: 테스트 실패 확인**

Run: `node tests/physics.test.js`
Expected: FAIL (`Cannot find module '../physics.js'`)

- [ ] **Step 4: 최소 구현 작성**

`physics.js`:
```js
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
```

- [ ] **Step 5: 테스트 통과 확인**

Run: `node tests/physics.test.js`
Expected: PASS (`6 passed, 0 failed`)

- [ ] **Step 6: 커밋**

```bash
git add tests/_harness.js physics.js tests/physics.test.js
git commit -m "기능: 물리 순수 함수(차단·영역·반사·모드) + 테스트 하니스"
```

---

## Task 2: 지수함수 피팅 (감쇠상수 측정용)

**Files:**
- Modify: `physics.js`
- Test: `tests/physics.test.js`

**Interfaces:**
- Produces: `fitExponential(xs, ys) -> { kappa: number, amp: number }`
  - `ys[i] = amp · exp(-kappa · xs[i])` 모델을 `ln y` 선형회귀로 적합. `ys` 는 모두 양수 가정.

- [ ] **Step 1: 실패하는 테스트 추가** (`tests/physics.test.js` 의 `done()` 위에 삽입)

```js
test("fitExponential: 합성 지수데이터 복원", () => {
  const trueAmp = 2.5, trueKappa = 0.3;
  const xs = [], ys = [];
  for (let x = 0; x <= 20; x++) { xs.push(x); ys.push(trueAmp * Math.exp(-trueKappa * x)); }
  const r = P.fitExponential(xs, ys);
  approx(r.kappa, trueKappa, 1e-6, "kappa");
  approx(r.amp, trueAmp, 1e-5, "amp");
});
```

- [ ] **Step 2: 테스트 실패 확인**

Run: `node tests/physics.test.js`
Expected: FAIL (`P.fitExponential is not a function`)

- [ ] **Step 3: 구현 추가** (`physics.js` 의 `const API` 위에 함수 추가, API 객체에 등록)

```js
  function fitExponential(xs, ys) {
    const n = xs.length;
    let sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (let i = 0; i < n; i++) {
      const lx = xs[i], ly = Math.log(ys[i]);
      sx += lx; sy += ly; sxx += lx * lx; sxy += lx * ly;
    }
    const slope = (n * sxy - sx * sy) / (n * sxx - sx * sx);
    const intercept = (sy - slope * sx) / n;
    return { kappa: -slope, amp: Math.exp(intercept) };
  }
```
그리고 `const API = { ... , fitExponential, };` 로 등록.

- [ ] **Step 4: 테스트 통과 확인**

Run: `node tests/physics.test.js`
Expected: PASS (`7 passed, 0 failed`)

- [ ] **Step 5: 커밋**

```bash
git add physics.js tests/physics.test.js
git commit -m "기능: 지수함수 피팅(감쇠상수 측정)"
```

---

## Task 3: FDTD 엔진 코어 (`FDTD2D`)

**Files:**
- Create: `fdtd.js`
- Test: `tests/fdtd.test.js`

**Interfaces:**
- Produces: `class FDTD2D` (분리장 PML, TMz)
  - `new FDTD2D(Nx, Ny, courant=0.5)` — 무손실 계수로 초기화
  - `idx(i, j) -> number` (= `i*Ny + j`, i=x열, j=y행)
  - `setupPML(D=10, sigmaMax=0.8, m=3)` — 네 경계에 두께 D의 분리장 PML 계수를 grading
  - `reset()` — 모든 필드(`Ez`,`Ezx`,`Ezy`,`Hx`,`Hy`)/위상자/스텝수 0
  - `setPECMask(maskUint8)` — 길이 `Nx*Ny`, 1=도체. 매 스텝 해당 셀 `Ez=Ezx=Ezy=0`
  - `setSources(list)` — `[{i, j, weight}]`
  - `step({amp, lambdaCells, rampSteps})` — 한 스텝 전진(H→E(Ezx/Ezy)→소스→PEC), `this.n++`
  - `beginPhasor(lambdaCells)` / `accumulate()` — 위상자 복소 누적
  - `phasorComplex() -> {re, im, factor}` / `phasorAmp(out)` — 진폭맵 채움
  - 필드 접근: `this.Ez` (= `Ezx+Ezy`, Float32Array)

- [ ] **Step 1: 실패하는 테스트 작성**

`tests/fdtd.test.js`:
```js
"use strict";
const assert = require("node:assert");
const { test, approx, done } = require("./_harness.js");
const { FDTD2D } = require("../fdtd.js");

test("초기화: 필드 0, 인덱스 규칙", () => {
  const g = new FDTD2D(8, 6);
  assert.strictEqual(g.Ez.length, 48);
  assert.strictEqual(g.idx(2, 3), 2 * 6 + 3);
  assert.strictEqual(g.Ez[0], 0);
});

test("점 소스 → 좌우/상하 대칭장(자유공간)", () => {
  const N = 41, g = new FDTD2D(N, N);
  const c = (N - 1) / 2;
  g.setSources([{ i: c, j: c, weight: 1 }]);
  for (let n = 0; n < 20; n++) g.step({ amp: 1, lambdaCells: 12, rampSteps: 5 });
  // 중심 기준 대칭: (c+5,c) 와 (c-5,c) 거의 동일
  approx(g.Ez[g.idx(c + 5, c)], g.Ez[g.idx(c - 5, c)], 1e-4, "좌우대칭");
  approx(g.Ez[g.idx(c, c + 5)], g.Ez[g.idx(c, c - 5)], 1e-4, "상하대칭");
});

test("안정성: 200스텝 후 NaN/발산 없음", () => {
  const N = 60, g = new FDTD2D(N, N);
  g.setupPML(10, 0.8, 3); // 흡수경계로 공진 누적 없이 유계
  g.setSources([{ i: 30, j: 30, weight: 1 }]);
  for (let n = 0; n < 200; n++) g.step({ amp: 1, lambdaCells: 16, rampSteps: 10 });
  let maxAbs = 0;
  for (let k = 0; k < g.Ez.length; k++) {
    assert.ok(Number.isFinite(g.Ez[k]), "유한값");
    maxAbs = Math.max(maxAbs, Math.abs(g.Ez[k]));
  }
  assert.ok(maxAbs < 100, "발산하지 않음: " + maxAbs);
});

test("PEC 마스크 셀은 Ez=0 유지", () => {
  const N = 40, g = new FDTD2D(N, N);
  const mask = new Uint8Array(N * N);
  mask[g.idx(20, 10)] = 1;
  g.setPECMask(mask);
  g.setSources([{ i: 20, j: 20, weight: 1 }]);
  for (let n = 0; n < 30; n++) g.step({ amp: 1, lambdaCells: 12, rampSteps: 5 });
  assert.strictEqual(g.Ez[g.idx(20, 10)], 0);
});

test("위상자 누적: 정상상태 진폭 > 0", () => {
  const N = 50, g = new FDTD2D(N, N);
  g.setSources([{ i: 25, j: 25, weight: 1 }]);
  const lambda = 14;
  for (let n = 0; n < 60; n++) g.step({ amp: 1, lambdaCells: lambda, rampSteps: 10 });
  g.beginPhasor(lambda);
  const period = Math.round(lambda / 0.5); // T_steps = λ/S
  for (let n = 0; n < period; n++) { g.step({ amp: 1, lambdaCells: lambda, rampSteps: 10 }); g.accumulate(); }
  const amp = new Float32Array(N * N);
  g.phasorAmp(amp);
  assert.ok(amp[g.idx(27, 25)] > 0, "소스 근처 진폭 양수");
});

test("PML: 경계 반사 감소 — PML이 PEC 박스보다 후기 잔향 작음", () => {
  // 2D 점원은 하이겐스 후행파(wake)가 남아 절대 잔향은 0이 안 된다.
  // 따라서 '절대 잔향'이 아니라, 반사하는 PEC 박스 대비 흡수하는 PML의
  // 후기 잔향 감소로 흡수를 검증한다(공통 wake는 상쇄됨).
  function run(usePML) {
    const N = 80, g = new FDTD2D(N, N);
    if (usePML) g.setupPML(10, 0.8, 3);
    const c = 40; g.setSources([{ i: c, j: c, weight: 1 }]);
    const probe = g.idx(c + 20, c); // 경계 쪽 탐침
    let early = 0, late = 0;
    for (let n = 0; n < 250; n++) {
      const amp = n < 40 ? 1 : 0;     // 짧은 톤버스트 후 소스 끔
      g.step({ amp: amp, lambdaCells: 12, rampSteps: 5 });
      const v = Math.abs(g.Ez[probe]);
      if (n < 120) early = Math.max(early, v);
      else late = Math.max(late, v);
    }
    return { early: early, late: late };
  }
  const pec = run(false), pml = run(true);
  assert.ok(pec.early > 1e-3 && pml.early > 1e-3, "초기 통과 신호 존재");
  assert.ok(pml.late < 0.5 * pec.late,
    "PML 후기 잔향이 PEC 박스보다 작음: pml=" + pml.late + " pec=" + pec.late);
});
done();
```
> 상대 임계(0.5)를 못 넘으면 `sigmaMax`(0.4~1.5)나 두께 D(8~16)만 조정하라. 격자/소스/물리는 그대로. (2D wake 때문에 '절대 잔향' 측정은 부적절 — 반드시 PEC 대비 상대 측정.)

- [ ] **Step 2: 테스트 실패 확인**

Run: `node tests/fdtd.test.js`
Expected: FAIL (`Cannot find module '../fdtd.js'`)

- [ ] **Step 3: 엔진 구현**

`fdtd.js`:
```js
"use strict";
(function (global) {
  function FDTD2D(Nx, Ny, courant) {
    this.Nx = Nx; this.Ny = Ny;
    this.S = (courant === undefined) ? 0.5 : courant;
    const sz = Nx * Ny;
    this.Ezx = new Float32Array(sz);    // 분리장 x성분
    this.Ezy = new Float32Array(sz);    // 분리장 y성분
    this.Ez  = new Float32Array(sz);    // = Ezx + Ezy
    this.Hx = new Float32Array(sz);
    this.Hy = new Float32Array(sz);
    this.mask = null;
    this.sources = [];
    this.n = 0;
    this.accRe = new Float32Array(sz);
    this.accIm = new Float32Array(sz);
    this.accCount = 0;
    this.phasorLambda = 0;
    // 셀별 업데이트 계수(기본=무손실)
    this.caEx = new Float32Array(sz); this.cbEx = new Float32Array(sz); // Ezx (x손실)
    this.caEy = new Float32Array(sz); this.cbEy = new Float32Array(sz); // Ezy (y손실)
    this.daHx = new Float32Array(sz); this.dbHx = new Float32Array(sz); // Hx (y손실)
    this.daHy = new Float32Array(sz); this.dbHy = new Float32Array(sz); // Hy (x손실)
    this._setLossless();
  }
  FDTD2D.prototype.idx = function (i, j) { return i * this.Ny + j; };

  FDTD2D.prototype._setLossless = function () {
    const S = this.S, sz = this.Nx * this.Ny;
    for (let k = 0; k < sz; k++) {
      this.caEx[k] = 1; this.cbEx[k] = S;
      this.caEy[k] = 1; this.cbEy[k] = S;
      this.daHx[k] = 1; this.dbHx[k] = S;
      this.daHy[k] = 1; this.dbHy[k] = S;
    }
  };

  // 네 경계에 두께 D의 분리장(Berenger) PML. PEC로 종단(외곽 링은 갱신 안 함).
  FDTD2D.prototype.setupPML = function (D, sigmaMax, m) {
    D = D || 10; sigmaMax = (sigmaMax === undefined) ? 0.8 : sigmaMax; m = m || 3;
    this._setLossless();
    const Nx = this.Nx, Ny = this.Ny, S = this.S;
    function eCoef(sig) { return { a: (1 - sig / 2) / (1 + sig / 2), b: S / (1 + sig / 2) }; }
    // x방향 PML(좌·우) → Ezx, Hy
    for (let i = 0; i < Nx; i++) {
      let depth = 0;
      if (i < D) depth = (D - i) / D;
      else if (i >= Nx - D) depth = (i - (Nx - D - 1)) / D;
      if (depth <= 0) continue;
      const c = eCoef(sigmaMax * Math.pow(depth, m));
      for (let j = 0; j < Ny; j++) {
        const k = i * Ny + j;
        this.caEx[k] = c.a; this.cbEx[k] = c.b;
        this.daHy[k] = c.a; this.dbHy[k] = c.b;
      }
    }
    // y방향 PML(상·하) → Ezy, Hx
    for (let j = 0; j < Ny; j++) {
      let depth = 0;
      if (j < D) depth = (D - j) / D;
      else if (j >= Ny - D) depth = (j - (Ny - D - 1)) / D;
      if (depth <= 0) continue;
      const c = eCoef(sigmaMax * Math.pow(depth, m));
      for (let i = 0; i < Nx; i++) {
        const k = i * Ny + j;
        this.caEy[k] = c.a; this.cbEy[k] = c.b;
        this.daHx[k] = c.a; this.dbHx[k] = c.b;
      }
    }
  };

  FDTD2D.prototype.reset = function () {
    this.Ezx.fill(0); this.Ezy.fill(0); this.Ez.fill(0);
    this.Hx.fill(0); this.Hy.fill(0);
    this.accRe.fill(0); this.accIm.fill(0); this.accCount = 0; this.n = 0;
  };
  FDTD2D.prototype.setPECMask = function (mask) { this.mask = mask; };
  FDTD2D.prototype.setSources = function (list) { this.sources = list || []; };

  FDTD2D.prototype.step = function (p) {
    const Nx = this.Nx, Ny = this.Ny;
    const Ez = this.Ez, Ezx = this.Ezx, Ezy = this.Ezy, Hx = this.Hx, Hy = this.Hy;
    // --- H 업데이트 ---
    for (let i = 0; i < Nx; i++) {
      for (let j = 0; j < Ny - 1; j++) {
        const k = i * Ny + j;
        Hx[k] = this.daHx[k] * Hx[k] - this.dbHx[k] * (Ez[k + 1] - Ez[k]); // ∂Ez/∂y
      }
    }
    for (let i = 0; i < Nx - 1; i++) {
      for (let j = 0; j < Ny; j++) {
        const k = i * Ny + j;
        Hy[k] = this.daHy[k] * Hy[k] + this.dbHy[k] * (Ez[k + Ny] - Ez[k]); // ∂Ez/∂x
      }
    }
    // --- E 업데이트(외곽 링 제외 → PEC 종단) ---
    for (let i = 1; i < Nx - 1; i++) {
      for (let j = 1; j < Ny - 1; j++) {
        const k = i * Ny + j;
        Ezx[k] = this.caEx[k] * Ezx[k] + this.cbEx[k] * (Hy[k] - Hy[k - Ny]);
        Ezy[k] = this.caEy[k] * Ezy[k] - this.cbEy[k] * (Hx[k] - Hx[k - 1]);
        Ez[k] = Ezx[k] + Ezy[k];
      }
    }
    // --- 소스(소프트, Ezx/Ezy 절반씩) ---
    const ramp = Math.min(1, this.n / Math.max(1, p.rampSteps));
    const phase = 2 * Math.PI * this.S * this.n / p.lambdaCells;
    const drive = p.amp * ramp * Math.sin(phase);
    for (let s = 0; s < this.sources.length; s++) {
      const src = this.sources[s], k = src.i * Ny + src.j;
      Ezx[k] += 0.5 * drive * src.weight;
      Ezy[k] += 0.5 * drive * src.weight;
      Ez[k] = Ezx[k] + Ezy[k];
    }
    // --- PEC(도체판) ---
    if (this.mask) {
      const m = this.mask;
      for (let k = 0; k < m.length; k++) if (m[k]) { Ezx[k] = 0; Ezy[k] = 0; Ez[k] = 0; }
    }
    this.n++;
  };

  FDTD2D.prototype.beginPhasor = function (lambdaCells) {
    this.accRe.fill(0); this.accIm.fill(0); this.accCount = 0;
    this.phasorLambda = lambdaCells;
  };
  FDTD2D.prototype.accumulate = function () {
    const ph = 2 * Math.PI * this.S * this.n / this.phasorLambda;
    const c = Math.cos(ph), s = Math.sin(ph), Ez = this.Ez;
    const re = this.accRe, im = this.accIm;
    for (let k = 0; k < Ez.length; k++) { re[k] += Ez[k] * c; im[k] += Ez[k] * s; }
    this.accCount++;
  };
  FDTD2D.prototype.phasorComplex = function () {
    const f = this.accCount ? 2 / this.accCount : 0;
    return { re: this.accRe, im: this.accIm, factor: f };
  };
  FDTD2D.prototype.phasorAmp = function (out) {
    const f = this.accCount ? 2 / this.accCount : 0;
    const re = this.accRe, im = this.accIm;
    for (let k = 0; k < out.length; k++) out[k] = f * Math.hypot(re[k], im[k]);
    return out;
  };

  const API = { FDTD2D };
  if (typeof module !== "undefined" && module.exports) module.exports = API;
  else { global.WaveSim = global.WaveSim || {}; global.WaveSim.FDTD2D = FDTD2D; }
})(typeof globalThis !== "undefined" ? globalThis : this);
```

- [ ] **Step 4: 테스트 통과 확인**

Run: `node tests/fdtd.test.js`
Expected: PASS (`6 passed, 0 failed`)

- [ ] **Step 5: 커밋**

```bash
git add fdtd.js tests/fdtd.test.js
git commit -m "기능: 2D TM_z FDTD 엔진(분리장 PML·PEC·소스·위상자)"
```

---

## Task 4: 시뮬레이션 컨트롤러 + 산란 분해 (`WaveSim`)

**Files:**
- Create: `sim.js`
- Test: `tests/sim.test.js`

**Interfaces:**
- Consumes: `physics.js`(전부), `fdtd.js`(`FDTD2D`)
- Produces: `class WaveSim`
  - `new WaveSim(cfg)` where `cfg = { Nx, Ny, aCells, mouthI, plateColStart, plateColEnd, courant }`
  - `setSourceCell(i, j)` — 판과 최소 2셀 클램프 후 두 격자에 동일 인덱스로 소스 설정. 클램프된 `{i,j}` 반환
  - `setLambda(lambdaCells)` / `setAmp(amp)`
  - `step()` — 두 격자 1스텝(같은 파라미터)
  - `measurePhasors(periods=1)` — 두 격자를 직접 스테핑하며 위상자 누적 **(테스트 전용 배치; 런타임 애니메이션 루프에서는 호출 금지)**
  - `beginMeasure()` / `accumulateMeasure()` / `periodSteps()` — 런타임용: 루프가 매 스텝 `accumulateMeasure()`로 누적하고, `periodSteps()`마다 스냅샷 후 `beginMeasure()`로 재시작
  - `incidentAmp(out)` / `totalAmp(out)` / `scatteredAmp(out)` — 진폭맵(scattered=복소차의 크기)
  - `incidentFrame()` / `totalFrame()` / `scatteredFrame()` — 순간장 `Ez`(scattered=total.Ez−incident.Ez, Float32Array)
  - `centerlineAmp()` — 중심행(j=중심) totalAmp 1D 배열
  - `regime()` → `physics.regimeOf(sourceI, sourceJ, {mouthI, jBot, jTop})`
  - `cutoffInfo()` → `{ lambdaCells, lambdaCCells, modeCount, evanescent, kappa }`

- [ ] **Step 1: 실패하는 테스트 작성**

`tests/sim.test.js`:
```js
"use strict";
const assert = require("node:assert");
const { test, approx, done } = require("./_harness.js");
const { WaveSim } = require("../sim.js");

function baseCfg() {
  return { Nx: 120, Ny: 60, aCells: 24, mouthI: 40,
           plateColStart: 40, plateColEnd: 119, courant: 0.5 };
}

test("소스 클램프: 판과 최소 2셀", () => {
  const s = new WaveSim(baseCfg());
  // 판은 중심행 ±a/2. 판 바로 위에 두려 해도 2셀 안으로 못 들어감
  const jBot = s.jBot, jTop = s.jTop;
  const c = s.setSourceCell(80, jBot + 1);
  assert.ok(c.j >= jBot + 2 && c.j <= jTop - 2, "클램프됨: " + c.j);
});

test("두 격자 소스 인덱스 동일", () => {
  const s = new WaveSim(baseCfg());
  s.setSourceCell(80, 30);
  assert.deepStrictEqual(s.incident.sources, s.total.sources);
});

test("판이 없으면 산란≈0", () => {
  const cfg = baseCfg(); cfg.plateColStart = 0; cfg.plateColEnd = -1; // 판 없음
  const s = new WaveSim(cfg);
  s.setSourceCell(60, 30); s.setLambda(14);
  for (let n = 0; n < 80; n++) s.step();
  const out = new Float32Array(cfg.Nx * cfg.Ny);
  s.measurePhasors(2);
  s.scatteredAmp(out);
  let mx = 0; for (let k = 0; k < out.length; k++) mx = Math.max(mx, out[k]);
  assert.ok(mx < 1e-3, "산란 거의 0: " + mx);
});

test("산란 순간장 = 전체 − 입사", () => {
  const s = new WaveSim(baseCfg());
  s.setSourceCell(80, 30); s.setLambda(14);
  for (let n = 0; n < 30; n++) s.step();
  const sc = s.scatteredFrame();
  const k = s.total.idx(50, 30);
  approx(sc[k], s.total.Ez[k] - s.incident.Ez[k], 1e-6);
});

test("cutoffInfo: λ>2a 면 소멸·모드0", () => {
  const s = new WaveSim(baseCfg());
  s.setLambda(60); // λ=60 > 2a=48
  const ci = s.cutoffInfo();
  assert.strictEqual(ci.modeCount, 0);
  assert.strictEqual(ci.evanescent, true);
});
done();
```

- [ ] **Step 2: 테스트 실패 확인**

Run: `node tests/sim.test.js`
Expected: FAIL (`Cannot find module '../sim.js'`)

- [ ] **Step 3: 컨트롤러 구현**

`sim.js`:
```js
"use strict";
(function (global) {
  const isNode = (typeof module !== "undefined" && module.exports);
  const FDTD2D = isNode ? require("./fdtd.js").FDTD2D : global.WaveSim.FDTD2D;
  const P = isNode ? require("./physics.js") : global.WaveSim.physics;

  function WaveSim(cfg) {
    this.cfg = cfg;
    this.Nx = cfg.Nx; this.Ny = cfg.Ny; this.aCells = cfg.aCells;
    this.mouthI = cfg.mouthI;
    this.lambdaCells = 14; this.amp = 1; this.rampSteps = 20;
    // 판 위치: 중심행 기준 위/아래로 a/2
    this.jMid = Math.round(cfg.Ny / 2);
    this.jBot = this.jMid - Math.round(cfg.aCells / 2);
    this.jTop = this.jMid + Math.round(cfg.aCells / 2);
    this.sourceI = Math.round(cfg.Nx / 2); this.sourceJ = this.jMid;

    this.incident = new FDTD2D(cfg.Nx, cfg.Ny, cfg.courant);
    this.total = new FDTD2D(cfg.Nx, cfg.Ny, cfg.courant);
    const pml = cfg.pml || 10;
    this.incident.setupPML(pml); this.total.setupPML(pml); // 두 격자 동일 PML
    // 전체 격자에만 PEC 평행판 마스크
    const mask = new Uint8Array(cfg.Nx * cfg.Ny);
    if (cfg.plateColEnd >= cfg.plateColStart) {
      for (let i = cfg.plateColStart; i <= cfg.plateColEnd; i++) {
        mask[this.total.idx(i, this.jBot)] = 1;
        mask[this.total.idx(i, this.jTop)] = 1;
      }
    }
    this.total.setPECMask(mask);
    this._tmp = new Float32Array(cfg.Nx * cfg.Ny);
  }

  WaveSim.prototype.setSourceCell = function (i, j) {
    i = Math.max(1, Math.min(this.Nx - 2, Math.round(i)));
    // 판 사이 구간(도파관 안)일 때만 j 클램프; 자유공간 바깥은 자유
    if (i >= this.cfg.plateColStart && i <= this.cfg.plateColEnd &&
        this.cfg.plateColEnd >= this.cfg.plateColStart) {
      j = Math.max(this.jBot + 2, Math.min(this.jTop - 2, Math.round(j)));
    } else {
      j = Math.max(1, Math.min(this.Ny - 2, Math.round(j)));
    }
    this.sourceI = i; this.sourceJ = j;
    const src = [{ i: i, j: j, weight: 1 }];
    this.incident.setSources(src);
    this.total.setSources(src.map(function (s) { return { i: s.i, j: s.j, weight: s.weight }; }));
    return { i: i, j: j };
  };
  WaveSim.prototype.setLambda = function (l) { this.lambdaCells = l; };
  WaveSim.prototype.setAmp = function (a) { this.amp = a; };

  WaveSim.prototype.step = function () {
    const p = { amp: this.amp, lambdaCells: this.lambdaCells, rampSteps: this.rampSteps };
    this.incident.step(p); this.total.step(p);
  };
  WaveSim.prototype.measurePhasors = function (periods) {
    periods = periods || 1;
    const period = Math.round(this.lambdaCells / this.incident.S);
    this.incident.beginPhasor(this.lambdaCells);
    this.total.beginPhasor(this.lambdaCells);
    for (let n = 0; n < periods * period; n++) {
      this.step(); this.incident.accumulate(); this.total.accumulate();
    }
  };
  // 런타임용(애니메이션 루프에서 사용): 별도 스테핑 없이 누적만
  WaveSim.prototype.periodSteps = function () {
    return Math.round(this.lambdaCells / this.incident.S);
  };
  WaveSim.prototype.beginMeasure = function () {
    this.incident.beginPhasor(this.lambdaCells);
    this.total.beginPhasor(this.lambdaCells);
  };
  WaveSim.prototype.accumulateMeasure = function () {
    this.incident.accumulate(); this.total.accumulate();
  };

  WaveSim.prototype.incidentAmp = function (out) { return this.incident.phasorAmp(out); };
  WaveSim.prototype.totalAmp = function (out) { return this.total.phasorAmp(out); };
  WaveSim.prototype.scatteredAmp = function (out) {
    const ci = this.incident.phasorComplex(), ct = this.total.phasorComplex();
    for (let k = 0; k < out.length; k++) {
      const re = ct.factor * ct.re[k] - ci.factor * ci.re[k];
      const im = ct.factor * ct.im[k] - ci.factor * ci.im[k];
      out[k] = Math.hypot(re, im);
    }
    return out;
  };
  WaveSim.prototype.incidentFrame = function () { return this.incident.Ez; };
  WaveSim.prototype.totalFrame = function () { return this.total.Ez; };
  WaveSim.prototype.scatteredFrame = function () {
    const t = this.total.Ez, i = this.incident.Ez, out = this._tmp;
    for (let k = 0; k < out.length; k++) out[k] = t[k] - i[k];
    return out;
  };
  WaveSim.prototype.centerlineAmp = function () {
    this.total.phasorAmp(this._tmp);
    const row = new Float32Array(this.Nx);
    for (let i = 0; i < this.Nx; i++) row[i] = this._tmp[this.total.idx(i, this.jMid)];
    return row;
  };
  WaveSim.prototype.regime = function () {
    return P.regimeOf(this.sourceI, this.sourceJ,
      { mouthI: this.mouthI, jBot: this.jBot, jTop: this.jTop });
  };
  WaveSim.prototype.cutoffInfo = function () {
    const kappa = P.attenuationConstant(this.aCells, this.lambdaCells);
    return {
      lambdaCells: this.lambdaCells,
      lambdaCCells: P.cutoffWavelength(this.aCells),
      modeCount: P.propagatingModeCount(this.aCells, this.lambdaCells),
      evanescent: kappa !== null,
      kappa: kappa,
    };
  };

  const API = { WaveSim };
  if (isNode) module.exports = API;
  else { global.WaveSim = global.WaveSim || {}; global.WaveSim.WaveSim = WaveSim; }
})(typeof globalThis !== "undefined" ? globalThis : this);
```

> 주의: 브라우저에서 `global.WaveSim` 네임스페이스 객체와 `WaveSim` 클래스 이름이 겹친다. 브라우저 빌드에서는 클래스를 `global.WaveSim.Sim` 으로도 부착하라(아래 한 줄을 else 분기에 추가):
> ```js
> global.WaveSim.Sim = WaveSim;
> ```
> `main.js` 는 `WaveSim.Sim` 을 사용한다.

- [ ] **Step 4: 테스트 통과 확인**

Run: `node tests/sim.test.js`
Expected: PASS (`5 passed, 0 failed`)

- [ ] **Step 5: 커밋**

```bash
git add sim.js tests/sim.test.js
git commit -m "기능: WaveSim 컨트롤러(두 격자·산란 분해·측정)"
```

---

## Task 5: 통합 검증 — 차단 시 감쇠상수 κ 일치

순수 모드(횡방향 `sin(πy/a)`) 소스를 도파관 안쪽 한 열에 넣어 깨끗한 소멸파를 만들고, 중심선 진폭의 지수감쇠가 이론 κ와 맞는지 검증한다(엔진+분해+피팅 전체 통합).

**Files:**
- Test: `tests/integration.test.js`

**Interfaces:**
- Consumes: `WaveSim`(`total` 격자에 직접 모드 소스 주입), `physics.fitExponential`, `physics.attenuationConstant`

- [ ] **Step 1: 통합 테스트 작성**

`tests/integration.test.js`:
```js
"use strict";
const assert = require("node:assert");
const { test, approx, done } = require("./_harness.js");
const { WaveSim } = require("../sim.js");
const P = require("../physics.js");

test("내부여기 λ>2a: 중심선 |E| 지수감쇠가 이론 κ와 일치", () => {
  const Nx = 200, Ny = 60, a = 20;
  const cfg = { Nx, Ny, aCells: a, mouthI: 0,
                plateColStart: 0, plateColEnd: Nx - 1, courant: 0.5 };
  const s = new WaveSim(cfg);
  const lambda = 56; // λ>2a=40 → 소멸
  s.setLambda(lambda);
  // 순수 n=1 모드 소스를 i=30 열 전체(판 사이)에 sin(πy/a) 가중으로 주입
  const i0 = 30, src = [];
  for (let j = s.jBot + 1; j < s.jTop; j++) {
    const y = j - s.jBot;
    src.push({ i: i0, j: j, weight: Math.sin(Math.PI * y / a) });
  }
  s.total.setSources(src); s.incident.setSources([]);
  // 정상상태까지
  for (let n = 0; n < 600; n++) s.step();
  s.total.beginPhasor(lambda);
  const period = Math.round(lambda / 0.5);
  for (let n = 0; n < period; n++) { s.step(); s.total.accumulate(); }
  const amp = new Float32Array(Nx * Ny);
  s.total.phasorAmp(amp);
  // 소스에서 충분히 떨어진 구간(근접장 제외)에서 중심선 샘플
  const xs = [], ys = [], jMid = s.jMid;
  for (let i = i0 + 10; i <= i0 + 40; i++) { xs.push(i - i0); ys.push(amp[s.total.idx(i, jMid)]); }
  const fit = P.fitExponential(xs, ys);
  const kappaTheory = P.attenuationConstant(a, lambda);
  // FDTD 수치분산 감안, 상대오차 20% 이내
  approx(fit.kappa / kappaTheory, 1, 0.2, "κ 상대오차");
});

test("내부여기 λ<2a: 전파(감쇠 거의 없음)", () => {
  const Nx = 200, Ny = 60, a = 20;
  const cfg = { Nx, Ny, aCells: a, mouthI: 0,
                plateColStart: 0, plateColEnd: Nx - 1, courant: 0.5 };
  const s = new WaveSim(cfg);
  const lambda = 28; // λ<2a=40 → 전파
  s.setLambda(lambda);
  const i0 = 30, src = [];
  for (let j = s.jBot + 1; j < s.jTop; j++) {
    const y = j - s.jBot;
    src.push({ i: i0, j: j, weight: Math.sin(Math.PI * y / a) });
  }
  s.total.setSources(src); s.incident.setSources([]);
  for (let n = 0; n < 600; n++) s.step();
  s.total.beginPhasor(lambda);
  const period = Math.round(lambda / 0.5);
  for (let n = 0; n < period; n++) { s.step(); s.total.accumulate(); }
  const amp = new Float32Array(Nx * Ny);
  s.total.phasorAmp(amp);
  const a1 = amp[s.total.idx(i0 + 20, s.jMid)];
  const a2 = amp[s.total.idx(i0 + 50, s.jMid)];
  assert.ok(a2 > 0.4 * a1, "전파면 멀리서도 진폭 유지: " + a1 + " -> " + a2);
});
done();
```

- [ ] **Step 2: 테스트 실행 (실패 시 샘플창 조정)**

Run: `node tests/integration.test.js`
Expected: PASS (`2 passed, 0 failed`)

만약 첫 테스트가 허용오차를 넘으면 샘플 구간(`i0+10 ~ i0+40`)을 근접장에서 더 멀리/짧게 조정하거나 격자/도메인을 키워라(소스 너무 가까우면 근접장 오염, 너무 멀면 Mur 경계 잔향). 이론 κ 자체는 절대 바꾸지 말 것.

- [ ] **Step 3: 커밋**

```bash
git add tests/integration.test.js
git commit -m "테스트: 차단 시 감쇠상수 κ 이론 일치 통합 검증"
```

---

## Task 6: 렌더링 (`render.js`) — 컬러맵·패널·도파관·그래프

**Files:**
- Create: `render.js`
- Test: `tests/render.test.js`

**Interfaces:**
- Produces (`window.WaveSim.render`):
  - `divergingColor(value, scale) -> [r,g,b]` (순수: 음수=청, 0=백, 양수=적)
  - `paintField(ctx, field, Nx, Ny, scale, plates)` — 순간장/진폭맵을 ImageData로 그림
  - `drawGuide(ctx, geom)` — 판·입구·도선 마커
  - `drawGraph(ctx, row, geom, opts)` — |E| vs x 포락선 + 차단/전파 색
  - 컬러맵 외에는 브라우저 캔버스에서 수동 검증

- [ ] **Step 1: 컬러맵 단위 테스트 작성**

`tests/render.test.js`:
```js
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
```
> `render.js` 도 UMD로 작성하되, 캔버스 함수(`paintField` 등)는 Node에서 호출하지 않는다(테스트는 `divergingColor`만).

- [ ] **Step 2: 테스트 실패 확인**

Run: `node tests/render.test.js`
Expected: FAIL (`Cannot find module '../render.js'`)

- [ ] **Step 3: 렌더 구현**

`render.js` (컬러맵 + 캔버스 드로잉. UMD):
```js
"use strict";
(function (global) {
  function clamp(v, lo, hi) { return v < lo ? lo : (v > hi ? hi : v); }

  function divergingColor(value, scale) {
    const t = clamp(value / scale, -1, 1); // -1..1
    if (t >= 0) return [255, Math.round(255 * (1 - t)), Math.round(255 * (1 - t))];
    const a = -t;
    return [Math.round(255 * (1 - a)), Math.round(255 * (1 - a)), 255];
  }

  // field: Float32Array(Nx*Ny) (i=x열, j=y행). 캔버스는 좌상단 원점.
  function paintField(ctx, field, Nx, Ny, scale, geom) {
    const img = ctx.createImageData(Nx, Ny);
    const d = img.data;
    for (let j = 0; j < Ny; j++) {
      for (let i = 0; i < Nx; i++) {
        const v = field[i * Ny + j];
        const c = divergingColor(v, scale);
        // 화면 y는 위가 0 → j를 뒤집어 도파관 위/아래가 직관적이게
        const px = (j * Nx + i) * 4;
        d[px] = c[0]; d[px + 1] = c[1]; d[px + 2] = c[2]; d[px + 3] = 255;
      }
    }
    // geom.scaleX/scaleY 로 패널 크기에 맞춰 확대
    ctx.imageSmoothingEnabled = false;
    ctx.putImageData(img, 0, 0);
    if (geom && geom.zoom && geom.zoom !== 1) {
      ctx.drawImage(ctx.canvas, 0, 0, Nx, Ny, 0, 0, geom.zoom * Nx, geom.zoom * Ny);
    }
  }

  function drawGuide(ctx, geom) {
    const z = geom.zoom, Ny = geom.Ny;
    const yTop = (Ny - 1 - geom.jTop) * z, yBot = (Ny - 1 - geom.jBot) * z;
    ctx.strokeStyle = "#444"; ctx.lineWidth = 3;
    ctx.beginPath();
    ctx.moveTo(geom.plateColStart * z, yTop); ctx.lineTo(geom.plateColEnd * z, yTop);
    ctx.moveTo(geom.plateColStart * z, yBot); ctx.lineTo(geom.plateColEnd * z, yBot);
    ctx.stroke();
    // 도선 마커
    const sx = geom.sourceI * z, sy = (Ny - 1 - geom.sourceJ) * z;
    ctx.fillStyle = "#111"; ctx.beginPath();
    ctx.arc(sx, sy, 5, 0, Math.PI * 2); ctx.fill();
  }

  function drawGraph(ctx, row, geom, opts) {
    const W = ctx.canvas.width, H = ctx.canvas.height;
    ctx.clearRect(0, 0, W, H);
    let mx = 1e-9; for (let i = 0; i < row.length; i++) mx = Math.max(mx, row[i]);
    ctx.strokeStyle = opts && opts.evanescent ? "#c0392b" : "#2f6feb";
    ctx.lineWidth = 2; ctx.beginPath();
    for (let i = 0; i < row.length; i++) {
      const x = i / (row.length - 1) * W;
      const y = H - (row[i] / mx) * (H - 8) - 4;
      if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    }
    ctx.stroke();
    // 입구 위치 점선
    if (geom && geom.mouthI != null) {
      const xm = geom.mouthI / (row.length - 1) * W;
      ctx.strokeStyle = "#999"; ctx.setLineDash([4, 4]);
      ctx.beginPath(); ctx.moveTo(xm, 0); ctx.lineTo(xm, H); ctx.stroke();
      ctx.setLineDash([]);
    }
  }

  const API = { divergingColor, paintField, drawGuide, drawGraph };
  if (typeof module !== "undefined" && module.exports) module.exports = API;
  else { global.WaveSim = global.WaveSim || {}; global.WaveSim.render = API; }
})(typeof globalThis !== "undefined" ? globalThis : this);
```

- [ ] **Step 4: 테스트 통과 확인**

Run: `node tests/render.test.js`
Expected: PASS (`2 passed, 0 failed`)

- [ ] **Step 5: 커밋**

```bash
git add render.js tests/render.test.js
git commit -m "기능: 렌더링(발산 컬러맵·패널·도파관·그래프)"
```

---

## Task 7: HTML/CSS 골격 + UI 헬퍼(라벨·배지·드래그 클램프)

**Files:**
- Create: `index.html`
- Create: `style.css`
- Create: `main.js`
- Test: `tests/ui.test.js`

**Interfaces:**
- Produces (`window.WaveSim.ui`, main.js와 공유하는 순수 헬퍼):
  - `panel1Label(regime) -> string`
  - `regimeBadge(regime) -> string`
  - `cutoffBadge(info, regime) -> string`
  - `clampDragToCell(px, py, geom) -> {i, j}` (캔버스 px→셀, 판 클램프 위임 전 단계)

- [ ] **Step 1: UI 헬퍼 테스트 작성**

`tests/ui.test.js`:
```js
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
```
> 순수 헬퍼는 `ui.js`(UMD)로 분리해 Node 테스트 가능하게 하고, `main.js`는 이를 사용한다. 브라우저 로드 순서에 `ui.js`를 `main.js` 앞에 추가.

- [ ] **Step 2: 테스트 실패 확인**

Run: `node tests/ui.test.js`
Expected: FAIL (`Cannot find module '../ui.js'`)

- [ ] **Step 3: `ui.js` 구현**

`ui.js`:
```js
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
    const i = Math.round(px / geom.zoom);
    const j = Math.round(geom.Ny - 1 - py / geom.zoom);
    return { i: i, j: j };
  }
  const API = { panel1Label, regimeBadge, cutoffBadge, clampDragToCell };
  if (typeof module !== "undefined" && module.exports) module.exports = API;
  else { global.WaveSim = global.WaveSim || {}; global.WaveSim.ui = API; }
})(typeof globalThis !== "undefined" ? globalThis : this);
```

- [ ] **Step 4: 테스트 통과 확인**

Run: `node tests/ui.test.js`
Expected: PASS (`3 passed, 0 failed`)

- [ ] **Step 5: `index.html` 작성**

```html
<!DOCTYPE html>
<html lang="ko">
<head>
<meta charset="UTF-8">
<title>도파관 = 입사파 + 산란파 중첩 (FDTD)</title>
<link rel="stylesheet" href="style.css">
</head>
<body>
<div id="app">
  <div id="stage">
    <div class="panel"><div class="label" id="label1">① 입사파 (원통파)</div><canvas id="cv1"></canvas></div>
    <div class="panel"><div class="label">② 산란파 (경계조건이 강제)</div><canvas id="cv2"></canvas></div>
    <div class="panel"><div class="label">③ 중첩 (전체장)</div><canvas id="cv3"></canvas></div>
    <div class="panel"><div class="label">|E| vs x (도파관 중심선)</div><canvas id="cvGraph"></canvas></div>
  </div>
  <aside id="controls">
    <h1>도파관 FDTD</h1>
    <div class="badge" id="regimeBadge">선원: 자유공간 (개구 결합)</div>
    <div class="badge" id="cutoffBadge">도파관 내부: —</div>
    <label>파장 λ <span id="lambdaVal"></span></label>
    <input type="range" id="lambda" min="18" max="120" step="1" value="28">
    <label>도파관 간격 a <span id="aVal"></span></label>
    <input type="range" id="aGap" min="8" max="48" step="2" value="24">
    <label>진폭 <span id="ampVal"></span></label>
    <input type="range" id="amp" min="0.2" max="2" step="0.1" value="1">
    <div class="row">
      <button id="playBtn">⏸ 일시정지</button>
      <button id="resetBtn">↻ 리셋</button>
    </div>
    <label>속도(스텝/프레임) <span id="speedVal"></span></label>
    <input type="range" id="speed" min="1" max="8" step="1" value="3">
    <div class="info" id="kappaInfo"></div>
    <div class="info" id="reflInfo"></div>
    <p class="caption" id="modeCaption"></p>
  </aside>
</div>
<script src="physics.js"></script>
<script src="fdtd.js"></script>
<script src="sim.js"></script>
<script src="render.js"></script>
<script src="ui.js"></script>
<script src="main.js"></script>
</body>
</html>
```

- [ ] **Step 6: `style.css` 작성**

```css
* { box-sizing: border-box; }
html, body { margin: 0; height: 100%; font-family: "Malgun Gothic", sans-serif; background: #f7f7f8; color: #1c1c1f; }
#app { display: flex; height: 100vh; }
#stage { flex: 1; display: flex; flex-direction: column; min-width: 0; padding: 8px; gap: 6px; }
.panel { position: relative; flex: 1; background: #fff; border: 1px solid #d8d8dc; border-radius: 6px; overflow: hidden; }
.panel canvas { display: block; width: 100%; height: 100%; }
.label { position: absolute; top: 4px; left: 8px; font-size: 12px; color: #333; background: rgba(255,255,255,.75); padding: 1px 6px; border-radius: 4px; z-index: 1; }
#controls { width: 320px; flex-shrink: 0; background: #fff; border-left: 1px solid #d8d8dc; padding: 14px; overflow-y: auto; font-size: 13px; }
h1 { font-size: 15px; margin: 0 0 10px; }
.badge { font-size: 12px; padding: 5px 8px; border-radius: 6px; margin-bottom: 6px; background: #eef3ff; border: 1px solid #d6e0fb; }
label { display: block; margin: 12px 0 4px; font-size: 12px; }
input[type=range] { width: 100%; }
.row { display: flex; gap: 8px; margin-top: 12px; }
button { flex: 1; padding: 6px; border: 1px solid #d8d8dc; border-radius: 6px; background: #f0f0f3; cursor: pointer; }
.info { margin-top: 12px; font-size: 12px; line-height: 1.6; background: #f6f6f8; border: 1px solid #e2e2e6; border-radius: 6px; padding: 8px; }
.caption { margin-top: 12px; font-size: 11px; color: #6b6b72; line-height: 1.6; }
```

- [ ] **Step 7: `main.js` 작성 (배선·루프·마우스)**

```js
"use strict";
(function () {
  const NS = window.WaveSim;
  const DX_CM = 0.5; // cm per cell
  const Nx = 480, Ny = 240;
  const aCells = 24, mouthI = Math.round(Nx * 0.25); // 입구: 왼쪽 1/4 지점
  const cfg = { Nx, Ny, aCells, mouthI, plateColStart: mouthI, plateColEnd: Nx - 1, courant: 0.5 };
  const sim = new NS.Sim(cfg);
  sim.setLambda(28); sim.setSourceCell(Math.round(Nx * 0.12), sim.jMid);

  const cv1 = document.getElementById("cv1"), cv2 = document.getElementById("cv2"),
        cv3 = document.getElementById("cv3"), cvG = document.getElementById("cvGraph");
  [cv1, cv2, cv3].forEach(function (c) { c.width = Nx; c.height = Ny; });
  cvG.width = Nx; cvG.height = 120;
  const g1 = cv1.getContext("2d"), g2 = cv2.getContext("2d"),
        g3 = cv3.getContext("2d"), gG = cvG.getContext("2d");
  const geom = { Nx, Ny, zoom: 1, jBot: sim.jBot, jTop: sim.jTop,
                 plateColStart: cfg.plateColStart, plateColEnd: cfg.plateColEnd,
                 sourceI: sim.sourceI, sourceJ: sim.sourceJ, mouthI: mouthI };

  let playing = true, speed = 3, scale = 0.6;
  let measCount = 0, period = sim.periodSteps();
  let graphRow = new Float32Array(Nx);
  sim.beginMeasure();

  // 파라미터(λ·소스·a) 변경 시 위상자 측정 창을 재시작
  function onParamChange() { sim.beginMeasure(); measCount = 0; period = sim.periodSteps(); }

  function syncLabels() {
    const r = sim.regime(), info = sim.cutoffInfo();
    document.getElementById("label1").textContent = NS.ui.panel1Label(r);
    document.getElementById("regimeBadge").textContent = NS.ui.regimeBadge(r);
    document.getElementById("cutoffBadge").textContent = NS.ui.cutoffBadge(info, r);
    document.getElementById("lambdaVal").textContent = (sim.lambdaCells * DX_CM).toFixed(0) + " cm";
    document.getElementById("aVal").textContent = (sim.aCells * DX_CM).toFixed(0) + " cm (λc=" + (2 * sim.aCells * DX_CM).toFixed(0) + ")";
    // κ 검증(내부여기 & 소멸): 마지막 스냅샷 graphRow에서 좌·우 양방향 fit 평균
    const kInfo = document.getElementById("kappaInfo");
    if (r === "internal" && info.evanescent) {
      const si = sim.sourceI;
      const xsR = [], ysR = [], xsL = [], ysL = [];
      for (let i = si + 8; i <= Math.min(Nx - 12, si + 40); i++) { xsR.push(i - si); ysR.push(graphRow[i] + 1e-9); }
      for (let i = Math.max(11, si - 40); i <= si - 8; i++) { xsL.push(si - i); ysL.push(graphRow[i] + 1e-9); }
      const kR = xsR.length > 1 ? NS.physics.fitExponential(xsR, ysR).kappa : 0;
      const kL = xsL.length > 1 ? NS.physics.fitExponential(xsL, ysL).kappa : kR;
      const kMeas = (kR + kL) / 2;
      kInfo.textContent = "감쇠상수 κ — 측정(양방향 평균) " + (kMeas / DX_CM).toFixed(3) +
        " /cm vs 이론 " + (info.kappa / DX_CM).toFixed(3) + " /cm";
    } else { kInfo.textContent = ""; }
    // 모드 여기 캡션
    const y = sim.sourceJ - sim.jBot;
    const strength = NS.physics.modeExcitationStrength(y, sim.aCells, 1);
    document.getElementById("modeCaption").textContent =
      "n=1 모드 여기 세기(선원 y위치): " + (strength * 100).toFixed(0) +
      "% — 중심선=최대, 판 근처=약함. (이 약해짐은 차단과 다른 현상입니다.)";
  }

  function frame() {
    if (playing) {
      for (let s = 0; s < speed; s++) {
        sim.step();
        sim.accumulateMeasure();        // 애니메이션과 동일한 스테핑에 위상자 누적(별도 스테핑 없음)
        if (++measCount >= period) {    // 한 주기마다 스냅샷 후 재시작
          graphRow = sim.centerlineAmp();
          sim.beginMeasure(); measCount = 0;
        }
      }
    }
    NS.render.paintField(g1, sim.incidentFrame(), Nx, Ny, scale, geom);
    NS.render.paintField(g2, sim.scatteredFrame(), Nx, Ny, scale, geom);
    NS.render.paintField(g3, sim.totalFrame(), Nx, Ny, scale, geom);
    [g1, g2, g3].forEach(function (c) { NS.render.drawGuide(c, geomLive()); });
    NS.render.drawGraph(gG, graphRow, geom, { evanescent: sim.cutoffInfo().evanescent });
    syncLabels();
    requestAnimationFrame(frame);
  }
  function geomLive() {
    geom.sourceI = sim.sourceI; geom.sourceJ = sim.sourceJ; return geom;
  }

  // --- 컨트롤 ---
  document.getElementById("lambda").addEventListener("input", function (e) { sim.setLambda(+e.target.value); onParamChange(); });
  document.getElementById("amp").addEventListener("input", function (e) { sim.setAmp(+e.target.value); });
  document.getElementById("speed").addEventListener("input", function (e) { speed = +e.target.value; document.getElementById("speedVal").textContent = speed; });
  document.getElementById("playBtn").addEventListener("click", function (e) { playing = !playing; e.target.textContent = playing ? "⏸ 일시정지" : "▶ 재생"; });
  document.getElementById("resetBtn").addEventListener("click", function () { sim.incident.reset(); sim.total.reset(); });
  document.getElementById("aGap").addEventListener("input", function () { /* a 변경은 격자 재구성 필요: 간단히 페이지 새 cfg로 재생성 */ });

  // --- 마우스 드래그(도선 이동) ---
  let dragging = false;
  cv3.addEventListener("mousedown", function () { dragging = true; });
  window.addEventListener("mouseup", function () { dragging = false; });
  cv3.addEventListener("mousemove", function (e) {
    if (!dragging) return;
    const r = cv3.getBoundingClientRect();
    const px = (e.clientX - r.left) / r.width * Nx;
    const py = (e.clientY - r.top) / r.height * Ny;
    const cell = NS.ui.clampDragToCell(px, py, geom);
    sim.setSourceCell(cell.i, cell.j);
    onParamChange();
  });

  document.getElementById("speedVal").textContent = speed;
  requestAnimationFrame(frame);
})();
```

> 참고: `a`(간격) 슬라이더는 격자/마스크 재구성이 필요하므로, 본 계획에서는 자리만 두고 Task 8에서 “`a` 변경 시 `WaveSim` 재생성” 한 가지로 연결한다(작은 함수). 우선 동작하는 화면을 먼저 확보한다.

- [ ] **Step 8: 브라우저 수동 확인**

`index.html`을 더블클릭해 연다. 확인 항목:
- 세 패널에 적·청 파동이 보이고, ③ 중첩 = ① + ② 모양으로 일관.
- λ를 차단(>2a)로 올리면 도파관 내부로 파동이 못 들어가고 그래프가 지수감소.
- ③ 패널에서 도선을 드래그해 입구 안/밖으로 이동 시 배지·① 라벨이 바뀜.

- [ ] **Step 9: 커밋**

```bash
git add index.html style.css ui.js main.js tests/ui.test.js
git commit -m "기능: HTML/CSS 골격 + UI 배선(라벨·배지·드래그·κ·캡션)"
```

---

## Task 8: 마무리 — a 슬라이더 재구성 + 반사율 표시 + 기본값 튜닝

**Files:**
- Modify: `main.js`
- Modify: `sim.js` (반사율 계산 메서드)
- Test: `tests/sim.test.js` (반사율 메서드 테스트 추가)

**Interfaces:**
- Produces: `WaveSim.prototype.reflectionPercent() -> number | null`
  - 개구 결합 & 전파일 때만: 도파관 내부 중심선 진폭의 max/min(SWR)→`swrToReflection`→%. 아니면 `null`.

- [ ] **Step 1: 반사율 테스트 추가** (`tests/sim.test.js`)

```js
test("reflectionPercent: 정재파비에서 % 산출", () => {
  const s = new WaveSim(baseCfg());
  // 합성 중심선: 내부 구간에 정재파(min/max) 주입해 SWR 검증만
  s._fakeCenterline = function () {
    const row = new Float32Array(s.Nx);
    for (let i = 0; i < s.Nx; i++) row[i] = (i > s.mouthI) ? (2 + Math.cos(i * 0.5)) : 0; // max3 min1 → SWR3
    return row;
  };
  const pct = s.reflectionPercent(s._fakeCenterline());
  // SWR=3 → |R|=0.5 → 50%
  approx(pct, 50, 1, "반사율%");
});
```

- [ ] **Step 2: 테스트 실패 확인**

Run: `node tests/sim.test.js`
Expected: FAIL (`reflectionPercent is not a function`)

- [ ] **Step 3: `sim.js`에 반사율 메서드 추가**

```js
  WaveSim.prototype.reflectionPercent = function (centerline) {
    const row = centerline || this.centerlineAmp();
    let mn = Infinity, mx = 0;
    for (let i = this.mouthI + 5; i < this.Nx - 5; i++) {
      if (row[i] <= 0) continue;
      mn = Math.min(mn, row[i]); mx = Math.max(mx, row[i]);
    }
    if (!isFinite(mn) || mn <= 0) return null;
    const swr = mx / mn;
    return P.swrToReflection(swr) * 100;
  };
```

- [ ] **Step 4: 테스트 통과 확인**

Run: `node tests/sim.test.js`
Expected: PASS (`6 passed, 0 failed`)

- [ ] **Step 5: `main.js`에 반사율 표시 + a 재구성 연결**

`syncLabels()` 끝에 추가:
```js
    const refl = document.getElementById("reflInfo");
    if (sim.regime() === "aperture" && !sim.cutoffInfo().evanescent) {
      const pct = sim.reflectionPercent(graphRow);   // 스냅샷 포락선 사용
      refl.textContent = pct == null ? "" : "입구 반사율 ≈ " + pct.toFixed(0) + "%";
    } else if (sim.cutoffInfo().evanescent) {
      refl.textContent = "입구 반사 ≈ 거의 전반사 (전파 모드 없음)";
    } else { refl.textContent = ""; }
```

`a` 슬라이더 핸들러를 재구성으로 교체:
```js
  document.getElementById("aGap").addEventListener("change", function (e) {
    rebuildSim(+e.target.value);
  });
  function rebuildSim(aCellsNew) {
    const si = sim.sourceI, lam = sim.lambdaCells, amp = sim.amp;
    cfg.aCells = aCellsNew;
    const fresh = new NS.Sim(cfg);          // 새 격자는 생성자에서 PML 적용됨
    fresh.setLambda(lam); fresh.setAmp(amp); fresh.setSourceCell(si, fresh.jMid);
    sim = fresh;                            // 전역 sim 재할당(아래 주의 참고)
    geom.jBot = sim.jBot; geom.jTop = sim.jTop;
    onParamChange();                        // 측정 창 재시작
  }
```
> 구현 주의: Task 7의 `const sim = new NS.Sim(cfg);` 선언을 **`let sim`**으로 바꿔야 `rebuildSim`에서 재할당할 수 있다. `geom.jBot/jTop`도 위처럼 함께 갱신한다.

- [ ] **Step 6: 전체 테스트 재실행**

Run: `node tests/physics.test.js && node tests/fdtd.test.js && node tests/sim.test.js && node tests/integration.test.js && node tests/render.test.js && node tests/ui.test.js`
Expected: 모든 파일 PASS

- [ ] **Step 7: 브라우저 최종 확인 + 기본값 튜닝**

`index.html` 더블클릭. 차단/전파 전환, 도선 in/out, 반사율·κ 표시, a 슬라이더 재구성, 그래프 지수감소가 모두 자연스러운지 확인하고 `scale`(컬러 대비)·기본 λ 등을 보기 좋게 조정.

- [ ] **Step 8: 커밋**

```bash
git add main.js sim.js tests/sim.test.js
git commit -m "기능: 반사율 표시 + a 슬라이더 재구성 + 기본값 튜닝"
```

---

## Self-Review (계획 작성자 점검 결과)

**Spec coverage:**
- 1절 물리모델/영역/차단 → Task 1·4·5 ✓
- 2절 FDTD/두 격자/산란분해/Mur/소스 양자화 → Task 3·4 ✓
- 3절 세로 3단+그래프/패널 라벨/배지/마우스/모드 캡션 → Task 6·7 ✓
- 4절 차단 배지/입구 반사/κ 검증 → Task 5·7·8 ✓
- 5절 단위·기본값 → Global Constraints + Task 7·8 ✓
- 6절 3파일→기능별 분리(더블클릭 유지) → 파일 구조 절 ✓ (스펙과의 차이는 사용자에게 고지)
- 7절 YAGNI(모달 풀정량·영상 오버레이) → 계획에서 제외 ✓

**Placeholder scan:** 모든 코드 스텝에 실제 코드 포함. "적절히 처리" 류 없음. `a` 슬라이더 재구성만 Task 7에서 자리표시→Task 8에서 실제 구현으로 연결(의도된 분할).

**Type consistency:** `FDTD2D`(메서드 idx/step/beginPhasor/accumulate/phasorAmp/phasorComplex), `WaveSim`(setSourceCell/step/measurePhasors/centerlineAmp/cutoffInfo/reflectionPercent), 브라우저 클래스명 `WaveSim.Sim`, 네임스페이스 `window.WaveSim.{physics,FDTD2D,Sim,render,ui}` 일관 확인.

## 알려진 리스크 / 구현 중 조정 가능 지점

- Task 3 PML 흡수: `sigmaMax`(0.4~1.5)·두께 D(8~16)로 흡수 테스트(후기 잔향<10%)를 맞춤. 산란 패널(②)에 경계 인공물이 보이면 같은 파라미터를 키운다.
- Task 5 κ 통합 테스트는 근접장·경계 잔향에 민감 → 샘플창/도메인 크기 조정으로 허용오차(20%) 맞춤(이론값은 불변).
- 실시간 성능: 위상자는 애니메이션 스테핑에 함께 누적하고 한 주기마다 스냅샷만 읽는다(루프 내 별도 스테핑 없음). 격자 480×240×2가 무거우면 Nx/Ny를 낮춘다.
