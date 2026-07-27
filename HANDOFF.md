# HANDOFF — RGD 소입자 산란 시뮬레이션

작성: 2026-07-27 · 기준 커밋 `56a26f8` (브랜치 `feature/waveguide-fdtd`)

> 다음 세션이 이 문서 하나만 읽고 이어받을 수 있게 쓴 인수인계 문서입니다.
> **§7(미해결)**과 **§8(구조 변경 분석)**이 핵심입니다.

---

## 1. 한 줄 요약

Rayleigh–Gans–Debye 근사로 소입자 산란을 계산하는 **물리 코어는 완성·검증 완료**(13항목 전부 PASS)이고,
그 위에 얹은 **탐색용 UI(탭 2개, 패널 6개)는 코드상 완성이지만 육안 검증을 한 번도 못 했습니다.**

무차원화 원칙: `a = 1` 고정, **λ가 유일한 조작 변인**, 조작 변수는 `x = ka` 하나뿐.

---

## 2. 파일 지도

| 파일 | 크기 | 역할 |
|---|---|---|
| `physics.js` | 24 KB | 물리 코어. 브라우저·Node 양쪽에서 쓰는 `RGD` 전역 IIFE. **DOM 의존 없음** |
| `verify.js` | 27 KB | 검증 13항목 + 회귀 체크. `console`/`fs`/`process`는 **여기서만** 사용 |
| `index.html` | 14 KB | 탭 2개 + 패널 6개 마크업 + ⓘ 설명 오버레이 본문 |
| `style.css` | 9 KB | 스크롤 없는 한 화면 레이아웃(`html,body{overflow:hidden}`) |
| `script.js` | 43 KB | 렌더링·상호작용 전부. 단일 IIFE, ES 모듈 아님 |
| `sweep.csv` | 5 KB | **논문 표 원본 데이터. 재생성만 하고 손으로 고치지 않음** |
| `verify-results.json` | 3 KB | 검증 결과 기록 (`generatedAt` 타임스탬프만 매 실행 갱신) |
| `DESIGN.md` / `README.md` | | 1단계 설계·개요 |
| `docs/superpowers/specs/2026-07-25-rgd-scattering-design.md` | | 1단계(물리 코어) 설계 전문 |
| `docs/superpowers/specs/2026-07-26-nearfield-panels-design.md` | | 2단계(근접장 패널) 설계 전문 |
| `RGD_1단계_프롬프트.md` / `RGD_2단계_프롬프트.md` | | 원 사양서. **2단계 쪽은 아직 git 미추적(untracked)** |

로딩 방식: `<link>` + 클래식 `<script src>`. **`index.html`을 더블클릭하면 그대로 열립니다**(ES 모듈·fetch 미사용).

---

## 3. 물리 코어 (`physics.js`) — 공개 API

```
chooseNd, buildDipoleGrid,
qVector, formFactor, formFactorPhi0Columns, analyticShapeFactor,
gaussLegendreNodes, dSigmaDOmega, sigma1D, sigma2D,
dipoleFieldKernel, scatteredFieldAt, scatteredFieldXZPlane
```

- **원거리 경로**(1단계): `formFactor` → `dSigmaDOmega` → `sigma1D`/`sigma2D` → σ, G
- **근접장 경로**(2단계): Jackson 9.18 완전 쌍극자 장을 **복사(1/R, k²) · 유도(1/R², k¹) · 정전(1/R³, k⁰)** 세 항으로 분해
  - `dipoleFieldKernel` — 쌍극자 하나의 세 항
  - `scatteredFieldAt` — 브루트포스 격자 합 (검증 기준)
  - `scatteredFieldXZPlane` — y 대칭을 이용한 최적화 경로 (**화면이 실제로 쓰는 것**)
- 시간 규약 `e^{-iωt}`, `E0 = 1`, `alphaTotal = 1`

### 검증 13항목 (`node verify.js`)

| # | 항목 | 측정 | 임계 |
|---|---|---|---|
| 1 | 형상 인자 수렴 (절대오차, θ 181점) | 2.69e-3 | 0.01 |
| 2 | 전방 산란 불변 (θ=0, F/N=1) | 0 | 1e-12 |
| 3 | 레일리 극한 — 기울기 | 3.999619 | \|기울기−4\| < 0.01 |
| 3′ | 레일리 극한 — 계수 (21cm, G→1) | G = 1.000000000 | \|G−1\| < 1e-3 |
| 4 | N 비의존성 (x=2.0, n_d=20/28/40/56) | 0.2376% | 0.5% |
| 5a | 격자 이방성 — 2D 대조 | 0.0280% | 0.5% |
| 5b | 격자 이방성 — φ=0 vs φ=45° | 0.1268% | 0.5% |
| 6 | 격자 중심 확인 (Im(F)/N ≈ 0) | 4.71e-16 | 1e-12 |
| 7 | 원거리 한계 (근접장 vs dσ/dΩ) | 2.55e-3 | 1e-2 |
| 8a | 정전기 극한 (kR=1e-4) | 4.46e-16 | 1e-12 |
| 8b | 유도항 k 스케일링 \|ind\|/\|sta\| = kR | 2.71e-16 | 1e-12 |
| 9a | y 대칭 최적화 vs 브루트포스 | 3.38e-15 | 1e-10 |
| 9b | x 거울 관계 (E_x 짝 / E_z 홀) | 2.37e-15 | 1e-10 |
| nd | sweep 전 구간 n_d 상수성 (계단 방지 회귀) | 상수(20) | 상수(20) |

**2026-07-27 재실행 결과: 전부 PASS.** `sweep.csv`는 재생성 후에도 git diff 없음(바이트 동일).

---

## 4. 현재 UI 구조

### 컨트롤 스트립 (탭 위 공용 줄)
- λ/a 슬라이더 (`#lambdaSlider`, 로그 매핑, λ ∈ [π, 1000] → x ∈ [0.0063, 2])
- 판독창 `#lambdaReadout` — **두 탭 공통(항상 표시)**
- 프리셋 5개: 21cm 전파 / 가시광선 500nm / 경계 12 / 경계 17 / 경계 40
- (B) 표시 항 `#termSelect` · 배율 `#scaleSelect`
- 프리셋·항 선택 그룹은 `.ctl-group[data-tab="1"]` → **탭 1에서만 표시**

### 탭 1 — "파동은 어떻게 통과하는가" (`#tab1`, `.grid2x2`)
읽기 순서 = 인과 순서로 배치:

| | 왼쪽 | 오른쪽 |
|---|---|---|
| 위 | **(A) 입사파** `Re(E_x), t=0` | **(D) 내부 전하의 진동** — 정사각 지도 + 범례 |
| 아래 | **(B) 산란파** — 항 토글·배율·수치 | **(C) 중첩 결과** — 배율 없이 원본값 합 |

### 탭 2 — "산란파는 어느 방향으로 얼마나 가는가" (`#tab2`, `.grid1x2`)
좌 = 원인, 우 = 결과:

- **왼쪽 산란 패턴** — dσ/dΩ (θ=0 규격화), 두 곡선(수직면 φ=90° / 편광면 φ=0)
  - `#polarBtn` "방향 분포로 보기" ↔ "그래프로 보기" 직교/극좌표 토글
  - **캔버스 클릭·드래그로 θ 선택** (`mousedown` → `pickTheta`, `window`의 `mousemove`/`mouseup`)
- **오른쪽 위상자 합** — 누적합 경로(코르뉘 나선), φ=0 고정, 길이 1/N 정규화

### 상호작용 규칙
- ⓘ 오버레이: 한 번에 하나만 열림, 탭 전환 시 전부 닫힘
- λ/a는 두 탭이 공유 (탭 바꿔도 유지) · θ는 탭 2 전용, 기본 90°
- "수치" 토글은 세션 동안 유지
- **숨은 탭의 캔버스는 `clientWidth = 0`** 이라 그릴 수 없음 → 활성 탭만 그리고 전환 시 다시 그림

---

## 5. `script.js` 렌더 파이프라인

```
슬라이더/프리셋 → recompute() → compute() → render()
                                              ├─ renderTab1() → paint/blit (A)(B)(C) + renderD
                                              └─ renderTab2() → renderPhasorSum + renderScatterPattern
θ 변경(클릭·드래그) → renderThetaOnly()   ← 장 재계산 없음
탭 전환 → setTab() → render()             ← 장 재계산 없음
resize → 120ms 디바운스 → render()        ← 장 재계산 없음
```

### 격자 3종 (혼동 주의)
| 이름 | 값 | 쓰는 곳 |
|---|---|---|
| `ND_DISPLAY` | 14 | (A)(B)(C) 장 지도 렌더링 **전용** |
| `ND_PHYS` | 20 | (D) 내부 위상 지도, (E) 위상자 합 — σ/G 계산과 같은 물리 격자 |
| 화면 격자 | NZ=120 × NX=72 | 표시 영역 z∈[−5a,5a], x∈[−3a,3a] (픽셀 정사각) |

### 주요 전역 상태 (파일 상단)
```js
var lambda, theta, activeTab, term, scaleMode, showMetrics;
var grid /*nd=14*/, physGrid /*nd=20*/;
var fieldRe = { radiation, induction, static };  // Float64Array(NZ*NX) ×3
var incidentRe, maskFlag, metrics;
```

### 성능·정확도 장치
- **x 거울 대칭**: 절반 행만 계산하고 나머지는 복사 (`E_x` 짝함수라 부호 그대로). 검증 9b가 관계식을, `DEBUG_SYMMETRY=true`가 렌더 루프의 적용을 각각 확인
- **마스크 반경 `MASK_R = 1.5a` 고정** — 축소 금지. 1.2a로 줄이면 1/R³ 발산 때문에 단일 쌍극자가 전체 장의 13%를 차지해 알갱이 무늬가 생김(1.5a에서는 1~3%)
- **세 항 캐싱**: (B)의 항 토글이 재계산 없이 표시만 바꿈. 재계산은 λ가 바뀔 때만
- `paint`/`blit`: 오프스크린 캔버스 key 기반 재사용 + letterbox

---

## 6. 절대 깨면 안 되는 것

1. **`sweep.csv`는 한 글자도 바뀌면 안 됨** — 논문 표 원본. 가시광선 행 `G = 0.5563911380033072` (63행).
   `verify.js`가 재생성하므로 **커밋 전 반드시 `git status`로 diff 없음을 확인**
2. **σ/G 계산 경로(`sigma1D`/`sigma2D`, nd=20)는 건드리지 않음**
3. **`node verify.js` 13항목 전부 PASS 유지**
4. 무차원화: `a = 1` 고정, λ가 조작 변인, 조작 변수는 `x = ka` 하나
5. 마스크 반경 1.5a · (E)의 φ=0 고정
6. 더블클릭으로 열리는 3파일 구조 유지 (ES 모듈·`fetch` 금지)

---

## 7. 미해결 — **UI 육안 검증 0회**

**Chrome 확장이 두 세션 연속 연결되지 않아 화면을 실제로 본 적이 없습니다.**
(`tabs_context_mcp` → "Browser extension is not connected")

코드상으로는 아래 요소의 배선을 확인했지만, **실제 동작은 미확인**입니다:

| 확인 대상 | 코드 배선 | 육안 |
|---|---|---|
| (1) 탭 전환 | ✅ `setTab()` — 패널·컨트롤그룹 `hidden` 토글, ⓘ 닫기, `render()` 재호출 | ❌ |
| (2) 산란 곡선 클릭·드래그 θ 선택 | ✅ `pickTheta()` — 직교/극좌표 모드 분기, `scatterGeom` 가드 | ❌ |
| (3) "방향 분포로 보기" 극좌표 전환 | ✅ `polarBtn` — 라벨 토글 + `renderScatterPattern()` | ❌ |
| (4) 스크롤 없이 한 화면 | ✅ `overflow:hidden` + flex/grid, 캡션 높이 3.9em 고정 | ❌ |
| 이미지 크기 일치 (A)(B)(C) | ✅ 캡션 높이 고정 + 셀 내 letterbox | ❌ |
| dpr > 1 에서 캔버스 넘침 | ✅ `.panel canvas { flex:none; max-width/height:100% }` | ❌ |

**다음 세션 첫 작업**: 확장을 연결하거나(Chrome 재시작 + claude.ai 로그인 계정 일치 확인),
사용자가 직접 `index.html`을 열어 스크린샷을 붙여 주는 방식으로 위 6가지를 확인.
레이아웃 관련 수정은 CLAUDE.md 규칙에 따라 **미리보기 HTML로 먼저** 보여주고 반영.

---

## 8. 다음 단계 — 구조 변경 분석

후보 두 가지가 논의되어 있었습니다.

### (a) 탭 1의 (A)(B)(C)(D) 라벨 제거
탭 2와 같은 원칙(제목만 남기고 알파벳 식별자 제거)으로 통일. **구조 변경은 거의 없음** — `index.html`의 `<h2>` 문자열과 ⓘ 오버레이 제목, `script.js`의 `paint('A'…)` key만 내용 기반으로 바꾸면 됨. 리스크 낮음.

> ⚠️ 인수인계 시점에 사용자 메시지가 잘려 있었습니다: *"지난 세션에 「탭1 확정 시…"* 이후가 유실.
> **(a)를 진행하기 전에 이 조건을 사용자에게 다시 확인해야 합니다.**

### (b) 논문 출력 모드 — **여기가 구조 변경이 필요한 지점**
목표: 탭 1 = 장 분포 묶음, 탭 2 = 기제 묶음으로 분리하고, **"행 = 물리량, 열 = λ" 전치 배치**, λ_A/λ_B 같은 색.

현재 코드가 이를 막는 지점 4가지:

1. **λ가 단일 전역 변수**
   `var lambda`를 `compute()`·`renderTab1()`·`renderTab2()`가 직접 읽습니다.
   여러 λ를 한 화면에 나란히 놓으려면 **λ가 인자로 흐르는 형태**로 바꿔야 합니다
   (`compute(lambda)` → 결과 객체 반환, `renderPanel(target, fieldSet)`).

2. **장 캐시가 λ 하나분만 존재**
   `fieldRe`/`incidentRe`/`maskFlag`가 모듈 스코프의 단일 `Float64Array` 세트입니다.
   열마다 다른 λ를 동시에 그리려면 **λ별 캐시 세트**(맵 또는 배열)로 승격해야 합니다.
   메모리는 문제없음 (120×72×3×8B ≈ 207 KB/λ).

3. **캔버스 id가 하드코딩된 단수형**
   `canvasA`/`canvasB`/`canvasC`/`canvasD`/`canvasScatter`/`canvasPhasor`를
   `document.getElementById`로 직접 잡습니다. 열이 N개가 되면 **패널을 DOM에서 생성**하거나
   최소한 id 대신 **요소 참조를 인자로 받는 형태**로 바꿔야 합니다.

4. **레이아웃이 `.grid2x2` / `.grid1x2` 두 클래스에 고정**
   전치 배치는 `grid-template-columns: repeat(λ개수, 1fr)`가 필요하므로
   CSS 클래스가 아니라 **인라인 grid-template 또는 CSS 변수**로 가야 합니다.
   더불어 화면 모드와 출력 모드는 캡션·ⓘ·툴버튼 노출 정책이 다르므로
   `body[data-mode="screen"|"print"]` 같은 **최상위 모드 스위치**가 자연스럽습니다.

**권장 순서**: (1)(2)를 먼저 하는 순수 리팩터링 커밋 → 화면 동작이 그대로임을 확인 →
그 위에 (3)(4)로 출력 모드를 얹기. 물리 코어(`physics.js`)는 손댈 필요 없습니다.

---

## 9. 사용자에게 확인이 필요한 열린 질문

1. **§7 브라우저 확인을 어떻게 진행할지** — 확장 재연결 시도 vs 사용자가 직접 열어 스크린샷
2. **§8의 (a)/(b) 중 무엇을 먼저 할지**
3. **(a)의 잘린 조건** — "탭1 확정 시 …" 이후 내용
4. **인수인계 메시지에서 잘린 두 곳**
   - *"(4) 스크롤 없이…"* → "한 화면에 들어오는지"로 해석했습니다. 맞는지 확인
   - *"커밋 전 git 인…"* → "git 인덱스/상태 확인"으로 해석했습니다. 맞는지 확인
5. **`RGD_2단계_프롬프트.md`를 git에 추가할지** (현재 untracked)

---

## 10. 현재 작업 트리 상태

```
 M verify-results.json        ← generatedAt 타임스탬프만 변경 (검증 재실행 부산물)
?? RGD_2단계_프롬프트.md      ← 미추적
```

`sweep.csv`는 검증 재실행 후에도 **변경 없음**(무결성 확인 완료).
