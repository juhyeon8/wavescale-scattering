// verify-ui.js — 화면(축척·렌더링) 회귀 검증
//
// verify.js 는 물리 코어 전용이다(DOM 없이 Node 에서 도는 순수 계산).
// 줄무늬 주기·물체 크기는 "실제로 그려진 픽셀"을 재야 의미가 있으므로
// 여기 Playwright 쪽에 둔다. 두 파일의 역할은 겹치지 않는다:
//   verify.js    — G, σ, 형상 인자, 근접장 커널 (13항목)
//   verify-ui.js — 화면 축척, 물체 표식 크기, 모드 복원
//
// 실행:
//   npm install playwright && npx playwright install chromium
//   node verify-ui.js
//
// 검증 항목
//   1. SCALE_MODE='a' 에서 줄무늬 주기가 λ에 비례한다 (±5%)
//   2. 물체 표식 지름 ≤ 패널 높이의 20%, λ에 무관하게 일정
//   3. 물리 판독값(Q, G, |E_산란|/|E_입사|)이 축척 모드와 무관하다
//   4. SCALE_MODE='lambda' 로 되돌리면 (A)가 세 λ에서 동일해진다(원래 성질)

const { chromium } = require('playwright');
const path = require('path');
const url = require('url');

const PAGE = url.pathToFileURL(path.join(__dirname, 'index.html')).href;
const TOL = 0.05;

// 슬라이더는 1/1000 로 양자화되어 있다(한 칸 ≈ 1.6%). λ/a = π 는 적용 범위의
// 경계 그 자체라 ceil 로 한 칸 위를 잡아야 "범위 밖"으로 떨어지지 않는다.
const setLambda = ([lam, up]) => {
  const s = document.getElementById('lambdaSlider');
  const f = 1000 * Math.log(lam / 0.1) / Math.log(1e6 / 0.1);
  s.value = up ? Math.ceil(f) : Math.round(f);
  s.dispatchEvent(new Event('input', { bubbles: true }));
  s.dispatchEvent(new Event('change', { bubbles: true }));
};

// (A)의 한 가로줄에서 영교차 간격으로 줄무늬 주기를 잰다.
// 마루를 세는 방식은 λ가 화면보다 길면 실패하지만, 영교차는 반주기마다
// 하나씩 나오므로 화면에 두 개만 있어도 주기를 얻는다.
const measurePeriod = () => {
  const c = document.getElementById('canvasA');
  const y = Math.round(c.height * 0.15); // λ 눈금 막대·모서리 라벨·물체를 피한 줄
  const d = c.getContext('2d').getImageData(0, y, c.width, 1).data;
  const val = [];
  for (let i = 0; i < c.width; i++) {
    val.push(d[i * 4] - (d[i * 4 + 1] + d[i * 4 + 2]) / 2); // 빨강 +, 파랑 −
  }
  const zeros = [];
  for (let i = 1; i < c.width; i++) {
    if ((val[i - 1] < 0) !== (val[i] < 0)) {
      const t = Math.abs(val[i - 1]) / (Math.abs(val[i - 1]) + Math.abs(val[i]));
      zeros.push(i - 1 + t);
    }
  }
  if (zeros.length < 2) return null;
  const gaps = [];
  for (let i = 1; i < zeros.length; i++) gaps.push(zeros[i] - zeros[i - 1]);
  gaps.sort((a, b) => a - b);
  return 2 * gaps[Math.floor(gaps.length / 2)]; // 중앙값 × 2 = 한 주기
};

// 중앙 세로줄에서 회색 표식이 차지하는 픽셀 수
const measureObject = () => {
  const c = document.getElementById('canvasA');
  const d = c.getContext('2d').getImageData(Math.round(c.width / 2), 0, 1, c.height).data;
  let n = 0;
  for (let i = 0; i < c.height; i++) {
    const r = d[i * 4], g = d[i * 4 + 1], b = d[i * 4 + 2];
    if (Math.abs(r - g) < 12 && Math.abs(g - b) < 12 && r < 215 && r > 120) n++;
  }
  return { px: n, frac: n / c.height };
};

const readNumbers = () => {
  const t = document.getElementById('lambdaReadout').textContent;
  const pick = (label) => {
    const m = t.match(new RegExp(label + '\\s*=\\s*([^·]+)'));
    return m ? m[1].trim() : null;
  };
  return { loa: pick('λ/a'), Q: pick('Q'), G: pick('G'), ratio: pick('\\|E_산란\\|/\\|E_입사\\|') };
};

const CASES = [
  ['3.19', [0.1 * Math.PI, true]],
  ['5.01', [0.5, false]],
  ['40.09', [4.0, false]]
];

(async () => {
  const browser = await chromium.launch();
  const page = await browser.newPage({ viewport: { width: 1920, height: 1000 } });
  const errors = [];
  page.on('console', (m) => { if (m.type() === 'error') errors.push(m.text()); });
  page.on('pageerror', (e) => errors.push(e.message));

  await page.goto(PAGE, { waitUntil: 'load' });
  await page.waitForTimeout(3500);

  const results = [];
  let failed = 0;
  const check = (name, ok, detail) => {
    results.push({ 항목: name, 결과: ok ? 'PASS' : 'FAIL', 측정: detail });
    if (!ok) failed++;
  };

  // --- SCALE_MODE = 'a' -------------------------------------------------
  const a = {};
  for (const [tag, lam] of CASES) {
    await page.evaluate(setLambda, lam);
    await page.waitForFunction(() => document.getElementById('status').textContent === '');
    await page.evaluate(() => window.RGDTest.freeze(0));
    await page.waitForTimeout(250);
    a[tag] = {
      period: await page.evaluate(measurePeriod),
      obj: await page.evaluate(measureObject),
      nums: await page.evaluate(readNumbers)
    };
  }

  // 1. 줄무늬 주기 ∝ λ
  const base = a[CASES[0][0]].period / parseFloat(CASES[0][0]);
  for (const [tag] of CASES) {
    const expected = base * parseFloat(tag);
    const dev = (a[tag].period - expected) / expected;
    check(
      `줄무늬 주기 ∝ λ (λ/a=${tag})`,
      Math.abs(dev) <= TOL,
      `${a[tag].period.toFixed(1)}px / 기대 ${expected.toFixed(1)}px, 편차 ${(dev * 100).toFixed(2)}%`
    );
  }

  // 2. 물체 표식 크기
  for (const [tag] of CASES) {
    check(
      `물체 지름 ≤ 패널 20% (λ/a=${tag})`,
      a[tag].obj.frac <= 0.20,
      `${(a[tag].obj.frac * 100).toFixed(1)}% (${a[tag].obj.px}px)`
    );
  }
  const sizes = CASES.map(([t]) => a[t].obj.px);
  check('물체 크기가 λ에 무관', new Set(sizes).size === 1, sizes.join(' / ') + ' px');

  // --- SCALE_MODE = 'lambda' 복원 ---------------------------------------
  await page.selectOption('#rangeSelect', 'lambda');
  await page.waitForFunction(() => document.getElementById('status').textContent === '');
  const lam = {};
  const hashes = [];
  for (const [tag, l] of CASES) {
    await page.evaluate(setLambda, l);
    await page.waitForFunction(() => document.getElementById('status').textContent === '');
    await page.evaluate(() => window.RGDTest.freeze(0));
    await page.waitForTimeout(250);
    lam[tag] = { nums: await page.evaluate(readNumbers) };
    await page.uncheck('#arrowBox');
    await page.waitForTimeout(200);
    hashes.push(await page.evaluate(() => {
      const c = document.getElementById('canvasA');
      const x0 = Math.round(c.width * 0.08), w = Math.round(c.width * 0.20);
      const d = c.getContext('2d').getImageData(x0, 0, w, Math.round(c.height * 0.60)).data;
      let x = 0;
      for (let i = 0; i < d.length; i++) x = (x * 31 + d[i]) >>> 0;
      return x;
    }));
    await page.check('#arrowBox');
    await page.waitForTimeout(150);
  }

  // 3. 물리 판독값은 축척과 무관해야 한다
  for (const [tag] of CASES) {
    const A = a[tag].nums, L = lam[tag].nums;
    check(
      `물리 판독값이 축척과 무관 (λ/a=${tag})`,
      A.Q === L.Q && A.G === L.G && A.ratio === L.ratio,
      `Q ${A.Q}=${L.Q} · G ${A.G}=${L.G} · |E| ${A.ratio}=${L.ratio}`
    );
  }

  // 4. λ 모드의 원래 성질 — (A)가 세 λ에서 동일
  check(
    "λ 모드에서 (A)가 세 λ에 동일",
    hashes[0] === hashes[1] && hashes[1] === hashes[2],
    hashes.join(' / ')
  );

  // 5. 레이아웃·에러
  const layout = await page.evaluate(() => ({
    scroll: document.body.scrollHeight, inner: window.innerHeight
  }));
  check('세로 스크롤 없음', layout.scroll === layout.inner, `${layout.scroll} / ${layout.inner}`);
  check('콘솔 에러 없음', errors.length === 0, errors.length ? errors.join(' | ') : '0건');

  console.log('=== 화면 회귀 검증 ===\n');
  console.table(results);
  console.log(failed === 0
    ? `\nPASS: 전체 ${results.length}항목 통과`
    : `\nFAIL: ${failed}/${results.length}항목 실패`);

  await browser.close();
  process.exit(failed === 0 ? 0 : 1);
})();
