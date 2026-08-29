// Headless render test: verify every <img> actually renders (换机可见 test)
const { chromium } = require('/Users/luoxiaowen/.workbuddy/binaries/node/workspace/node_modules/playwright');
const path = require('path');
const fs = require('fs');

const DIR = '/Users/luoxiaowen/Desktop/LKDock/LKina论文';
const TMP = '/tmp/lkina_render_test';
fs.mkdirSync(TMP, { recursive: true });

(async () => {
  const browser = await chromium.launch({
    executablePath: '/Users/luoxiaowen/Library/Caches/ms-playwright/chromium_headless_shell-1217/chrome-headless-shell-mac-arm64/chrome-headless-shell',
  });
  const files = ['LKina论文.html', 'LKina论文_中文.html', 'LKina_paper_EN.html'];
  let allOk = true;

  for (const f of files) {
    const page = await browser.newPage({ viewport: { width: 1000, height: 800 } });
    const errors = [];
    page.on('pageerror', e => errors.push('pageerror: ' + e.message));
    page.on('console', m => { if (m.type() === 'error') errors.push('console: ' + m.text()); });

    await page.goto('file://' + path.join(DIR, f), { waitUntil: 'networkidle' });
    await page.waitForTimeout(800);

    const stats = await page.evaluate(() => {
      const imgs = [...document.querySelectorAll('img')];
      const loaded = imgs.filter(i => i.complete && i.naturalWidth > 0);
      const broken = imgs.filter(i => !(i.complete && i.naturalWidth > 0));
      return {
        total: imgs.length,
        loaded: loaded.length,
        broken: broken.map(i => ({ src: i.src.slice(0, 40), w: i.naturalWidth })),
        figures: document.querySelectorAll('figure').length,
      };
    });

    // scroll through whole page to force lazy/any deferred rendering, re-check
    await page.evaluate(async () => {
      for (let y = 0; y < document.body.scrollHeight; y += 600) {
        window.scrollTo(0, y);
        await new Promise(r => setTimeout(r, 60));
      }
    });
    const after = await page.evaluate(() => {
      const imgs = [...document.querySelectorAll('img')];
      const broken = imgs.filter(i => !(i.complete && i.naturalWidth > 0));
      return broken.length;
    });

    await page.screenshot({ path: path.join(TMP, f.replace('.html', '_full.png')), fullPage: true });

    const ok = stats.total === stats.loaded && after === 0 && errors.length === 0;
    if (!ok) allOk = false;
    console.log(`== ${f} ==`);
    console.log(`  imgs: ${stats.loaded}/${stats.total} loaded | broken-after-scroll: ${after} | figures: ${stats.figures} | js-errors: ${errors.length}`);
    errors.slice(0, 3).forEach(e => console.log('  ERR:', e));
    console.log(ok ? '  RENDER: PASS' : '  RENDER: FAIL');
    await page.close();
  }
  await browser.close();
  console.log(allOk ? '\n=== ALL PASS ===' : '\n=== HAS FAILURES ===');
  process.exit(allOk ? 0 : 1);
})();
