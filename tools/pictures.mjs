// The website's pictures: the sections of the example report, for the report
// page, and the card a link to the site shows where it is shared.
//
// Opens docs/examples/demo_report.html in headless Chrome, scrolls to each
// top-level section and keeps the top of its column as a small WebP in
// docs/assets/report/. Then draws the card, 1200 by 630, from the front
// page's own drawing, into docs/assets/social-card.png, the icon a phone
// keeps on its home screen into docs/assets/touch-icon.png, and the README's
// pictures, light and dark, into docs/assets/readme/: the metro map and the
// five outcomes, as the site draws them. Run it after the example report, a
// drawing or the logo changes:
//
//   node tools/pictures.mjs
//
// Needs Node 22 or newer (for its WebSocket) and Google Chrome; CHROME points
// at another Chrome or Chromium binary.
import { spawn } from 'node:child_process'
import { mkdirSync, mkdtempSync, rmSync, writeFileSync } from 'node:fs'
import { tmpdir } from 'node:os'
import { dirname, join, resolve } from 'node:path'
import { fileURLToPath, pathToFileURL } from 'node:url'

const ROOT = resolve(dirname(fileURLToPath(import.meta.url)), '..')
const REPORT = join(ROOT, 'docs', 'examples', 'demo_report.html')
const OUT = join(ROOT, 'docs', 'assets', 'report')
// The report's top-level sections, in its order (the ids pandoc gives them).
export const SECTIONS = ['overview', 'amap-annotation-results', 'functional-annotation', 'gene-duplications', 'pangenome-analysis', 'amap-cds-stats', 'master-table-explorer']
const CHROME = process.env.CHROME || '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome'
const sleep = (ms) => new Promise((r) => setTimeout(r, ms))

const profile = mkdtempSync(join(tmpdir(), 'amap-shots-'))
const chrome = spawn(CHROME, ['--headless=new', '--remote-debugging-port=0', `--user-data-dir=${profile}`, '--hide-scrollbars',
  '--no-first-run', '--allow-file-access-from-files', '--window-size=1200,900', 'about:blank'], { stdio: ['ignore', 'ignore', 'pipe'] })
const wsUrl = await new Promise((ok, fail) => {
  let seen = ''
  const timer = setTimeout(() => fail(new Error('Chrome did not start')), 20000)
  chrome.stderr.on('data', (d) => {
    seen += d
    const m = seen.match(/DevTools listening on (ws:\/\/\S+)/)
    if (m) { clearTimeout(timer); ok(m[1]) }
  })
})
const ws = new WebSocket(wsUrl)
await new Promise((ok, fail) => { ws.onopen = ok; ws.onerror = fail })
let id = 0
const waiting = new Map()
ws.onmessage = (e) => {
  const msg = JSON.parse(e.data)
  if (msg.id && waiting.has(msg.id)) {
    const { ok, fail } = waiting.get(msg.id)
    waiting.delete(msg.id)
    msg.error ? fail(new Error(msg.error.message)) : ok(msg.result)
  }
}
const send = (method, params = {}, sessionId) => new Promise((ok, fail) => {
  waiting.set(++id, { ok, fail })
  ws.send(JSON.stringify({ id, method, params, ...(sessionId ? { sessionId } : {}) }))
})

try {
  const { targetId } = await send('Target.createTarget', { url: 'about:blank' })
  const { sessionId } = await send('Target.attachToTarget', { targetId, flatten: true })
  const page = (method, params) => send(method, params, sessionId)
  await page('Emulation.setDeviceMetricsOverride', { width: 1200, height: 900, deviceScaleFactor: 1, mobile: false })
  await page('Page.navigate', { url: pathToFileURL(REPORT).href })
  await sleep(5000)
  // The report's header (two logos on a white band, body::before) stays at the
  // top of the window; it would cover the title of every section, so it is
  // hidden for the pictures.
  await page('Runtime.evaluate', {
    expression: `document.querySelectorAll('body *').forEach((el) => {
      const cs = getComputedStyle(el)
      if ((cs.position === 'fixed' || cs.position === 'sticky') && el.getBoundingClientRect().top < 20) el.style.visibility = 'hidden'
    })
    document.head.insertAdjacentHTML('beforeend', '<style>body::before { display: none !important; }</style>')`,
  })
  mkdirSync(OUT, { recursive: true })
  for (const section of SECTIONS) {
    const { result } = await page('Runtime.evaluate', {
      returnByValue: true,
      expression: `(() => {
        const el = document.getElementById(${JSON.stringify(section)})
        if (!el) return null
        el.scrollIntoView({ block: 'start' })
        const r = el.getBoundingClientRect()
        return { x: r.left + scrollX, y: r.top + scrollY, w: r.width }
      })()`,
    })
    if (!result.value) throw new Error(`the example report has no section #${section}`)
    await sleep(1500)
    const box = result.value
    const { data } = await page('Page.captureScreenshot', {
      format: 'webp', quality: 72, captureBeyondViewport: true,
      clip: { x: box.x - 12, y: box.y - 12, width: box.w + 24, height: 400, scale: 0.8 },
    })
    writeFileSync(join(OUT, `${section}.webp`), Buffer.from(data, 'base64'))
    console.log(`docs/assets/report/${section}.webp`)
  }

  // The card and the icon: small pages of their own, on the site's stylesheet
  // and drawings, written next to the site so its paths resolve.
  const DOCS = join(ROOT, 'docs')
  const shoot = async (name, html, width, height, file) => {
    const tmp = join(DOCS, `.${name}.html`)
    writeFileSync(tmp, html)
    try {
      await page('Emulation.setDeviceMetricsOverride', { width, height, deviceScaleFactor: 1, mobile: false })
      await page('Page.navigate', { url: pathToFileURL(tmp).href })
      await sleep(2500)
      const { data } = await page('Page.captureScreenshot', { format: 'png', clip: { x: 0, y: 0, width, height, scale: 1 } })
      writeFileSync(join(DOCS, 'assets', file), Buffer.from(data, 'base64'))
      console.log(`docs/assets/${file}`)
    } finally {
      rmSync(tmp, { force: true })
    }
  }
  await shoot('social-card', `<!doctype html><html lang="en" data-theme="dark"><head><meta charset="utf-8">
<link rel="stylesheet" href="assets/site.css"><style>
body { margin: 0; width: 1200px; height: 630px; overflow: hidden;
  background: radial-gradient(900px 420px at 85% 0%, color-mix(in srgb, var(--accent) 16%, transparent), transparent 70%), var(--bg); }
.card { box-sizing: border-box; width: 1200px; height: 630px; padding: 50px 64px 40px; display: flex; flex-direction: column; }
.c-top { display: flex; align-items: center; gap: 16px; font-size: 30px; font-weight: 750; }
.c-top span { color: var(--muted); font-weight: 500; font-size: 22px; }
h1 { font-size: 62px; line-height: 1.04; letter-spacing: -0.025em; margin: 26px 0 0; max-width: 1000px; }
.draw { margin-top: auto; background: var(--bg-elev); border: 1px solid var(--line); border-radius: 18px; padding: 16px 20px 6px; }
.c-foot { display: flex; justify-content: space-between; color: var(--muted); font-size: 19px; margin-top: 16px; }
.c-foot b { color: var(--ink); font-weight: 650; }
</style></head><body><div class="card">
<div class="c-top"><img src="assets/amap-logo.png" width="86" height="60" alt="">AMAP<span>Assisted Microbial Annotator Pipeline</span></div>
<h1>Bakta calls the genes. <span class="grad">AMAP checks every one.</span></h1>
<div class="draw"><div class="demo-stage" data-card></div></div>
<div class="c-foot"><span>Real genes: <b><i>M. tuberculosis</i> CDC1551</b>, corrected against H37Rv</span><span>pathogenomics-lab.github.io/AMAP</span></div>
</div><script src="assets/figure-data.js"></script><script src="assets/figures.js"></script><script>
AMAPFigures.render(document.querySelector('[data-card]'), function () { return AMAPFigures.windowSpec('hero') }, function (svg) { svg.setAttribute('data-phase', '4') })
</script></body></html>`, 1200, 630, 'social-card.png')
  // The README's pictures: drawings of the site, in either theme, at twice
  // the size for sharp screens.
  const README = join(DOCS, 'assets', 'readme')
  mkdirSync(README, { recursive: true })
  const element = async (url, selector, file, theme, width) => {
    await page('Emulation.setDeviceMetricsOverride', { width, height: 900, deviceScaleFactor: 2, mobile: false })
    await page('Page.navigate', { url: pathToFileURL(join(DOCS, url)).href })
    await sleep(2500)
    // The theme is all in the stylesheet, so it can change once the page is in.
    const { result } = await page('Runtime.evaluate', {
      returnByValue: true,
      expression: `(() => { document.documentElement.dataset.theme = '${theme}'; const el = document.querySelector(${JSON.stringify(selector)}); el.scrollIntoView({ block: 'center' }); const r = el.getBoundingClientRect(); return { x: r.left + scrollX, y: r.top + scrollY, w: r.width, h: r.height } })()`,
    })
    await sleep(800)
    const b = result.value
    const { data } = await page('Page.captureScreenshot', { format: 'png', captureBeyondViewport: true, clip: { x: b.x, y: b.y, width: b.w, height: b.h, scale: 1 } })
    writeFileSync(join(README, file), Buffer.from(data, 'base64'))
    console.log(`docs/assets/readme/${file}`)
  }
  // The five outcomes, as the overview page draws them, laid out wider.
  const overview = (await import('node:fs')).readFileSync(join(DOCS, 'how-it-works', 'index.html'), 'utf8')
  const gallery = overview.match(/<figure class="fig outcomes"[\s\S]*?<\/figure>/)[0].replace(/<figcaption>[\s\S]*?<\/figcaption>/, '')
  const wide = join(DOCS, 'how-it-works', '.readme-outcomes.html')
  writeFileSync(wide, `<!doctype html><html lang="en" data-theme="dark"><head><meta charset="utf-8">
<link rel="stylesheet" href="../assets/site.css"><style>
body { margin: 0; padding: 16px; width: 1100px; }
.prose { max-width: none; }
.prose .fig { margin: 0; }
.oc-grid { grid-template-columns: repeat(3, minmax(0, 1fr)) !important; }
/* A picture has no links: nothing in it should look like one. */
.oc-more { display: none !important; }
.prose .oc-harder a { color: var(--ink-soft); text-decoration: none; }
</style></head><body><main class="prose">${gallery}</main>
<script src="../assets/figure-data.js"></script><script src="../assets/figures.js"></script></body></html>`)
  try {
    for (const theme of ['light', 'dark']) {
      await element('index.html', '#pipeline [data-metro-stage]', `metro-${theme}.png`, theme, 1240)
      await element('how-it-works/.readme-outcomes.html', '.fig.outcomes', `outcomes-${theme}.png`, theme, 1140)
    }
  } finally {
    rmSync(wide, { force: true })
  }

  await shoot('touch-icon', `<!doctype html><html lang="en" data-theme="dark"><head><meta charset="utf-8">
<link rel="stylesheet" href="assets/site.css"><style>
body { margin: 0; width: 180px; height: 180px; display: grid; place-items: center; background: var(--bg); }
img { width: 140px; height: auto; }
</style></head><body><img src="assets/amap-logo.png" alt=""></body></html>`, 180, 180, 'touch-icon.png')
} finally {
  // Chrome writes to its profile until it has exited, so wait for that.
  const exited = new Promise((ok) => chrome.once('exit', ok))
  try { await send('Browser.close') } catch {}
  ws.close()
  await Promise.race([exited, sleep(5000)])
  chrome.kill()
  rmSync(profile, { recursive: true, force: true, maxRetries: 5, retryDelay: 200 })
}
