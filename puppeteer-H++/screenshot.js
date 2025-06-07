// screenshot.js
const puppeteer = require('puppeteer');

(async () => {
  const url = 'http://newbiophysics.cs.vt.edu/H++';
  const browser = await puppeteer.launch({
    headless: 'shell',  // for Puppeteer v20+; use true for older versions
    args: ['--ignore-certificate-errors', '--unsafely-treat-insecure-origin-as-secure', '--disable-web-security']
  });

  const page = await browser.newPage();
  await page.goto(url, { waitUntil: 'networkidle2' });

  // Optional: set viewport for consistent screenshot size
  await page.setViewport({ width: 1280, height: 1024 });

  // Take the screenshot
  await page.screenshot({ path: 'hplusplus_homepage.png', fullPage: true });

  await browser.close();
})();
