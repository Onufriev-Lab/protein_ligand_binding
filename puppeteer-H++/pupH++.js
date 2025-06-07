const puppeteer = require('puppeteer');
require('dotenv').config();
const fs = require('fs');
const path = require('path');
const fetch = require('node-fetch'); 

(async () => {
  const USERNAME = process.env.HPP_USERNAME;
  const PASSWORD = process.env.HPP_PASSWORD;
  const pdbCode = process.argv[2];

  if (!USERNAME || !PASSWORD || !pdbCode) {
    console.error('Usage: node screenshot_process_file.js <pdb_code>');
    process.exit(1);
  }

  // Convert to WSL-friendly path if needed (use actual WSL file mount path!)
  const wslFilePath = `/home/johann/vt_ml/pdb txt and LIGAND files/${pdbCode}_LIGAND.pdb`;

  if (!fs.existsSync(wslFilePath)) {
    console.error(`File not found: ${wslFilePath}`);
    process.exit(1);
  }

  const browser = await puppeteer.launch({
    headless: 'shell',
    args: ['--no-sandbox', '--disable-setuid-sandbox','--safebrowsing-disable-download-protection','--disable-extensions','--disable-popup-blocking']
  });

  const page = await browser.newPage();
  page.setDefaultTimeout(60000); // increase timeout safety

  // 1. Go to login page
  await page.goto('http://newbiophysics.cs.vt.edu/H++/', { waitUntil: 'networkidle2' });

  // 2. Login
  await page.type('input[name="username"]', USERNAME);
  await page.type('input[name="pass"]', PASSWORD);
  await Promise.all([
    page.click('input[name="submit"]'),
    page.waitForNavigation({ waitUntil: 'networkidle2' })
  ]);

  // 3. Click “Process Structure”
  await Promise.all([
    page.click('a[href="uploadpdb.php"]'),
    page.waitForNavigation({ waitUntil: 'networkidle2' })
  ]);

  // 4. Upload file
  const input = await page.$('input[type="file"][name="userfile"]');
  if (!input) {
    console.error('Upload input not found!');
    await browser.close();
    return;
  }

  await input.uploadFile(wslFilePath);

  // 5. Click “Process File” and wait for load or delay fallback
  await Promise.all([
    page.click('input[type="submit"][value="Process File"]')
  ]);

// Wait indefinitely until the "Process" button appears (lower right)
  await page.waitForSelector('input[type="image"][src="images/processthis.jpg"]', { timeout: 0 });

// Wait for the pH input field to exist
await page.waitForSelector('input[name="phlevel"]', { timeout: 0 });

// Set the pH value to 7.0
await page.evaluate(() => {
  const phInput = document.querySelector('input[name="phlevel"]');
  if (phInput) phInput.value = '7.0';
});

// Wait until the DOM reflects the updated pH value
await page.waitForFunction(
  () => {
    const el = document.querySelector('input[name="phlevel"]');
    return el && el.value === '7.0';
  },
  { timeout: 0 }
);

// Click the "Process" button and wait for processing page
await Promise.all([
  page.click('input[type="image"][src="images/processthis.jpg"]'),
  page.waitForNavigation({ waitUntil: 'domcontentloaded', timeout: 0 })
]);

// Wait for "VIEW RESULTS" link
await page.waitForSelector('a[href^="display_results.php"]', { timeout: 0 });

// Optional screenshot before clicking "VIEW RESULTS"
await page.screenshot({ path: `${pdbCode}_before_view_results.png`, fullPage: true });

// Click "VIEW RESULTS"
await Promise.all([
  page.click('a[href^="display_results.php"]'),
  page.waitForNavigation({ waitUntil: 'domcontentloaded', timeout: 0 })
]);

// Wait for the result PDB file link to confirm completion
await page.waitForSelector('a[href$=".result.pdb"]', { timeout: 0 });

// Final screenshot of completed results
await page.screenshot({ path: `${pdbCode}_final_results.png`, fullPage: true });

// Create download directory
const downloadPath = path.resolve(__dirname, 'downloads');
fs.mkdirSync(downloadPath, { recursive: true });

// Extract .top and .crd URLs from the page
const resultPrefix = `0.15_80_10_pH7.0_${pdbCode}_LIGAND`;
const fileLinks = await page.evaluate((prefix) => {
  const anchors = Array.from(document.querySelectorAll('a'));
  const result = { top: null, crd: null };

  for (const a of anchors) {
    const href = a.getAttribute('href');
    if (!href) continue;
    const parts = href.split('/');
    const fileName = parts[parts.length - 1];
    if (fileName === `${prefix}.top`) {
      result.top = new URL(href, window.location.origin).href;
    } else if (fileName === `${prefix}.crd`) {
      result.crd = new URL(href, window.location.origin).href;
    }
  }

  return result;
}, resultPrefix);

// Function to download a file using fetch and save locally
const downloadFile = async (url, filename) => {
  if (!url) {
    console.warn(`No URL found for ${filename}`);
    return;
  }
  const res = await fetch(url);
  if (!res.ok) {
    console.error(`Failed to fetch ${url}: ${res.statusText}`);
    return;
  }
  const fileStream = fs.createWriteStream(path.join(downloadPath, filename));
  await new Promise((resolve, reject) => {
    res.body.pipe(fileStream);
    res.body.on('error', reject);
    fileStream.on('finish', resolve);
  });
  console.log(`Downloaded: ${filename}`);
};


// Download the .top and .crd files
await downloadFile(fileLinks.top, `${pdbCode}_LIGAND.top`);
await downloadFile(fileLinks.crd, `${pdbCode}_LIGAND.crd`);


await browser.close();
})();
