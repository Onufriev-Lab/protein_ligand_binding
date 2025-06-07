const puppeteer = require('puppeteer');
require('dotenv').config();
const fs = require('fs');
const path = require('path');


(async () => {
  const USERNAME = process.env.HPP_USERNAME;
  const PASSWORD = process.env.HPP_PASSWORD;
  const pdbCode = process.argv[2];

  if (!USERNAME || !PASSWORD || !pdbCode) {
    console.error('Usage: node screenshot_process_file.js <pdb_code>');
    process.exit(1);
  }

  const wslFilePath = `/home/johann/vt_ml/pdb txt and LIGAND files/${pdbCode}_LIGAND.pdb`;

  if (!fs.existsSync(wslFilePath)) {
    console.error(`File not found: ${wslFilePath}`);
    process.exit(1);
  }

  const downloadPath = path.resolve(__dirname, 'downloads');
  fs.mkdirSync(downloadPath, { recursive: true });

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
  for await (const chunk of res.body) {
    fileStream.write(chunk);
  }
  fileStream.end();
  console.log(`Downloaded: ${filename}`);
};

  const browser = await puppeteer.launch({
    headless: 'shell',
    args: [
      '--no-sandbox',
      '--disable-setuid-sandbox',
      '--safebrowsing-disable-download-protection',
      '--disable-extensions',
      '--disable-popup-blocking'
    ]
  });

  const page = await browser.newPage();
  page.setDefaultTimeout(60000);

  await page.goto('http://newbiophysics.cs.vt.edu/H++/', { waitUntil: 'networkidle2' });

  await page.type('input[name="username"]', USERNAME);
  await page.type('input[name="pass"]', PASSWORD);
  await Promise.all([
    page.click('input[name="submit"]'),
    page.waitForNavigation({ waitUntil: 'networkidle2' })
  ]);

  await Promise.all([
    page.click('a[href="uploadpdb.php"]'),
    page.waitForNavigation({ waitUntil: 'networkidle2' })
  ]);

  const input = await page.$('input[type="file"][name="userfile"]');
  if (!input) {
    console.error('Upload input not found!');
    await browser.close();
    return;
  }

  await input.uploadFile(wslFilePath);

  await Promise.all([
    page.click('input[type="submit"][value="Process File"]')
  ]);

  await page.waitForSelector('input[type="image"][src="images/processthis.jpg"]', { timeout: 0 });
  await page.waitForSelector('input[name="phlevel"]', { timeout: 0 });

  await page.evaluate(() => {
    const phInput = document.querySelector('input[name="phlevel"]');
    if (phInput) phInput.value = '7.0';
  });

  await page.waitForFunction(
    () => {
      const el = document.querySelector('input[name="phlevel"]');
      return el && el.value === '7.0';
    },
    { timeout: 0 }
  );

  await Promise.all([
    page.click('input[type="image"][src="images/processthis.jpg"]'),
    page.waitForNavigation({ waitUntil: 'domcontentloaded', timeout: 0 })
  ]);

  await page.waitForSelector('a[href^="display_results.php"]', { timeout: 0 });

  await page.screenshot({ path: `${pdbCode}_before_view_results.png`, fullPage: true });

  await downloadFile(`http://newbiophysics.cs.vt.edu/H++/uploads/johann/0.15_80_10_pH7_${pdbCode}_LIGAND/0.15_80_10_pH7_${pdbCode}_LIGAND.top`, `${pdbCode}_LIGAND.top`);
  await downloadFile(`http://newbiophysics.cs.vt.edu/H++/uploads/johann/0.15_80_10_pH7_${pdbCode}_LIGAND/0.15_80_10_pH7_${pdbCode}_LIGAND.crd`, `${pdbCode}_LIGAND.crd`);

  await browser.close();
})();
