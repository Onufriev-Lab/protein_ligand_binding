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

  const wslFilePath = `../divide-charge-remark-dev/Remarked-files/${pdbCode}_LIGAND_remarked.pdb`;

  if (!fs.existsSync(wslFilePath)) {
    console.error(`File not found: ${wslFilePath}`);
    process.exit(1);
  }

  const downloadPath = path.resolve(__dirname, 'top-crd-files');
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
  console.log(fileStream);
  for await (const chunk of res.body) {
    fileStream.write(chunk);
  }
  fileStream.end();
  console.log(`Downloaded: ${filename}`);
};

  const browser = await puppeteer.launch({
    headless: 'shell',
    protocolTimeout: 0,
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

await page.click('input[type="submit"][value="Process File"]');

console.log('Waiting for either validation form or early failure...');

// Wait for either the phlevel field (success) or failure message (stop)
const result = await Promise.race([
  page.waitForSelector('input[name="phlevel"]', { timeout: 300000 }).then(() => 'phlevel'),
  page.waitForFunction(
    () => document.body && document.body.innerText.includes("THE CALCULATION HAS STOPPED"),
    { timeout: 300000 }
  ).then(() => 'calculation_stopped')
]);


if (result === 'calculation_stopped') {
  console.error(`H++ halted before processing for ${pdbCode}. Reason: THE CALCULATION HAS STOPPED.`);
  await browser.close();
  process.exit(1);
}

console.log('Validation form loaded — continuing with processing...');


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

  console.log('Process this pressed')

try {
  // Wait indefinitely for either the success or failure link to appear
  const successSelector = 'a[href^="display_results.php"]';
  const failureSelector = 'a[href*="debug.php"]';

  const foundSelector = await Promise.race([
    page.waitForSelector(successSelector, { timeout: 0 }).then(() => successSelector),
    page.waitForSelector(failureSelector, { timeout: 0 }).then(() => failureSelector)
  ]);

  if (foundSelector === successSelector) {
    await Promise.all([
      page.click(successSelector),
      console.log('Success clicked'),
      page.waitForNavigation({ waitUntil: 'domcontentloaded', timeout: 0 })
    ]);
  } else if (foundSelector === failureSelector) {
    console.error(`H++ processing failed for ${pdbCode}. "Files generated so far" link detected.`);
    await browser.close();
    process.exit(1);
  }
} catch (err) {
  console.error(`Error waiting for result or failure link: ${err}`);
  await browser.close();
  process.exit(1);
}
console.log('Before setTimeout');
await new Promise(resolve => setTimeout(resolve, 10000));
console.log('After setTimeout');

  await downloadFile(`http://newbiophysics.cs.vt.edu/H++/uploads/johann/0.15_80_10_pH7.0_${pdbCode}_LIGAND_remarked/0.15_80_10_pH7.0_${pdbCode}_LIGAND_remarked.top`, `${pdbCode}.top`);
  await downloadFile(`http://newbiophysics.cs.vt.edu/H++/uploads/johann/0.15_80_10_pH7.0_${pdbCode}_LIGAND_remarked/0.15_80_10_pH7.0_${pdbCode}_LIGAND_remarked.crd`, `${pdbCode}.crd`);

  await browser.close();
  console.log('Just closed the browser');
})();
