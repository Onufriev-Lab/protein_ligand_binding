const puppeteer = require('puppeteer');
require('dotenv').config();
const fs = require('fs');
const path = require('path');
const csv = require('csv-parser');

const USERNAME = process.env.HPP_USERNAME;
const PASSWORD = process.env.HPP_PASSWORD;
const csvFilePath = './pdb_codes.csv';
const topCrdOutputDir = path.resolve(__dirname, 'top-crd-files');
const summaryFile = path.join(topCrdOutputDir, 'hpp_summary.txt');

// Ensure output directory exists
fs.mkdirSync(topCrdOutputDir, { recursive: true });

// Read PDB codes from CSV
const readPdbCodes = () => {
  return new Promise((resolve, reject) => {
    const codes = [];
    fs.createReadStream(csvFilePath)
      .pipe(csv())
      .on('data', (row) => {
        if (row['PDB code']) codes.push(row['PDB code'].trim());
      })
      .on('end', () => resolve(codes))
      .on('error', reject);
  });
};

// Download helper
const downloadFile = async (url, filename) => {
  const res = await fetch(url);
  if (!res.ok) throw new Error(`Failed to fetch ${url}: ${res.statusText}`);
  const fileStream = fs.createWriteStream(path.join(topCrdOutputDir, filename));
  for await (const chunk of res.body) fileStream.write(chunk);
  fileStream.end();
};

(async () => {
  if (!USERNAME || !PASSWORD) {
    console.error('Please set HPP_USERNAME and HPP_PASSWORD in your .env file');
    process.exit(1);
  }

  const pdbCodes = await readPdbCodes();
  const successes = [];
  const failures = [];

  for (const pdbCode of pdbCodes) {
    const pdbPath = `/home/johann/protein_ligand_binding/divide-charge-remark-dev/Remarked-files/${pdbCode}_LIGAND_remarked.pdb`;

    if (!fs.existsSync(pdbPath)) {
      console.error(`[${pdbCode}] File not found. Skipping.`);
      failures.push(pdbCode);
      continue;
    }

    console.log(`[${pdbCode}] Processing...`);
    const browser = await puppeteer.launch({
      headless: 'shell',
      protocolTimeout: 0,
      args: ['--no-sandbox', '--disable-setuid-sandbox']
    });

    try {
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
      if (!input) throw new Error('Upload input not found');

      await input.uploadFile(pdbPath);
      await page.click('input[type="submit"][value="Process File"]');

      const result = await Promise.race([
        page.waitForSelector('input[name="phlevel"]', { timeout: 300000 }).then(() => 'success'),
        page.waitForFunction(() =>
          document.body && document.body.innerText.includes("THE CALCULATION HAS STOPPED"),
          { timeout: 300000 }
        ).then(() => 'halted')
      ]);

      if (result === 'halted') {
        console.error(`[${pdbCode}] H++ halted: THE CALCULATION HAS STOPPED`);
        failures.push(pdbCode);
        await browser.close();
        continue;
      }

      await page.evaluate(() => {
        const phInput = document.querySelector('input[name="phlevel"]');
        if (phInput) phInput.value = '7.0';
      });

      await page.waitForFunction(() =>
        document.querySelector('input[name="phlevel"]')?.value === '7.0'
      );

      await Promise.all([
        page.click('input[type="image"][src="images/processthis.jpg"]'),
        page.waitForNavigation({ waitUntil: 'domcontentloaded' })
      ]);

      const successSelector = 'a[href^="display_results.php"]';
      const failureSelector = 'a[href*="debug.php"]';

      const finalResult = await Promise.race([
        page.waitForSelector(successSelector).then(() => 'success'),
        page.waitForSelector(failureSelector).then(() => 'failure')
      ]);

      if (finalResult === 'failure') {
        console.error(`[${pdbCode}] Processing failed.`);
        failures.push(pdbCode);
        await browser.close();
        continue;
      }

      await Promise.all([
        page.click(successSelector),
        page.waitForNavigation({ waitUntil: 'domcontentloaded' })
      ]);

      await new Promise(res => setTimeout(res, 10000)); // 10 second wait

      const topUrl = `http://newbiophysics.cs.vt.edu/H++/uploads/johann/0.15_80_10_pH7_${pdbCode}_LIGAND/0.15_80_10_pH7_${pdbCode}_LIGAND.top`;
      const crdUrl = `http://newbiophysics.cs.vt.edu/H++/uploads/johann/0.15_80_10_pH7_${pdbCode}_LIGAND/0.15_80_10_pH7_${pdbCode}_LIGAND.crd`;

      await downloadFile(topUrl, `${pdbCode}.top`);
      await downloadFile(crdUrl, `${pdbCode}.crd`);

      successes.push(pdbCode);
      console.log(`[${pdbCode}] Success.`);
    } catch (err) {
      console.error(`[${pdbCode}] Error: ${err.message}`);
      failures.push(pdbCode);
    } finally {
      await browser.close();
    }
  }

  const summary = [
    `Total PDB codes: ${pdbCodes.length}`,
    `Successful: ${successes.length}`,
    `Failed: ${failures.length}`,
    '',
    'Failed PDB codes:',
    ...failures
  ].join('\n');

  console.log('\n' + summary);
  fs.writeFileSync(summaryFile, summary);
})();
