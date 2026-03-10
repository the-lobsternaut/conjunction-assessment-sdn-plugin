#!/usr/bin/env node
/**
 * SOCRATES Scraper
 *
 * Scrapes CelesTrak SOCRATES Plus for conjunction data + GP/TLE data.
 * Outputs JSON test fixtures for validating the conjunction assessment plugin.
 *
 * Usage:
 *   node scrape-socrates.mjs [--max 10] [--order maxProb|minRange|relSpeed]
 */

import { writeFileSync, mkdirSync } from 'fs';
import { join, dirname } from 'path';
import { fileURLToPath } from 'url';

const __dirname = dirname(fileURLToPath(import.meta.url));
const DATA_DIR = join(__dirname, '..', 'tests', 'data');

const BASE_URL = 'https://celestrak.org/SOCRATES';
const MAX = parseInt(process.argv.find((a, i) => process.argv[i-1] === '--max') || '10');
const ORDER = process.argv.find((a, i) => process.argv[i-1] === '--order') || 'MAXPROB';

async function fetchText(url) {
  const res = await fetch(url);
  if (!res.ok) throw new Error(`HTTP ${res.status}: ${url}`);
  return res.text();
}

async function scrapeSOCRATES() {
  console.log(`Scraping SOCRATES Plus: top ${MAX} by ${ORDER}...`);

  // Fetch the HTML table page
  const tableUrl = `${BASE_URL}/table-socrates.php?NAME=,&ORDER=${ORDER}&MAX=${MAX}`;
  const html = await fetchText(tableUrl);

  // Parse metadata from HTML
  const dataCurrentMatch = html.match(/Data current as of ([\d\- :A-Z]+)/);
  const intervalMatch = html.match(/Start = ([\d\- :A-Z]+), Stop = ([\d\- :A-Z]+)/);
  const conjMatch = html.match(/([\d,]+) Conjunctions/);
  const primMatch = html.match(/([\d,]+) Primaries/);
  const secMatch = html.match(/([\d,]+) Secondaries/);

  // Parse conjunction rows using regex on table rows
  // Each conjunction is 2 table rows
  // Row 1: GP Data | NORAD1 | NAME1 | DSE1 | TCA | MIN_RANGE | REL_SPEED
  // Row 2: graphs  | NORAD2 | NAME2 | DSE2 | MAX_PROB | DILUTION
  const gpPattern = /data\.php\?CATNR=([\d,]+)/g;
  const catnrPairs = [...html.matchAll(gpPattern)].map(m => m[1]);

  // Parse table cells more carefully
  const rowPattern = /<tr[^>]*>[\s\S]*?<\/tr>/gi;
  const cellPattern = /<td[^>]*>([\s\S]*?)<\/td>/gi;
  const rows = [...html.matchAll(rowPattern)].map(m => m[0]);

  const conjunctions = [];
  let pairIndex = 0;

  for (let i = 0; i < rows.length - 1; i++) {
    const cells1 = [...rows[i].matchAll(cellPattern)].map(m => m[1].replace(/<[^>]+>/g, '').trim());
    const cells2 = [...rows[i + 1].matchAll(cellPattern)].map(m => m[1].replace(/<[^>]+>/g, '').trim());

    // Check if this is a data row (has GP Data)
    if (!rows[i].includes('data.php?CATNR=')) continue;

    const catnrs = catnrPairs[pairIndex++];
    if (!catnrs) break;

    // Find the actual data cells (skip the GP Data cell)
    const norad1 = cells1.find(c => /^\d{5,6}$/.test(c));
    const norad2 = cells2.find(c => /^\d{5,6}$/.test(c));
    if (!norad1 || !norad2) continue;

    const idx1 = cells1.indexOf(norad1);
    const idx2 = cells2.indexOf(norad2);

    const name1 = cells1[idx1 + 1] || '';
    const dse1 = parseFloat(cells1[idx1 + 2]) || 0;
    const tca = cells1[idx1 + 3] || '';
    const minRange = parseFloat(cells1[idx1 + 4]) || 0;
    const relSpeed = parseFloat(cells1[idx1 + 5]) || 0;

    const name2 = cells2[idx2 + 1] || '';
    const dse2 = parseFloat(cells2[idx2 + 2]) || 0;
    const maxProb = cells2[idx2 + 3] || '';
    const dilution = parseFloat(cells2[idx2 + 4]) || 0;

    // Parse ops status from name
    const parseStatus = (name) => {
      const m = name.match(/\[([^\]]+)\]/);
      return m ? m[1] : '';
    };

    conjunctions.push({
      obj1_norad: parseInt(norad1),
      obj1_name: name1.replace(/\s*\[[^\]]+\]/, '').trim(),
      obj1_status: parseStatus(name1),
      obj1_dse: dse1,
      obj2_norad: parseInt(norad2),
      obj2_name: name2.replace(/\s*\[[^\]]+\]/, '').trim(),
      obj2_status: parseStatus(name2),
      obj2_dse: dse2,
      tca: tca.replace(' ', 'T') + 'Z',
      min_range_km: minRange,
      rel_speed_kms: relSpeed,
      max_prob: parseFloat(maxProb) || 0,
      dilution_km: dilution,
      gp_file: `gp_${catnrs}.txt`
    });

    i++; // Skip the second row
  }

  // Fetch GP data for each conjunction
  console.log(`Found ${conjunctions.length} conjunctions. Fetching GP data...`);
  mkdirSync(DATA_DIR, { recursive: true });

  for (const conj of conjunctions) {
    const catnrs = `${conj.obj1_norad},${conj.obj2_norad}`;
    const gpUrl = `${BASE_URL}/data.php?CATNR=${catnrs}`;
    try {
      const gpData = await fetchText(gpUrl);
      const gpPath = join(DATA_DIR, conj.gp_file);
      writeFileSync(gpPath, gpData);
      console.log(`  ✓ ${catnrs} (${gpData.split('\n').length} lines)`);
    } catch (e) {
      console.error(`  ✗ ${catnrs}: ${e.message}`);
    }
    // Be polite to CelesTrak
    await new Promise(r => setTimeout(r, 500));
  }

  // Build output
  const now = new Date().toISOString();
  const output = {
    source: 'CelesTrak SOCRATES Plus',
    scraped_at: now,
    data_current_as_of: dataCurrentMatch?.[1] || '',
    computation_interval: {
      start: intervalMatch?.[1] || '',
      stop: intervalMatch?.[2] || ''
    },
    threshold_km: 5.0,
    primaries: parseInt((primMatch?.[1] || '0').replace(/,/g, '')),
    secondaries: parseInt((secMatch?.[1] || '0').replace(/,/g, '')),
    total_conjunctions: parseInt((conjMatch?.[1] || '0').replace(/,/g, '')),
    order: ORDER,
    max_results: MAX,
    conjunctions
  };

  const outPath = join(DATA_DIR, `socrates_top${MAX}_${ORDER.toLowerCase()}_${now.slice(0, 10)}.json`);
  writeFileSync(outPath, JSON.stringify(output, null, 2));
  console.log(`\nWritten: ${outPath}`);
  console.log(`${conjunctions.length} conjunctions with GP data saved.`);
}

scrapeSOCRATES().catch(e => {
  console.error('Scraper failed:', e);
  process.exit(1);
});
