#!/usr/bin/env node
/**
 * SOCRATES CSV Scraper
 *
 * Downloads SOCRATES CSV data from CelesTrak and fetches GP/TLE data
 * for each conjunction pair. Outputs test fixtures.
 *
 * Usage:
 *   node scrape-socrates.mjs [--top N] [--sort minRange|maxProb]
 *   node scrape-socrates.mjs --fetch-gp tests/data/socrates_maxprob.csv --top 20
 *
 * CSV endpoints (full catalog, ~120K conjunctions):
 *   https://celestrak.org/SOCRATES/sort-minRange.csv
 *   https://celestrak.org/SOCRATES/sort-maxProb.csv
 *
 * CSV format (RFC 4180):
 *   NORAD_CAT_ID_1,OBJECT_NAME_1,DSE_1,NORAD_CAT_ID_2,OBJECT_NAME_2,DSE_2,
 *   TCA,TCA_RANGE,TCA_RELATIVE_SPEED,MAX_PROB,DILUTION
 *
 * GP data for each pair:
 *   https://celestrak.org/SOCRATES/data.php?CATNR={id1},{id2}
 */

import { writeFileSync, readFileSync, mkdirSync, existsSync } from 'fs';
import { join, dirname } from 'path';
import { fileURLToPath } from 'url';

const __dirname = dirname(fileURLToPath(import.meta.url));
const DATA_DIR = join(__dirname, '..', 'tests', 'data');

const BASE_URL = 'https://celestrak.org/SOCRATES';

// Parse args
const args = process.argv.slice(2);
const getArg = (name, def) => {
  const i = args.indexOf(name);
  return i >= 0 && args[i + 1] ? args[i + 1] : def;
};

const TOP = parseInt(getArg('--top', '10'));
const SORT = getArg('--sort', 'maxProb');
const FETCH_GP = getArg('--fetch-gp', null);

function parseCSV(text) {
  const lines = text.trim().split('\n');
  const header = lines[0].split(',');
  return lines.slice(1).map(line => {
    // Handle quoted fields (object names can contain commas in brackets)
    const fields = [];
    let current = '';
    let inQuotes = false;
    for (const ch of line) {
      if (ch === '"') { inQuotes = !inQuotes; continue; }
      if (ch === ',' && !inQuotes) { fields.push(current.trim()); current = ''; continue; }
      current += ch;
    }
    fields.push(current.trim());

    const obj = {};
    header.forEach((h, i) => obj[h] = fields[i] || '');
    return obj;
  });
}

async function fetchCSV() {
  const sortMap = { maxProb: 'sort-maxProb.csv', minRange: 'sort-minRange.csv' };
  const csvFile = sortMap[SORT] || sortMap.maxProb;
  const url = `${BASE_URL}/${csvFile}`;

  console.log(`Fetching ${url}...`);
  const res = await fetch(url);
  if (!res.ok) throw new Error(`HTTP ${res.status}`);
  const text = await res.text();

  mkdirSync(DATA_DIR, { recursive: true });
  const outPath = join(DATA_DIR, `socrates_${SORT}.csv`);
  writeFileSync(outPath, text);
  console.log(`Saved ${outPath} (${text.split('\n').length} rows)`);
  return text;
}

async function fetchGPData(conjunctions) {
  console.log(`\nFetching GP data for ${conjunctions.length} conjunctions...`);
  mkdirSync(DATA_DIR, { recursive: true });

  for (let i = 0; i < conjunctions.length; i++) {
    const c = conjunctions[i];
    const id1 = c.NORAD_CAT_ID_1;
    const id2 = c.NORAD_CAT_ID_2;

    // Fetch in all 3 formats: TLE (txt), JSON, CSV
    const formats = [
      { ext: 'txt', param: '' },           // Default TLE text
      { ext: 'json', param: '&FORMAT=json' },
      { ext: 'csv', param: '&FORMAT=csv' },
    ];

    for (const fmt of formats) {
      const gpFile = `gp_${id1},${id2}.${fmt.ext}`;
      const gpPath = join(DATA_DIR, gpFile);

      if (existsSync(gpPath)) {
        continue; // cached
      }

      const gpUrl = `${BASE_URL}/data.php?CATNR=${id1},${id2}${fmt.param}`;
      try {
        const res = await fetch(gpUrl);
        if (!res.ok) throw new Error(`HTTP ${res.status}`);
        const gpText = await res.text();
        writeFileSync(gpPath, gpText);
      } catch (e) {
        console.error(`  [${i + 1}] ${id1},${id2}.${fmt.ext} — ✗ ${e.message}`);
      }
      await new Promise(r => setTimeout(r, 300));
    }

    console.log(`  [${i + 1}/${conjunctions.length}] ${id1},${id2} — ✓ (txt+json+csv)`);
  }
}

async function main() {
  let csvText;

  if (FETCH_GP) {
    // Load existing CSV and fetch GP data
    csvText = readFileSync(FETCH_GP, 'utf8');
  } else {
    // Download fresh CSV
    csvText = await fetchCSV();
  }

  const all = parseCSV(csvText);
  console.log(`Parsed ${all.length} total conjunctions`);

  const top = all.slice(0, TOP);
  console.log(`Processing top ${top.length}`);

  // Fetch GP data for each
  await fetchGPData(top);

  // Write JSON fixture
  const fixture = {
    source: 'CelesTrak SOCRATES Plus (CSV)',
    scraped_at: new Date().toISOString(),
    sort: SORT,
    total: all.length,
    top_n: TOP,
    conjunctions: top.map(c => ({
      obj1_norad: parseInt(c.NORAD_CAT_ID_1),
      obj1_name: c.OBJECT_NAME_1,
      obj1_dse: parseFloat(c.DSE_1),
      obj2_norad: parseInt(c.NORAD_CAT_ID_2),
      obj2_name: c.OBJECT_NAME_2,
      obj2_dse: parseFloat(c.DSE_2),
      tca: c.TCA.replace(' ', 'T') + 'Z',
      min_range_km: parseFloat(c.TCA_RANGE),
      rel_speed_kms: parseFloat(c.TCA_RELATIVE_SPEED),
      max_prob: parseFloat(c.MAX_PROB),
      dilution_km: parseFloat(c.DILUTION),
      gp_file: `gp_${c.NORAD_CAT_ID_1},${c.NORAD_CAT_ID_2}.txt`
    }))
  };

  const fixturePath = join(DATA_DIR, `socrates_top${TOP}_${SORT}.json`);
  writeFileSync(fixturePath, JSON.stringify(fixture, null, 2));
  console.log(`\nFixture: ${fixturePath}`);
}

main().catch(e => { console.error(e); process.exit(1); });
