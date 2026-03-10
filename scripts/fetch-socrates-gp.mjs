#!/usr/bin/env node
/**
 * Batch-download SOCRATES GP data (JSON) for top N conjunctions.
 * 
 * SOCRATES data.php returns the EXACT GP data used for the computation,
 * not the current catalog. This is critical for validation.
 *
 * Usage: node fetch-socrates-gp.mjs [--top N] [--rate MS]
 */

import { readFileSync, writeFileSync, existsSync, mkdirSync } from 'fs';
import { join, dirname } from 'path';
import { fileURLToPath } from 'url';

const __dirname = dirname(fileURLToPath(import.meta.url));
const DATA_DIR = join(__dirname, '..', 'tests', 'data', 'socrates_gp');

const args = process.argv.slice(2);
const getArg = (name, def) => { const i = args.indexOf(name); return i >= 0 && args[i+1] ? args[i+1] : def; };
const TOP = parseInt(getArg('--top', '1000'));
const RATE_MS = parseInt(getArg('--rate', '200')); // 5/sec

mkdirSync(DATA_DIR, { recursive: true });

// Parse SOCRATES CSV
// Use the most recent SOCRATES CSV
const csvPath = existsSync(join(__dirname, '..', 'tests', 'data', 'socrates_current.csv'))
  ? join(__dirname, '..', 'tests', 'data', 'socrates_current.csv')
  : join(__dirname, '..', 'tests', 'data', 'socrates_maxprob.csv');
console.log(`Using CSV: ${csvPath}`);
const csv = readFileSync(csvPath, 'utf8');
const lines = csv.trim().split('\n');
const header = lines[0].split(',');

const pairs = [];
const seen = new Set();

for (let i = 1; i < lines.length && pairs.length < TOP; i++) {
  const fields = [];
  let current = '', inQ = false;
  for (const ch of lines[i]) {
    if (ch === '"') { inQ = !inQ; continue; }
    if (ch === ',' && !inQ) { fields.push(current.trim()); current = ''; continue; }
    current += ch;
  }
  fields.push(current.trim());

  const id1 = fields[0], id2 = fields[3];
  const key = `${id1},${id2}`;
  if (!seen.has(key)) {
    seen.add(key);
    pairs.push({ id1, id2 });
  }
}

console.log(`Fetching GP data for ${pairs.length} SOCRATES pairs...`);
console.log(`Output: ${DATA_DIR}`);
console.log(`Rate: ${RATE_MS}ms between requests (${(1000/RATE_MS).toFixed(0)}/sec)\n`);

let downloaded = 0, cached = 0, errors = 0;

for (let i = 0; i < pairs.length; i++) {
  const { id1, id2 } = pairs[i];
  const outFile = join(DATA_DIR, `gp_${id1},${id2}.json`);

  if (existsSync(outFile)) {
    cached++;
    continue;
  }

  const url = `https://celestrak.org/SOCRATES/data.php?CATNR=${id1},${id2}&FORMAT=json`;
  try {
    const res = await fetch(url);
    if (!res.ok) throw new Error(`HTTP ${res.status}`);
    const text = await res.text();
    const data = JSON.parse(text);
    if (!data || data.length < 2) throw new Error(`Only ${data?.length || 0} objects`);
    writeFileSync(outFile, text);
    downloaded++;
  } catch (e) {
    errors++;
    if (errors <= 10) console.error(`  ✗ ${id1},${id2}: ${e.message}`);
  }

  if ((downloaded + errors) % 100 === 0 || i === pairs.length - 1) {
    console.log(`  [${i+1}/${pairs.length}] ${downloaded} downloaded, ${cached} cached, ${errors} errors`);
  }

  await new Promise(r => setTimeout(r, RATE_MS));
}

console.log(`\nDone: ${downloaded} downloaded, ${cached} cached, ${errors} errors`);
console.log(`Total GP files: ${downloaded + cached}`);
