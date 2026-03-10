#!/usr/bin/env node
/**
 * Batch-download SOCRATES GP data via Tor with circuit rotation.
 * Rotates Tor circuit on HTTP 429 / 500 / timeout.
 *
 * Usage: node fetch-socrates-gp-tor.mjs [--top N] [--rate MS] [--start OFFSET]
 */

import { readFileSync, writeFileSync, existsSync, mkdirSync } from 'fs';
import { join, dirname } from 'path';
import { fileURLToPath } from 'url';
import net from 'net';
import { SocksProxyAgent } from 'socks-proxy-agent';

const __dirname = dirname(fileURLToPath(import.meta.url));
const DATA_DIR = join(__dirname, '..', 'tests', 'data', 'socrates_gp');

const args = process.argv.slice(2);
const getArg = (name, def) => { const i = args.indexOf(name); return i >= 0 && args[i+1] ? args[i+1] : def; };
const TOP = parseInt(getArg('--top', '50000'));
const RATE_MS = parseInt(getArg('--rate', '80'));
const START = parseInt(getArg('--start', '0'));

mkdirSync(DATA_DIR, { recursive: true });

// Rotate Tor circuit via control port
async function rotateTorCircuit() {
  return new Promise((resolve, reject) => {
    const client = net.connect(9051, '127.0.0.1', () => {
      client.write('AUTHENTICATE ""\r\n');
    });
    let buf = '';
    client.on('data', (data) => {
      buf += data.toString();
      if (buf.includes('250 OK') && !buf.includes('NEWNYM')) {
        client.write('SIGNAL NEWNYM\r\n');
      } else if (buf.includes('250 OK') && buf.includes('NEWNYM')) {
        client.end();
        resolve();
      }
    });
    client.on('error', reject);
    setTimeout(() => { client.end(); resolve(); }, 3000);
  });
}

// Parse CSV
const csvPath = existsSync(join(__dirname, '..', 'tests', 'data', 'socrates_current.csv'))
  ? join(__dirname, '..', 'tests', 'data', 'socrates_current.csv')
  : join(__dirname, '..', 'tests', 'data', 'socrates_maxprob.csv');
console.log(`CSV: ${csvPath}`);
const csv = readFileSync(csvPath, 'utf8');
const lines = csv.trim().split('\n');

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
  if (!seen.has(key)) { seen.add(key); pairs.push({ id1, id2 }); }
}

console.log(`Total pairs: ${pairs.length}, starting at offset ${START}`);

let downloaded = 0, cached = 0, errors = 0, rotations = 0;
let consecutiveErrors = 0;

const agent = new SocksProxyAgent('socks5h://127.0.0.1:9050');

for (let i = START; i < pairs.length; i++) {
  const { id1, id2 } = pairs[i];
  const outFile = join(DATA_DIR, `gp_${id1},${id2}.json`);

  if (existsSync(outFile)) { cached++; continue; }

  const url = `https://celestrak.org/SOCRATES/data.php?CATNR=${id1},${id2}&FORMAT=json`;
  try {
    const controller = new AbortController();
    const timeout = setTimeout(() => controller.abort(), 15000);
    
    const res = await fetch(url, { 
      agent,
      signal: controller.signal,
      headers: { 'User-Agent': 'Mozilla/5.0 (compatible)' }
    });
    clearTimeout(timeout);

    if (res.status === 429 || res.status === 503) {
      // Rate limited — rotate circuit and retry
      rotations++;
      console.log(`  ⟳ Rate limited at ${i}, rotating circuit... (${rotations})`);
      await rotateTorCircuit();
      await new Promise(r => setTimeout(r, 5000)); // Wait for new circuit
      i--; // Retry this one
      continue;
    }

    if (!res.ok) {
      errors++;
      consecutiveErrors++;
      if (consecutiveErrors > 10) {
        console.log(`  ⟳ ${consecutiveErrors} consecutive errors, rotating...`);
        await rotateTorCircuit();
        await new Promise(r => setTimeout(r, 3000));
        consecutiveErrors = 0;
      }
      continue;
    }

    const text = await res.text();
    const data = JSON.parse(text);
    if (!data || data.length < 2) { errors++; continue; }
    writeFileSync(outFile, text);
    downloaded++;
    consecutiveErrors = 0;
  } catch (e) {
    errors++;
    consecutiveErrors++;
    if (consecutiveErrors > 5) {
      rotations++;
      console.log(`  ⟳ Error streak at ${i}: ${e.message}, rotating... (${rotations})`);
      await rotateTorCircuit();
      await new Promise(r => setTimeout(r, 3000));
      consecutiveErrors = 0;
    }
  }

  const total = downloaded + errors;
  if (total % 500 === 0 || i === pairs.length - 1) {
    console.log(`  [${i+1}/${pairs.length}] ${downloaded} new, ${cached} cached, ${errors} err, ${rotations} rotations`);
  }

  await new Promise(r => setTimeout(r, RATE_MS));
}

console.log(`\nDone: ${downloaded} downloaded, ${cached} cached, ${errors} errors, ${rotations} rotations`);
console.log(`Total GP files: ${downloaded + cached}`);
