// Compares tmPrimer3() with reference values produced by primer3-py (tools/check-tm.js <reference.json>)
const fs = require('fs'), vm = require('vm'), path = require('path');
const ctx = vm.createContext({ console, localStorage: undefined });
vm.runInContext(fs.readFileSync(path.join(__dirname, '..', 'src', 'js', '01-core.js'), 'utf8').replace(/^'use strict';/, ''), ctx);
const ref = JSON.parse(fs.readFileSync(process.argv[2], 'utf8'));
let worst = 0, bad = 0;
for (const r of ref) {
  const mine = vm.runInContext(`tmPrimer3(${JSON.stringify(r.seq)}, ${JSON.stringify(r.conf)})`, ctx), d = Math.abs(mine - r.tm);
  worst = Math.max(worst, d); if (d > 0.05) { bad++; console.log('DIFF', r.seq, r.seq.length, JSON.stringify(r.conf), mine.toFixed(2), r.tm.toFixed(2)); }
}
console.log(`${ref.length} comparisons, max |Δ| = ${worst.toFixed(3)} °C, ${bad} over 0.05 °C`);
