// Cross-checks the built-in enzyme table against NEB's recognition/cut notation.  usage: node tools/check-enzymes.js
const fs = require('fs'), vm = require('vm'), path = require('path');
const ctx = vm.createContext({ console, localStorage: undefined });
for (const f of ['01-core', '03-library', '03b-neb-data', '03c-neb']) {
  const code = fs.readFileSync(path.join(__dirname, '..', 'src', 'js', f + '.js'), 'utf8').replace(/^'use strict';/, '');
  vm.runInContext(code, ctx, { filename: f });
}
const out = vm.runInContext(`(() => {
  const bad = [], ok = [], extra = [];
  const mine = new Map(ENZYMES.map(e => [e.name, e]));
  const own = new Map(); for (const line of ENZYME_DEFS.split('\\n')) { const [n, d] = line.trim().split(/\\s+/); if (n && d) own.set(n, d); }
  for (const [n, d] of nebExtraDefs()) extra.push(n);
  for (const k of Object.keys(NEB)) {
    const base = nebBase(k); if (!own.has(base)) continue;
    const nd = nebToDef(NEB[k].site); if (!nd) continue;
    const a = compileEnzyme(base, own.get(base)), b = compileEnzyme(base, nd);
    const same = a.site === b.site && a.top === b.top && a.bot === b.bot;
    (same ? ok : bad).push(same ? base : base + ': mine ' + own.get(base) + ' | NEB ' + nd);
  }
  return JSON.stringify({ total: ENZYMES.length, agree: ok.length, bad, extraCount: extra.length, extra });
})()`, ctx);
console.log(JSON.stringify(JSON.parse(out), null, 1));
