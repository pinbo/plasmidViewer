'use strict';
/* ---------- NEB enzyme information: parsing, extra enzymes, pop-up content ---------- */

const NEB = {};
for (const line of NEB_RAW.split('\n')) {
  const f = line.split('~'); if (f.length < 10) continue;
  const [name, site, buf, acts, heat, temp, dil, meth, sub, notes, flags] = f, a = acts.split('|');
  NEB[name] = {
    name, site, buf, act: { 'r1.1': a[0], 'r2.1': a[1], 'r3.1': a[2], 'rCutSmart': a[3] }, heat, temp, dil,
    meth: { Dam: meth[0], Dcm: meth[1], CpG: meth[2] }, sub, notes: notes ? notes.split(',').map(s => s.trim()) : [], flags: flags || '',
  };
}
const nebBase = n => n.replace(/-(HFv2|HF|v2)$/, '');

/* returns {main, variants} for an enzyme name (HF / v2 versions are listed as variants) */
function nebFor(name) {
  const variants = Object.keys(NEB).filter(k => k !== name && nebBase(k) === name).map(k => NEB[k]);
  let main = NEB[name];
  if (!main && variants.length) main = variants.shift();
  return main ? { main, variants } : null;
}

/* NEB notation ("G/GATCC", "GGTCTC(1/5)", "CCGC(-3/-1)") → our "G^GATCC" / "GGTCTCN^NNNN_" notation */
function nebToDef(seq) {
  let site, top, bot = null, m;
  if ((m = seq.match(/^([A-Z]+)\((-?\d+)\/(-?\d+)\)$/))) { site = m[1]; top = site.length + +m[2]; bot = site.length + +m[3]; }
  else if (/^[A-Z]*\/[A-Z]*$/.test(seq)) { top = seq.indexOf('/'); site = seq.replace('/', ''); }
  else return null;
  if (top < 0 || (bot !== null && bot < 0) || site.length > 30) return null;
  const len = Math.max(site.length, top, bot || 0), chars = site.padEnd(len, 'N');
  let out = '';
  for (let i = 0; i <= len; i++) { if (i === top) out += '^'; if (i === bot) out += '_'; if (i < len) out += chars[i]; }
  return out;
}
const NEB_SKIP = new Set(['MspJI', 'DpnI']);
/* NEB enzymes that are not yet in the built-in table (base names only) */
function nebExtraDefs() {
  const out = [];
  for (const k of Object.keys(NEB)) {
    const base = nebBase(k); if (NEB_SKIP.has(base) || /^(Nb|Nt|I-|PI-)/.test(base)) continue;
    const def = nebToDef(NEB[k].site) || nebToDef(NEB[base] ? NEB[base].site : '') || (NEB[k].site.match(/^[A-Z]+$/) && NEB[base] ? nebToDef(NEB[base].site) : null);
    if (def) out.push([base, def]);
  }
  return out;
}
buildEnzymes();

/* ---------- pop-up content ---------- */
const METH = { N: 'not sensitive', I: 'impaired', B: 'blocked', O: 'impaired if overlapping' };
const NOTE_TXT = {
  a: 'Ligation is <10 %', b: 'Ligation 25–75 %', c: 'Recutting after ligation <5 %', d: 'Recutting after ligation 50–75 %', e: 'Ligation/recutting not applicable (nicking, methylation-affected or cuts outside its site)',
  1: 'Star activity possible with extended digestion, high enzyme or >5 % glycerol', 2: 'Star activity possible with extended digestion', 3: 'Star activity possible with >5 % glycerol',
};
const bufName = b => (/^r\d/.test(b) ? 'NEBuffer ' + b : /^rCutSmart/.test(b) ? b.replace('rCutSmart', 'rCutSmart Buffer') : /Buffer/.test(b) ? b : 'NEBuffer ' + b);
function overhangInfo(e) {
  const d = e.bot - e.top;
  return d === 0 ? 'Blunt' : d > 0 ? `5′ overhang (${d} nt)` : `3′ overhang (${-d} nt)`;
}
function isTypeIIS(e) { return e.top < e.lead || e.top > e.len - e.trail || e.bot < e.lead || e.bot > e.len - e.trail; }

function enzymeTipHTML(doc, name) {
  const e = ENZYMES.find(x => x.name === name); if (!e) return '';
  const res = doc ? analyzeEnzymes(doc).get(name) : null, cuts = res ? res.cuts : [], k = cuts.length;
  const row = (a, b, cls = '') => `<div class="tr ${cls}"><span>${a}</span><b>${b}</b></div>`;
  let h = `<div class="th"><i style="background:var(--enz)"></i>${esc(name)}${GOLDEN_GATE.has(name) ? '<em class="tag">Golden Gate</em>' : ''}</div>`;
  h += row('Sites in this plasmid', k ? `${k} (${cuts.slice(0, 10).map(c => c.top).join(', ')}${k > 10 ? ', …' : ''})` : 'none', 'hl');
  h += row('Recognition', `<code>${esc(e.def)}</code>`);
  h += row('Cuts', `${overhangInfo(e)} · ${isTypeIIS(e) ? 'Type IIS (cuts outside site)' : 'Type II'}`);
  const n = nebFor(name);
  if (n) {
    const m = n.main;
    h += `<div class="sec">NEB data${m.name !== name ? ` (${esc(m.name)})` : ''}</div>`;
    h += row('Supplied buffer', esc(bufName(m.buf)));
    h += `<div class="acts">${['r1.1', 'r2.1', 'r3.1', 'rCutSmart'].map(b => { const v = m.act[b] || ''; const num = parseFloat(v.replace('<', '')) || 0, cls = /^<10$/.test(v) || num < 10 ? 'lo' : num >= 75 ? 'hi' : 'mid'; return `<span class="${cls}"><i>${b.replace('rCutSmart', 'CutSmart')}</i>${esc(v)}</span>`; }).join('')}</div>`;
    h += row('Incubate / inactivate', `${esc(m.temp)} °C · ${m.heat === 'No' ? 'not heat-inactivated' : esc(m.heat) + ' °C'}`);
    const sens = Object.entries(m.meth).filter(([, v]) => v !== 'N').map(([a, v]) => `${a} ${METH[v] || v}`);
    h += row('Methylation', sens.length ? esc(sens.join(' · ')) : 'not sensitive to Dam / Dcm / CpG');
    if (m.notes.length) h += `<div class="tn">${m.notes.map(x => esc(NOTE_TXT[x] || x)).join('<br>')}</div>`;
    if (/\*/.test(Object.values(m.act).join(''))) h += `<div class="tn">* may show star activity in that buffer</div>`;
    if (n.variants.length) h += row('High-fidelity / updated', esc(n.variants.map(v => `${v.name} (${bufName(v.buf)})`).join(', ')));
    h += row('Unit assay substrate', esc(m.sub));
  } else h += `<div class="tn">${e.user ? 'Custom enzyme' : 'No NEB data for this enzyme'}</div>`;
  return h;
}
