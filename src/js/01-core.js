'use strict';
/* ---------- tiny helpers ---------- */
const $ = (s, r = document) => r.querySelector(s);
const $$ = (s, r = document) => Array.from(r.querySelectorAll(s));
const esc = s => String(s ?? '').replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));

function el(tag, attrs, ...kids) {
  const e = document.createElement(tag);
  let value;
  for (const [k, v] of Object.entries(attrs || {})) {
    if (k === 'class') e.className = v;
    else if (k === 'style') e.style.cssText = v;
    else if (k === 'value') value = v;
    else if (k.startsWith('on')) e.addEventListener(k.slice(2), v);
    else if (v === true) e.setAttribute(k, '');
    else if (v !== false && v != null) e.setAttribute(k, v);
  }
  for (const kid of kids.flat()) if (kid != null && kid !== false) e.append(kid.nodeType ? kid : document.createTextNode(kid));
  if (value !== undefined) e.value = value;
  return e;
}

const store = {
  mem: {},   // in-memory fallback when localStorage is unavailable (private mode, sandboxed pages)
  get(k, d) {
    if (k in this.mem) return JSON.parse(this.mem[k]);
    try { const v = localStorage.getItem('pv.' + k); return v == null ? d : JSON.parse(v); } catch (e) { return d; }
  },
  set(k, v) {
    const j = JSON.stringify(v); this.mem[k] = j;
    try { localStorage.setItem('pv.' + k, j); } catch (e) { /* ignore */ }
  },
};

const clamp = (v, a, b) => Math.max(a, Math.min(b, v));
const debounce = (fn, ms) => { let t; return (...a) => { clearTimeout(t); t = setTimeout(() => fn(...a), ms); }; };

/* ---------- sequence helpers ---------- */
const COMP = { A: 'T', T: 'A', G: 'C', C: 'G', U: 'A', R: 'Y', Y: 'R', K: 'M', M: 'K', S: 'S', W: 'W', B: 'V', V: 'B', D: 'H', H: 'D', N: 'N' };
/* complement of one base; keeps the case of the input (a → t, A → T); unknown letters become N/n */
function compChar(c) {
  const u = COMP[c]; if (u) return u;
  const l = COMP[c.toUpperCase()];
  return l ? l.toLowerCase() : (c === c.toLowerCase() ? 'n' : 'N');
}
function revcomp(s) { let o = ''; for (let i = s.length - 1; i >= 0; i--) o += compChar(s[i]); return o; }
function complement(s) { let o = ''; for (let i = 0; i < s.length; i++) o += compChar(s[i]); return o; }
const IUPAC = { A: 'A', C: 'C', G: 'G', T: 'T', U: 'T', R: 'AG', Y: 'CT', S: 'CG', W: 'AT', K: 'GT', M: 'AC', B: 'CGT', D: 'AGT', H: 'ACT', V: 'ACG', N: 'ACGTN' };
const iupacRegexSrc = site => site.split('').map(c => '[' + (IUPAC[c] || c) + ']').join('');
/* Keeps letters only (FASTA headers, digits, spaces dropped) and PRESERVES CASE; RNA U/u becomes T/t.
   Letters other than A/C/G/T are kept so they can be flagged and highlighted – see nonATGCStats(). */
function cleanSeq(text) {
  return String(text).replace(/^>.*$/gm, '').replace(/[^A-Za-z]/g, '').replace(/U/g, 'T').replace(/u/g, 't');
}
function nonATGCStats(s) {
  const m = String(s).match(/[^ACGTacgt]/g) || [], c = {};
  for (const ch of m) c[ch.toUpperCase()] = (c[ch.toUpperCase()] || 0) + 1;
  const detail = Object.entries(c).sort((a, b) => b[1] - a[1]).slice(0, 6).map(([k, v]) => `${k}×${v}`).join(', ');
  return { count: m.length, detail };
}

const CODON_AA = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG';
const CIDX = { T: 0, C: 1, A: 2, G: 3 };
function codonToAA(c) {
  const a = CIDX[c[0].toUpperCase()], b = CIDX[c[1].toUpperCase()], d = CIDX[c[2].toUpperCase()];
  return (a === undefined || b === undefined || d === undefined) ? 'X' : CODON_AA[a * 16 + b * 4 + d];
}
function translate(s) { let o = ''; for (let i = 0; i + 2 < s.length; i += 3) o += codonToAA(s.substr(i, 3)); return o; }

function gcPercent(s) { if (!s.length) return 0; let g = 0; for (const c of s) if (c === 'G' || c === 'C' || c === 'g' || c === 'c') g++; return 100 * g / s.length; }
/* ---------- melting temperature: primer3 method ----------
   SantaLucia 1998 nearest-neighbour thermodynamics + SantaLucia salt correction, with primer3's conversion of Mg2+/dNTP into an
   equivalent monovalent concentration; oligos longer than 60 nt use primer3's long-sequence formula (as in oligotm.c / primer3-py). */
const NN_PARAMS = { AA: [-7.9, -22.2], TT: [-7.9, -22.2], AT: [-7.2, -20.4], TA: [-7.2, -21.3], CA: [-8.5, -22.7], TG: [-8.5, -22.7], GT: [-8.4, -22.4], AC: [-8.4, -22.4],
  CT: [-7.8, -21.0], AG: [-7.8, -21.0], GA: [-8.2, -22.2], TC: [-8.2, -22.2], CG: [-10.6, -27.2], GC: [-9.8, -24.4], GG: [-8.0, -19.9], CC: [-8.0, -19.9] };
const TM_DEFAULTS = { mv: 50, dv: 1.5, dntp: 0.6, dna: 50 };    // mM Na+/K+, mM Mg2+, mM dNTP, nM oligo (primer3 defaults)
const tmConf = () => Object.assign({}, TM_DEFAULTS, store.get('tmConf', {}));
function tmPrimer3(seq, c = tmConf()) {
  const s = String(seq).toUpperCase(), n = s.length;
  if (n < 2 || /[^ACGT]/.test(s)) return NaN;
  let dv = c.dv, dntp = c.dntp; if (dv === 0) dntp = 0; if (dv < dntp) dv = dntp;
  const mv = c.mv + (dv > 0 ? 120 * Math.sqrt(dv - dntp) : 0);                     // equivalent monovalent salt, mM
  if (n > 60) { let gc = 0; for (const ch of s) if (ch === 'G' || ch === 'C') gc++; return 81.5 + 16.6 * Math.log10(mv / 1000) + 41 * gc / n - 600 / n; }
  let H = 0, S = 0;
  for (let i = 0; i < n - 1; i++) { const v = NN_PARAMS[s.substr(i, 2)]; H += v[0]; S += v[1]; }
  for (const ch of [s[0], s[n - 1]]) { if (ch === 'G' || ch === 'C') { H += 0.1; S -= 2.8; } else { H += 2.3; S += 4.1; } }
  const sym = s === revcomp(s); if (sym) S -= 1.4;
  S += 0.368 * (n - 1) * Math.log(mv / 1000);
  return H * 1000 / (S + 1.9872 * Math.log(c.dna * 1e-9 / (sym ? 1 : 4))) - 273.15;
}
const meltingTemp = s => tmPrimer3(s);
const fmtTm = t => (Number.isFinite(t) ? t.toFixed(1) + ' °C' : '—');
const tmConfText = (c = tmConf()) => `primer3 method · ${c.mv} mM Na⁺, ${c.dv} mM Mg²⁺, ${c.dntp} mM dNTP, ${c.dna} nM oligo`;

/* ---------- feature colours / types ---------- */
const TYPE_COLORS = {
  promoter: '#2eb872', CDS: '#f29b34', terminator: '#e0525f', rep_origin: '#8f96a3', primer_bind: '#6b7ae0',
  protein_bind: '#a66bd6', misc_feature: '#7f9fb5', enhancer: '#35b0c9', polyA_signal: '#d9534f', regulatory: '#2eb872',
  gene: '#f2c14e', mRNA: '#c9a227', sig_peptide: '#c97bd1', misc_recomb: '#e08fb5', LTR: '#6f8f3d', "5'UTR": '#8cc4a3',
  "3'UTR": '#8cc4a3', mat_peptide: '#e9a15a', misc_binding: '#a66bd6', ORF: '#b9a8d9',
};
const FEATURE_TYPES = ['CDS', 'promoter', 'terminator', 'rep_origin', 'enhancer', 'polyA_signal', 'primer_bind', 'protein_bind', 'regulatory',
  'misc_feature', 'gene', 'mRNA', 'sig_peptide', 'mat_peptide', 'misc_recomb', 'misc_binding', 'LTR', "5'UTR", "3'UTR"];
const typeColor = t => TYPE_COLORS[t] || '#7f9fb5';

function textOn(hex) {
  const m = /^#?([0-9a-f]{6})$/i.exec(hex || ''); if (!m) return '#111';
  const n = parseInt(m[1], 16), r = n >> 16, g = (n >> 8) & 255, b = n & 255;
  return (0.299 * r + 0.587 * g + 0.114 * b) > 150 ? '#14181f' : '#ffffff';
}

/* ---------- misc ---------- */
let _mctx;
function textWidth(text, font = '12px system-ui, sans-serif') {
  _mctx = _mctx || document.createElement('canvas').getContext('2d');
  _mctx.font = font; return _mctx.measureText(text).width;
}
function niceStep(n, target = 10) {
  const raw = n / target, p = Math.pow(10, Math.floor(Math.log10(raw)));
  for (const m of [1, 2, 2.5, 5, 10]) if (m * p >= raw) return Math.max(1, m * p);
  return 10 * p;
}
function fmtBp(n) { return n.toLocaleString('en-US') + ' bp'; }

function toast(msg, ms = 2600, kind = '') {
  const t = el('div', { class: 'toast ' + kind }, msg);
  $('#toasts').append(t);
  setTimeout(() => t.classList.add('out'), ms);
  setTimeout(() => t.remove(), ms + 400);
}

function downloadBlob(name, data, type = 'text/plain') {
  const blob = new Blob([data], { type });
  const a = el('a', { href: URL.createObjectURL(blob), download: name });
  document.body.append(a); a.click();
  setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 1000);
}

/* Greedy lane assignment. items: [{id, ivs:[[a,b],...]}] (sorted internally). Returns {lane:Map, count} */
function assignLanes(items, gap = 0) {
  const sorted = items.slice().sort((x, y) => Math.min(...x.ivs.map(v => v[0])) - Math.min(...y.ivs.map(v => v[0])));
  const lanes = [], out = new Map();
  for (const it of sorted) {
    let k = 0;
    for (; k < lanes.length; k++) {
      if (!it.ivs.some(([a, b]) => lanes[k].some(([c, d]) => a < d + gap && c < b + gap))) break;
    }
    if (k === lanes.length) lanes.push([]);
    lanes[k].push(...it.ivs); out.set(it.id, k);
  }
  return { lane: out, count: lanes.length };
}

const AA_INFO = {
  A: ['Alanine', 'Ala'], R: ['Arginine', 'Arg'], N: ['Asparagine', 'Asn'], D: ['Aspartic acid', 'Asp'], C: ['Cysteine', 'Cys'], Q: ['Glutamine', 'Gln'],
  E: ['Glutamic acid', 'Glu'], G: ['Glycine', 'Gly'], H: ['Histidine', 'His'], I: ['Isoleucine', 'Ile'], L: ['Leucine', 'Leu'], K: ['Lysine', 'Lys'],
  M: ['Methionine', 'Met'], F: ['Phenylalanine', 'Phe'], P: ['Proline', 'Pro'], S: ['Serine', 'Ser'], T: ['Threonine', 'Thr'], W: ['Tryptophan', 'Trp'],
  Y: ['Tyrosine', 'Tyr'], V: ['Valine', 'Val'], '*': ['Stop codon', 'Stop'], X: ['Unknown (ambiguous codon)', 'Xaa'],
};
