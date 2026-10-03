'use strict';
/* ---------- restriction-site search, feature detection, ORFs, find ---------- */

/* returns Map name -> {enzyme, cuts:[{top, bot, s, e, strand}]}  (positions are 0-based boundaries) */
function analyzeEnzymes(doc) {
  const c = doc._cache;
  if (c.enz && c.enzV === enzVersion) return c.enz;
  const seq = doc.seq, n = seq.length, res = new Map();
  for (const e of ENZYMES) {
    const ext = doc.circular ? seq + seq.slice(0, e.len - 1) : seq;
    const cuts = [];
    const add = (i, strand) => {
      let top, bot;
      let cs = i + (strand === 1 ? e.lead : e.rcLead), ce = i + e.len - (strand === 1 ? e.trail : e.rcTrail);   // recognition core (without flanking N spacers)
      if (doc.circular && cs >= n) { cs -= n; ce -= n; }
      if (strand === 1) { top = i + e.top; bot = i + e.bot; } else { top = i + e.len - e.bot; bot = i + e.len - e.top; }
      if (doc.circular) { top = ((top % n) + n) % n; bot = ((bot % n) + n) % n; }
      else if (top < 0 || top > n || bot < 0 || bot > n) return;
      cuts.push({ top, bot, s: cs, e: ce, strand });
    };
    for (const m of ext.matchAll(e.regex)) if (m.index < n) add(m.index, 1);
    if (!e.palin) for (const m of ext.matchAll(e.rcRegex)) if (m.index < n) add(m.index, -1);
    cuts.sort((a, b) => a.top - b.top);
    res.set(e.name, { enzyme: e, cuts });
  }
  c.enz = res; c.enzV = enzVersion;
  return res;
}

/* enzymes to show on maps, based on display mode */
function enzymeView(doc) {
  const mode = App.settings.enzMode, mine = new Set(App.settings.enzShow), out = [];
  if (mode === 'none') return out;
  for (const [name, r] of analyzeEnzymes(doc)) {
    const k = r.cuts.length; if (!k) continue;
    const show = mode === 'all' || (mode === 'unique' && k === 1) || (mode === 'dual' && k <= 2) || (mode === 'gg' && GOLDEN_GATE.has(name)) || (mode === 'custom' && mine.has(name));
    if (!show) continue;
    for (const cut of r.cuts) out.push({ name, count: k, ...cut });
  }
  return out;
}

/* ---------- feature detection ---------- */
function scanDNA(seq, circular, pat, maxMis) {
  const L = pat.length, n = seq.length, hits = [];
  if (L > n) return hits;
  const s = circular ? seq + seq.slice(0, L - 1) : seq, lim = circular ? n : n - L + 1;
  const p0 = pat.charCodeAt(0);
  for (let i = 0; i < lim; i++) {
    let mis = 0, j = 0;
    for (; j < L; j++) { if (s.charCodeAt(i + j) !== pat.charCodeAt(j) && ++mis > maxMis) break; }
    if (j === L) hits.push(i);
  }
  return hits;
}

function segsFor(start, len, n) {
  return start + len > n ? [[start, n], [0, start + len - n]] : [[start, start + len]];
}

function detectFeatures(doc, lib) {
  const seq = doc.seq, n = seq.length, found = [];
  for (const entry of lib) {
    if (entry.seq) {
      const pat = entry.seq.toUpperCase(), L = pat.length;
      const maxMis = L < 30 ? 0 : Math.floor(L * 0.04);
      const rc = revcomp(pat), palin = rc === pat;
      for (const i of scanDNA(seq, doc.circular, pat, maxMis)) found.push({ entry, locs: segsFor(i, L, n), strand: 1 });
      if (!palin) for (const i of scanDNA(seq, doc.circular, rc, maxMis)) found.push({ entry, locs: segsFor(i, L, n), strand: -1 });
    } else if (entry.aa) {
      const aa = entry.aa.toUpperCase(), L3 = aa.length * 3;
      const ext = doc.circular ? seq + seq.slice(0, L3 - 1) : seq;
      for (const strand of [1, -1]) {
        const t = strand === 1 ? ext : revcomp(ext), m = t.length;
        for (let f = 0; f < 3; f++) {
          const prot = translate(t.slice(f));
          let idx = prot.indexOf(aa);
          while (idx !== -1) {
            const r = f + idx * 3;
            const s = strand === 1 ? r : m - (r + L3);
            if (s >= 0 && s < n) found.push({ entry, locs: segsFor(s, L3, n), strand });
            idx = prot.indexOf(aa, idx + 1);
          }
        }
      }
    }
  }
  // de-duplicate against existing annotations and among themselves
  const overlap = (a, b) => a.some(([x, y]) => b.some(([p, q]) => x < q && p < y));
  const keep = [];
  for (const h of found) {
    const nm = h.entry.name.toLowerCase();
    const dup = doc.features.concat(keep).some(f => f.name.toLowerCase() === nm && overlap(f.locs, h.locs));
    if (!dup) keep.push({ name: h.entry.name, type: h.entry.type || 'misc_feature', strand: h.entry.dir === false ? 0 : h.strand, locs: h.locs, color: h.entry.color || typeColor(h.entry.type), quals: { note: ['auto-detected'] } });
  }
  return keep;
}

/* ---------- ORFs ---------- */
function findORFs(doc, minAA) {
  const c = doc._cache; const key = 'orf' + minAA;
  if (c[key]) return c[key];
  const seq = doc.seq, n = seq.length, out = [];
  for (const strand of [1, -1]) {
    const t = strand === 1 ? seq : revcomp(seq);
    for (let f = 0; f < 3; f++) {
      let start = -1;
      for (let i = f; i + 2 < n; i += 3) {
        const aa = codonToAA(t.substr(i, 3));
        if (start < 0 && aa === 'M') start = i;
        else if (aa === '*') {
          if (start >= 0 && (i - start) / 3 >= minAA) {
            const a = start, b = i + 3;
            out.push({ id: 'orf' + out.length, name: `ORF ${Math.floor((b - a) / 3) - 1} aa`, type: 'ORF', color: typeColor('ORF'), strand, orf: true, quals: {}, locs: [strand === 1 ? [a, b] : [n - b, n - a]] });
          }
          start = -1;
        }
      }
    }
  }
  c[key] = out; return out;
}

/* ---------- find ---------- */
function findAll(doc, query) {
  const q = cleanSeq(query); if (!q) return [];
  const n = doc.seq.length, ext = doc.circular ? doc.seq + doc.seq.slice(0, q.length - 1) : doc.seq;
  const res = [], seen = new Set();
  const src = iupacRegexSrc(q), rcSrc = iupacRegexSrc(revcomp(q));
  const run = (re, strand) => { for (const m of ext.matchAll(new RegExp('(?=(' + re + '))', 'g'))) if (m.index < n && !seen.has(m.index + ':' + strand)) { seen.add(m.index + ':' + strand); res.push({ s: m.index, e: Math.min(m.index + q.length, n), strand }); } };
  run(src, 1); if (revcomp(q) !== q) run(rcSrc, -1);
  return res.sort((a, b) => a.s - b.s);
}

/* codons of a CDS feature: [{aa, pos:[p0,p1,p2]}] in biological order */
function cdsCodons(doc, f) {
  const key = 'cd' + f.id + ':' + f.locs.map(l => l.join('-')).join(',') + f.strand;
  if (doc._cache[key]) return doc._cache[key];
  const pos = [];
  const ordered = f.strand === -1 ? f.locs.slice().reverse() : f.locs;
  for (const [a, b] of ordered) {
    if (f.strand === -1) for (let i = b - 1; i >= a; i--) pos.push(i); else for (let i = a; i < b; i++) pos.push(i);
  }
  let bases = '';
  for (const p of pos) bases += f.strand === -1 ? (COMP[doc.seq[p]] || 'N') : doc.seq[p];
  const out = [];
  for (let i = 0; i + 2 < pos.length; i += 3) out.push({ aa: codonToAA(bases.substr(i, 3)), pos: [pos[i], pos[i + 1], pos[i + 2]] });
  doc._cache[key] = out; return out;
}
