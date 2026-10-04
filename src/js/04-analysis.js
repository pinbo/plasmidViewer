'use strict';
/* ---------- restriction-site search, feature detection, ORFs, find ---------- */

/* upper-case copy of the sequence for all analyses (doc.seq itself keeps the user's case) */
function seqU(doc) {
  const c = doc._cache;
  if (c.Uf !== doc.seq) { c.U = doc.seq.toUpperCase(); c.Uf = doc.seq; }
  return c.U;
}
/* positions of characters that are not A/C/G/T (either case) */
function badBases(doc) {
  const c = doc._cache;
  if (!c.bad) { c.bad = []; for (const m of doc.seq.matchAll(/[^ACGTacgt]/g)) c.bad.push(m.index); }
  return c.bad;
}

/* returns Map name -> {enzyme, cuts:[{top, bot, s, e, strand}]}  (positions are 0-based boundaries) */
function analyzeEnzymes(doc) {
  const c = doc._cache;
  if (c.enz && c.enzV === enzVersion) return c.enz;
  const seq = seqU(doc), n = seq.length, res = new Map();
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
  for (let i = 0; i < lim; i++) {
    let mis = 0, j = 0;
    for (; j < L; j++) { if (s.charCodeAt(i + j) !== pat.charCodeAt(j) && ++mis > maxMis) break; }
    if (j === L) hits.push([i, mis]);
  }
  return hits;
}

function segsFor(start, len, n) {
  return start + len > n ? [[start, n], [0, start + len - n]] : [[start, start + len]];
}

const locsOverlap = (a, b) => a.some(([x, y]) => b.some(([p, q]) => x < q && p < y));

/* All library matches with >= thr identity (mismatch-only alignment, both strands; protein motifs in all 6 frames).
   Each hit: {entry, locs, strand, mis, len, identity (0-1), present (an annotation with the same name already overlaps)} */
function detectCandidates(doc, lib, thr = 0.96) {
  const seq = seqU(doc), n = seq.length, found = [];
  for (const entry of lib) {
    if (entry.seq) {
      const pat = entry.seq.toUpperCase(), L = pat.length;
      const maxMis = Math.floor(L * (1 - thr) + 1e-9);
      const rc = revcomp(pat), palin = rc === pat;
      for (const [i, mis] of scanDNA(seq, doc.circular, pat, maxMis)) found.push({ entry, locs: segsFor(i, L, n), strand: 1, mis, len: L });
      if (!palin) for (const [i, mis] of scanDNA(seq, doc.circular, rc, maxMis)) found.push({ entry, locs: segsFor(i, L, n), strand: -1, mis, len: L });
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
            if (s >= 0 && s < n) found.push({ entry, locs: segsFor(s, L3, n), strand, mis: 0, len: L3 });
            idx = prot.indexOf(aa, idx + 1);
          }
        }
      }
    }
  }
  // keep the best hit when the same library entry matches the same place more than once
  found.sort((a, b) => a.mis / a.len - b.mis / b.len);
  const keep = [];
  for (const h of found) {
    const nm = h.entry.name.toLowerCase();
    if (keep.some(k => k.entry.name.toLowerCase() === nm && locsOverlap(k.locs, h.locs))) continue;
    h.identity = 1 - h.mis / h.len;
    h.present = doc.features.some(f => f.name.toLowerCase() === nm && locsOverlap(f.locs, h.locs));
    keep.push(h);
  }
  // entries sharing a `group` (e.g. the near-identical Gateway att sites) are mutually exclusive where they overlap: the best match wins
  const groups = keep.filter(h => h.entry.group).sort((a, b) => b.identity - a.identity || b.len - a.len), won = [];
  for (const h of groups) {
    const clash = won.some(w => w.entry.group === h.entry.group && h.locs.some(([x, y]) => w.locs.some(([p, q]) => { const ov = Math.min(y, q) - Math.max(x, p); return ov > 0.3 * Math.min(y - x, q - p); })));
    if (!clash) won.push(h);
  }
  const lose = new Set(groups.filter(h => !won.includes(h)));
  return keep.filter(h => !lose.has(h)).sort((a, b) => a.locs[0][0] - b.locs[0][0]);
}

function candidateToFeature(h) {
  const note = `auto-detected from library (${(h.identity * 100).toFixed(1)}% identity)`;
  return { name: h.entry.name, type: h.entry.type || 'misc_feature', strand: h.entry.dir === false ? 0 : h.strand, locs: h.locs, color: h.entry.color || typeColor(h.entry.type), quals: { note: [note] } };
}

/* new features only (used for automatic annotation when opening un-annotated files) */
function detectFeatures(doc, lib, thr = App.settings.detectThr) {
  return detectCandidates(doc, lib, thr).filter(h => !h.present).map(candidateToFeature);
}

/* ---------- ORFs ---------- */
function findORFs(doc, minAA) {
  const c = doc._cache; const key = 'orf' + minAA;
  if (c[key]) return c[key];
  const seq = seqU(doc), n = seq.length, out = [];
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
  const q = cleanSeq(query).toUpperCase(); if (!q) return [];
  const U = seqU(doc), n = U.length, ext = doc.circular ? U + U.slice(0, q.length - 1) : U;
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
  const U = seqU(doc);
  for (const p of pos) bases += f.strand === -1 ? (COMP[U[p]] || 'N') : U[p];
  const out = [];
  for (let i = 0; i + 2 < pos.length; i += 3) out.push({ aa: codonToAA(bases.substr(i, 3)), pos: [pos[i], pos[i + 1], pos[i + 2]] });
  doc._cache[key] = out; return out;
}


/* CDS without a terminal stop codon: continue translating in frame (downstream, same strand) until the next stop codon,
   the start of another CDS that is in frame with it, or the end of the sequence.
   Returns a pseudo-feature {ext:true, parent, locs, aaOffset, stopped, endCds} or null. */
function cdsExtension(doc, f) {
  const key = 'ext' + f.id + ':' + f.locs.map(l => l.join('-')).join(',') + f.strand;
  if (key in doc._cache) return doc._cache[key];
  let res = null;
  const cod = cdsCodons(doc, f), total = f.locs.reduce((s, l) => s + l[1] - l[0], 0);
  if (cod.length && total % 3 === 0 && cod[cod.length - 1].aa !== '*') {
    const n = doc.seq.length, U = seqU(doc), dir = f.strand === -1 ? -1 : 1, limit = Math.min(n - total, 30000);
    let p = cod[cod.length - 1].pos[2], stopped = false, endCds = '';
    // first base (in reading direction) of every other CDS on the same strand: reaching one at a codon boundary means it is in frame
    const starts = new Map();
    for (const g of doc.features) if (g !== f && g.type === 'CDS' && g.strand === f.strand) starts.set(dir === 1 ? g.locs[0][0] : g.locs[g.locs.length - 1][1] - 1, g.name);
    const pos = [];
    while (pos.length < limit) {
      const trio = [];
      for (let k = 0; k < 3; k++) {
        p += dir;
        if (doc.circular) p = ((p % n) + n) % n;
        else if (p < 0 || p >= n) break;
        trio.push(p);
      }
      if (trio.length < 3) break;
      if (starts.has(trio[0])) { endCds = starts.get(trio[0]); break; }
      const bases = trio.map(q => (dir === -1 ? COMP[U[q]] || 'N' : U[q])).join('');
      pos.push(...trio);
      if (codonToAA(bases) === '*') { stopped = true; break; }
    }
    if (pos.length) {
      const runs = []; let cur = [pos[0], pos[0]];
      for (let i = 1; i < pos.length; i++) { if (pos[i] === pos[i - 1] + dir) cur[1] = pos[i]; else { runs.push(cur); cur = [pos[i], pos[i]]; } }
      runs.push(cur);
      const locs = runs.map(([a, b]) => (dir === 1 ? [a, b + 1] : [b, a + 1]));
      if (dir === -1) locs.reverse();
      res = { id: 'ext' + f.id, name: 'in frame with ' + f.name, type: 'CDS', ext: true, parent: f, strand: f.strand, color: f.color, locs, aaOffset: cod.length, stopped, endCds, quals: {} };
    }
  }
  doc._cache[key] = res; return res;
}
