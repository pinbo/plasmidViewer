'use strict';
/* ---------- cloning design tools: homology-based (In-Fusion / Gibson), Golden Gate, Gateway (BP / LR) ----------
   Every design works on opened plasmids and returns the final plasmid (sequence + carried-over features) plus the primers. */

/* ===== primer helpers ===== */
const primerTm = seq => { const t = tmPrimer3(seq); return Number.isFinite(t) ? t : 0; };   // primer3-method Tm (conditions: see Tm settings)
/* gene-specific 5' part of a primer: shortest 18–40 nt stretch reaching the target Tm, preferring a G/C 3' end */
function pickGS(seq, target) {
  seq = seq.toUpperCase();
  const max = Math.min(40, seq.length); let L = Math.min(18, max);
  while (L < max && primerTm(seq.slice(0, L)) < target) L++;
  for (let k = 0; k < 3 && L < max && !/[GC]$/.test(seq.slice(0, L)); k++) L++;
  return seq.slice(0, L);
}

/* ===== sequence / feature plumbing ===== */
const cloneFeat = f => { const c = JSON.parse(JSON.stringify(f)); delete c.id; return c; };
function circSlice(vec, a, b) {
  const U = vec.seq, n = U.length; if (b <= a) return '';
  if (!vec.circular) { if (a < 0 || b > n) throw new Error('Not enough vector sequence next to the junction for the requested homology arm.'); return U.slice(a, b); }
  let o = ''; for (let i = a; i < b; i++) o += U[((i % n) + n) % n]; return o;
}
const circLocs = (a, b, N) => (b > N ? [[a, N], [0, b - N]] : a < 0 ? [[N + a, N], [0, b]] : [[a, b]]);

/* features of `doc` carried into new coordinates; ranges: [{a, b, off}] (doc coords → off). Features not fully inside a retained range are dropped. */
function carryFeatures(doc, ranges) {
  const out = []; let dropped = 0;
  for (const f of doc.features) {
    const locs = []; let ok = true;
    for (const [x, y] of f.locs) { const r = ranges.find(r => x >= r.a && y <= r.b); if (!r) { ok = false; break; } locs.push([x - r.a + r.off, y - r.a + r.off]); }
    if (ok) { const c = cloneFeat(f); c.locs = locs; out.push(c); } else dropped++;
  }
  return { feats: out, dropped };
}
/* new plasmid = retained backbone (forward path keepFrom → keepTo) + `mid` inserted between keepTo and keepFrom */
function backboneAssembly(vec, keepFrom, keepTo, mid) {
  const U = vec.seq, n = U.length; let seq, ranges, midStart;
  if (keepTo <= keepFrom) {
    seq = U.slice(0, keepTo) + mid.text + U.slice(keepFrom);
    ranges = [{ a: 0, b: keepTo, off: 0 }, { a: keepFrom, b: n, off: keepTo + mid.text.length }]; midStart = keepTo;
  } else {
    if (!vec.circular) throw new Error('The chosen cut positions would remove the vector ends – use a circular vector, or pick the other fragment.');
    seq = U.slice(keepFrom, keepTo) + mid.text; ranges = [{ a: keepFrom, b: keepTo, off: 0 }]; midStart = keepTo - keepFrom;
  }
  const { feats, dropped } = carryFeatures(vec, ranges);
  for (const f of mid.features) { const c = cloneFeat(f); c.locs = c.locs.map(([a, b]) => [a + midStart, b + midStart]); feats.push(c); }
  return { seq, features: feats, midStart, dropped };
}
function primerFeature(name, seq, strand, locs, note, unmatched5) {
  const quals = { note: [`${note || name}: 5'-${seq}-3'`] };
  if (unmatched5) { quals.unmatched5 = [String(unmatched5)]; quals.note.push(`the first ${unmatched5} nt at the 5' end do not match this plasmid (drawn as a zigzag)`); }
  return { name, type: 'primer_bind', strand, locs, color: '#6b7ae0', quals };
}
function joinParts(parts) {
  let text = '', features = []; const offs = [];
  for (const p of parts) { offs.push(text.length); for (const f of p.features) { const c = cloneFeat(f); c.locs = c.locs.map(([a, b]) => [a + text.length, b + text.length]); features.push(c); } text += p.text; }
  return { text, features, offs };
}
/* a plasmid doc view in the orientation where att1 … att2 reads forward (used by Gateway) */
function rcDoc(doc) {
  const clip = revcompClip({ text: doc.seq, features: doc.features.map(cloneFeat) });
  return makeDoc({ name: doc.name, seq: clip.text, circular: doc.circular, features: clip.features, meta: {} });
}

/* resolve a part spec {src, name, rc, text} into {name, text, features} */
function resolvePart(spec, docs) {
  let text, features, name = spec.name || 'insert';
  if (spec.src === 'text') { text = cleanSeq(spec.text || ''); features = []; }
  else {
    const m = spec.src.split(':'), doc = docs.find(d => d.id === +m[1]); if (!doc) throw new Error('Source plasmid is no longer open.');
    let a, b;
    if (m[0] === 'sel') { const r = selRange(doc); if (!r) throw new Error(`Nothing is selected in “${doc.name}”.`); [a, b] = r; }
    else { const f = doc.features.find(x => x.id === +m[2]); if (!f) throw new Error('Feature not found.'); [a, b] = featBounds(f); }
    text = doc.seq.slice(a, b); features = clipFeatures(doc, a, b);
  }
  if (/[^ACGTacgt]/.test(text)) throw new Error(`Insert “${name}” contains non-ATGC characters – clean it up first.`);
  if (text.length < 30) throw new Error(`Insert “${name}” is too short (${text.length} bp).`);
  if (spec.rc) { const c = revcompClip({ text, features }); text = c.text; features = c.features; }
  return { name, text, features };
}

/* ===================== homology-based: In-Fusion / Gibson / HiFi ===================== */
/* opts: vec, leftEnd, rightStart (insert goes between), parts, H (homology length), tm,
   addL / addR (optional): vector bases re-added next to the junctions so the restriction sites are regenerated */
function designHomology(o) {
  const { vec, parts, H } = o, notes = [], warnings = [];
  const addL = (o.addL || '').toUpperCase(), addR = (o.addR || '').toUpperCase();
  const leftFlank = circSlice(vec, o.leftEnd - H, o.leftEnd).toUpperCase() + addL, rightFlank = addR + circSlice(vec, o.rightStart, o.rightStart + H).toUpperCase();
  const mid = joinParts([{ text: o.addL || '', features: [] }, ...parts, { text: o.addR || '', features: [] }]), primers = [];
  parts.forEach((p, i) => {
    const U = p.text.toUpperCase(), tailF = i === 0 ? leftFlank : parts[i - 1].text.toUpperCase().slice(-H);
    const gsF = pickGS(U, o.tm), gsR = pickGS(revcomp(U), o.tm), tailR = i === parts.length - 1 ? revcomp(rightFlank) : '';
    if (i > 0 && parts[i - 1].text.length < H) warnings.push(`“${parts[i - 1].name}” is shorter than the ${H} nt overlap.`);
    primers.push({ name: `${p.name}_F`, seq: tailF + gsF, tail: tailF.length, gs: gsF, strand: 1, part: i }, { name: `${p.name}_R`, seq: tailR + gsR, tail: tailR.length, gs: gsR, strand: -1, part: i });
  });
  const asm = backboneAssembly(vec, o.rightStart, o.leftEnd, mid), N = asm.seq.length, feats = asm.features;
  primers.forEach(pr => {
    const s = asm.midStart + mid.offs[pr.part + 1], e = s + parts[pr.part].text.length;
    const locs = pr.strand === 1 ? circLocs(s - pr.tail, s + pr.gs.length, N) : circLocs(e - pr.gs.length, e + pr.tail, N);
    feats.push(primerFeature(pr.name, pr.seq, pr.strand, locs));
  });
  if (o.siteNote) notes.push(o.siteNote);
  notes.push(`Each insert is amplified with ${H}-nt homology arms: the first forward primer matches the vector upstream of the junction, each further forward primer matches the end of the previous insert, and the last reverse primer matches the vector downstream.`);
  if (asm.dropped) warnings.push(`${asm.dropped} vector feature(s) lie in the removed region or span the junction and were not carried over.`);
  return { seq: asm.seq, circular: vec.circular, features: feats, primers, notes, warnings, method: 'homology cloning' };
}

/* ===================== Golden Gate ===================== */
/* a = nt between the end of the recognition core and the top-strand cut; b = same for the bottom strand */
function ggParams(e) { const end = e.len - e.trail; return { a: e.top - end, b: e.bot - end, core: e.site.slice(e.lead, end) }; }
function goldenGateEligible(e) { const { a, b } = ggParams(e); return !e.palin && a >= 0 && b > a && b - a <= 5; }
const GG_OVERHANGS = ['AATG', 'GCTT', 'GGAG', 'TACT', 'CAGC', 'CGCT', 'AGGT', 'GCAG', 'TGCC', 'ACTA', 'TTAC', 'CCAA', 'GTTC', 'ATCC', 'CTCG', 'TCGA'];

/* opts: vec, enzyme (name), parts, tm, pad, keep ('auto'|'A'|'B'), internal (array of user overhangs, may be empty) */
function designGoldenGate(o) {
  const { vec, parts } = o, notes = [], warnings = [];
  const e = ENZYMES.find(x => x.name === o.enzyme); if (!e || !goldenGateEligible(e)) throw new Error('Choose a Type IIS enzyme that cuts downstream of its site and leaves a 5′ overhang.');
  const { a, b, core } = ggParams(e), ovLen = b - a;
  const res = analyzeEnzymes(vec).get(e.name), cuts = res ? res.cuts : [];
  if (cuts.length !== 2) throw new Error(`${e.name} cuts this vector ${cuts.length} time${cuts.length === 1 ? '' : 's'}; a Golden Gate destination vector needs exactly 2 sites flanking the dropout region.`);
  const [j1, j2] = cuts.slice().sort((x, y) => x.top - y.top);
  const lo = j => Math.min(j.top, j.bot), hi = j => Math.max(j.top, j.bot), n = vec.seq.length;
  const inA = c => c.s >= lo(j1) && c.e <= hi(j2) + 0;   // fragment A = j1 → j2 (no origin crossing)
  const sitesInA = cuts.filter(inA).length, sitesInB = cuts.length - sitesInA;
  let keep = o.keep;
  if (keep === 'auto') {
    if (sitesInA === cuts.length) keep = 'B'; else if (sitesInB === cuts.length) keep = 'A';
    else throw new Error(`Cannot tell which fragment is the dropout (${e.name} sites are not both inside one fragment) – choose the fragment to keep.`);
  }
  // retained path runs rs → le; the insert is placed between le and rs
  const rs = keep === 'A' ? j1 : j2, le = keep === 'A' ? j2 : j1;
  const U = vec.seq.toUpperCase(), ovL = U.slice(lo(le), hi(le)), ovR = U.slice(lo(rs), hi(rs));
  if (ovL.length !== ovLen || ovR.length !== ovLen) throw new Error('Unexpected overhang length at a vector cut.');
  const keepFrom = hi(rs), keepTo = lo(le);
  const retainedLen = keepTo >= keepFrom ? keepTo - keepFrom : n - keepFrom + keepTo;
  if (retainedLen < 100) warnings.push('The retained backbone is very short – you may have chosen the wrong fragment.');
  // junction overhangs: ovL, internal…, ovR
  const k = parts.length, junc = [ovL];
  const used = new Set([ovL, ovR, revcomp(ovL), revcomp(ovR)]);
  for (let i = 1; i < k; i++) {
    let ov = (o.internal[i - 1] || '').toUpperCase().replace(/[^ACGT]/g, '');
    if (!ov) { ov = GG_OVERHANGS.map(x => x.slice(0, ovLen)).find(x => !used.has(x) && !used.has(revcomp(x)) && x !== revcomp(x)); if (!ov) throw new Error('No unused overhang left – supply internal overhangs manually.'); }
    if (ov.length !== ovLen) throw new Error(`Internal overhangs must be ${ovLen} nt for ${e.name}.`);
    if (used.has(ov) || used.has(revcomp(ov))) warnings.push(`Overhang ${ov} (junction ${i}) is not unique – assembly may misligate.`);
    if (ov === revcomp(ov)) warnings.push(`Overhang ${ov} is palindromic.`);
    used.add(ov); used.add(revcomp(ov)); junc.push(ov);
  }
  junc.push(ovR);
  if (ovL === ovR || ovL === revcomp(ovR)) warnings.push('The two vector overhangs are identical – the insert could ligate in either direction.');
  // primers
  const spacer = 'ACGT'.slice(0, a), pad = (o.pad || '').toUpperCase().replace(/[^ACGT]/g, ''), head = pad + core + spacer, primers = [];
  parts.forEach((p, i) => {
    const Up = p.text.toUpperCase(), gsF = pickGS(Up, o.tm), gsR = pickGS(revcomp(Up), o.tm);
    const hits = [...Up.matchAll(e.regex)].length + (e.palin ? 0 : [...Up.matchAll(e.rcRegex)].length);
    if (hits) warnings.push(`Insert “${p.name}” contains ${hits} ${e.name} site${hits > 1 ? 's' : ''} – it would be cut during assembly. Remove it (silent mutation) first.`);
    primers.push({ name: `${p.name}_F`, seq: head + junc[i] + gsF, tail: junc[i].length, gs: gsF, strand: 1, part: i, ov: junc[i] },
      { name: `${p.name}_R`, seq: head + revcomp(junc[i + 1]) + gsR, tail: junc[i + 1].length, gs: gsR, strand: -1, part: i, ov: junc[i + 1] });
  });
  // final sequence: backbone + ovL + p1 + ov12 + p2 … + pk + ovR
  const pieces = [{ text: junc[0], features: [] }];
  parts.forEach((p, i) => { pieces.push(p); pieces.push({ text: junc[i + 1], features: [] }); });
  const mid = joinParts(pieces), partOff = parts.map((p, i) => mid.offs[1 + 2 * i]);
  const asm = backboneAssembly(vec, keepFrom, keepTo, mid), N = asm.seq.length, feats = asm.features;
  primers.forEach(pr => {
    const s = asm.midStart + partOff[pr.part], e2 = s + parts[pr.part].text.length;
    const locs = pr.strand === 1 ? circLocs(s - pr.tail, s + pr.gs.length, N) : circLocs(e2 - pr.gs.length, e2 + pr.tail, N);
    feats.push(primerFeature(pr.name, pr.seq, pr.strand, locs, null, head.length));   // pad + site + spacer are cut off, only overhang + binding part remain
  });
  if (asm.dropped) warnings.push(`${asm.dropped} vector feature(s) lie in the dropout region or span a cut and were not carried over.`);
  notes.push(`${e.name} (${core}, ${a} nt spacer, ${ovLen}-nt 5′ overhang). Vector overhangs: ${ovL} → insert → ${ovR}. Primers = ${pad ? pad + ' + ' : ''}${core} + ${spacer} + overhang + gene-specific sequence; the pad, site and spacer are removed by the digestion and are not in the final plasmid. Junction overhangs: ${junc.join(' | ')}.`);
  return { seq: asm.seq, circular: vec.circular, features: feats, primers, notes, warnings, method: 'Golden Gate (' + e.name + ')' };
}

/* ===================== Gateway (BP / LR) ===================== */
const ATT1_CORE = 'TTTGTACAAAAAAGC', ATT2_CORE = 'GCTTTCTTGTACAAA';        // crossover cores, att1 and att2 (as read with the gene)
const ATTB1_PRIMER = 'GGGGACAAGTTTGTACAAAAAAGCAGGCT', ATTB2_PRIMER = 'GGGGACCACTTTGTACAAGAAAGCTGGGT';

/* locate the two crossover cores in a plasmid; returns a doc in the orientation att1 … att2 plus the core start positions */
function findAtt(doc, label) {
  const U = seqU(doc), find = s => { const out = []; let i = U.indexOf(s); while (i !== -1) { out.push(i); i = U.indexOf(s, i + 1); } return out; };
  const f1 = find(ATT1_CORE), r1 = find(revcomp(ATT1_CORE)), f2 = find(ATT2_CORE), r2 = find(revcomp(ATT2_CORE));
  const n1 = f1.length + r1.length, n2 = f2.length + r2.length;
  if (n1 !== 1 || n2 !== 1) {
    // explain near-misses: the 15-nt crossover core must match exactly
    const near = (core, label) => {
      const hits = [...scanDNA(U, doc.circular, core, 3), ...scanDNA(U, doc.circular, revcomp(core), 3)].filter(([, m]) => m > 0);
      return hits.length ? ` Near-matches for the ${label} core: ${hits.slice(0, 3).map(([i, m]) => `${i + 1} (${m} mismatch${m > 1 ? 'es' : ''})`).join(', ')} – the recombination core must be exactly ${core}.` : '';
    };
    throw new Error(`“${doc.name}” must contain exactly one att1 and one att2 crossover core (found ${n1} and ${n2}).` + (n1 !== 1 ? near(ATT1_CORE, 'att1') : '') + (n2 !== 1 ? near(ATT2_CORE, 'att2') : ''));
  }
  if (f1.length && f2.length) return { doc, p1: f1[0], p2: f2[0] };
  if (r1.length && r2.length) { const d = rcDoc(doc), N = d.seq.length; return { doc: d, p1: N - r1[0] - 15, p2: N - r2[0] - 15 }; }
  throw new Error(`The att1 and att2 sites in “${doc.name}” point in opposite directions – not a valid Gateway vector.`);
}

/* text + features of the forward path a → b of a doc; the path may run across the origin of a circular plasmid */
function pathSegment(doc, a, b) {
  const U = doc.seq, n = U.length;
  if (a <= b) return { text: U.slice(a, b), features: clipFeatures(doc, a, b) };
  if (!doc.circular) throw new Error(`The att sites of “${doc.name}” are in the wrong order for a linear sequence.`);
  return { text: U.slice(a) + U.slice(0, b), features: carryFeatures(doc, [{ a, b: n, off: 0 }, { a: 0, b, off: n - a }]).feats };
}

/* BP: opts {donor, part, extraF, extraR, tm}.  LR: opts {entry, dest} */
function designGateway(o) {
  const notes = [], warnings = [], primers = [];
  let asm, name;
  if (o.reaction === 'BP') {
    const part = o.part, U = part.text.toUpperCase(), tm = o.tm;
    const tailF = ATTB1_PRIMER + o.extraF, tailR = ATTB2_PRIMER + o.extraR, gsF = pickGS(U, tm), gsR = pickGS(revcomp(U), tm);
    primers.push({ name: `${part.name}_attB1_F`, seq: tailF + gsF, tail: tailF.length, gs: gsF, strand: 1 }, { name: `${part.name}_attB2_R`, seq: tailR + gsR, tail: tailR.length, gs: gsR, strand: -1 });
    const product = tailF + part.text + revcomp(tailR), P = product.toUpperCase();
    const c1 = P.indexOf(ATT1_CORE), c2 = P.indexOf(ATT2_CORE);
    if (c1 < 0 || c2 < 0 || c2 < c1) throw new Error('Internal error: att cores not found in the PCR product.');
    const donor = findAtt(o.donor), D = donor.doc, shift = tailF.length - (c1 + 15);
    const mid = { text: product.slice(c1 + 15, c2), features: part.features.map(f => { const c = cloneFeat(f); c.locs = c.locs.map(([a, b]) => [a + shift, b + shift]); return c; }) };
    asm = backboneAssembly(D, donor.p2, donor.p1 + 15, mid);
    notes.push('BP reaction: the attB1/attB2 PCR product recombines with the attP sites of the donor vector; the region between the attP crossover cores (ccdB cassette) is replaced by your insert, giving an entry clone with attL1/attL2.');
    if (o.extraF.length % 3 !== 0 || o.extraR.length % 3 !== 0) notes.push('Check the reading frame: the attB1 primer’s extra bases should keep the insert in frame for downstream LR reactions.');
    name = `${D.name.replace(/_rc$/, '')}_${part.name}_entry`;
    const N = asm.seq.length, s = asm.midStart + (tailF.length - (c1 + 15)), e = s + part.text.length;
    // the attB primers only match the final plasmid from the crossover core inwards; the rest of each tail is drawn as an unmatched zigzag
    const hangF = c1, hangR = P.length - (c2 + 15), mR = tailR.length - hangR;
    asm.features.push(primerFeature(primers[0].name, primers[0].seq, 1, circLocs(s - (tailF.length - hangF), s + gsF.length, N), null, hangF), primerFeature(primers[1].name, primers[1].seq, -1, circLocs(e - gsR.length, e + mR, N), null, hangR));
  } else {
    const entry = findAtt(o.entry), dest = findAtt(o.dest);
    const E = entry.doc, D = dest.doc;
    const seg = pathSegment(E, entry.p1 + 15, entry.p2);    // attL1 core end → attL2 core start (may cross the origin)
    asm = backboneAssembly(D, dest.p2, dest.p1 + 15, seg);
    notes.push('LR reaction: the insert between the attL1/attL2 crossover cores of the entry clone replaces the region between the attR1/attR2 cores of the destination vector (e.g. the ccdB cassette), producing an expression clone with attB1/attB2.');
    name = `${D.name.replace(/_rc$/, '')}_${E.name.replace(/_rc$/, '').replace(/_entry$/, '')}`;
  }
  // annotate attB sites in the product
  const t = makeDoc({ name: 'tmp', seq: asm.seq, circular: true, features: JSON.parse(JSON.stringify(asm.features)) });
  const attLib = [{ name: 'attB1', type: 'misc_recomb', color: '#e08fb5', seq: 'ACAAGTTTGTACAAAAAAGCAGGCT' }, { name: 'attB2', type: 'misc_recomb', color: '#e08fb5', seq: 'ACCCAGCTTTCTTGTACAAAGTGGT' }];
  for (const h of detectCandidates(t, attLib, 1)) if (!h.present) asm.features.push(candidateToFeature(h));
  if (asm.dropped) warnings.push(`${asm.dropped} backbone feature(s) lay in the replaced region or spanned the recombination site and were not carried over (e.g. ccdB, attP/attR annotations).`);
  return { seq: asm.seq, circular: true, features: asm.features, primers, notes, warnings, method: 'Gateway ' + o.reaction, name };
}

/* ===================== dialog ===================== */
function openCloningDialog() {
  if (!App.docs.length) return toast('Open a plasmid first');
  const docs = App.docs.slice(), cur = App.cur;
  let method = 'homology', design = null, modal = null;
  const sel = (opts, value) => { const s = el('select', {}, opts.map(([v, t]) => el('option', { value: v }, t))); if (value !== undefined) s.value = value; return s; };
  const docOpts = docs.map(d => [d.id, `${d.name} (${d.seq.length.toLocaleString('en-US')} bp)`]);
  const num = (v, min, max, w = 70) => el('input', { type: 'number', min, max, value: v, style: `width:${w}px` });
  const txt = (v, ph, w) => el('input', { value: v, placeholder: ph || '', style: w ? `width:${w}px` : '' });
  const getDoc = s => docs.find(d => d.id === +s.value);

  /* ----- insert parts ----- */
  const sourceOptions = () => {
    const o = [];
    for (const d of docs) { const r = selRange(d); if (r) o.push([`sel:${d.id}`, `Selection in “${d.name}” (${r[0] + 1}..${r[1]}, ${(r[1] - r[0]).toLocaleString('en-US')} bp)`]); }
    for (const d of docs) for (const f of d.features) o.push([`feat:${d.id}:${f.id}`, `Feature: ${f.name} — ${d.name}`]);
    o.push(['text', 'Pasted sequence…']); return o;
  };
  const partsHost = el('div', { class: 'parts' }); const rows = [];
  const addPart = () => {
    const src = sel(sourceOptions()), name = txt('insert', 'name', 110), rc = el('input', { type: 'checkbox' }), ta = el('textarea', { rows: 2, class: 'mono', placeholder: 'Paste the insert DNA (5′→3′)…', style: 'display:none;width:100%' });
    const row = { src, name, rc, ta };
    const sync = () => {
      ta.style.display = src.value === 'text' ? 'block' : 'none';
      if (src.value.startsWith('feat:')) { const m = src.value.split(':'), d = docs.find(x => x.id === +m[1]), f = d && d.features.find(x => x.id === +m[2]); if (f) { name.value = f.name.replace(/[^\w.-]+/g, '_'); rc.checked = f.strand === -1; } }
      else if (src.value.startsWith('sel:')) name.value = 'insert' + (rows.indexOf(row) + 1 || 1);
    };
    src.addEventListener('change', sync);
    const first = !rows.length; if (first) { const d0 = docs.find(d => selRange(d)); if (d0) src.value = `sel:${d0.id}`; else if (src.options.length > 1) src.selectedIndex = 0; }
    row.el = el('div', { class: 'partrow' }, el('div', { class: 'prow' }, el('span', { class: 'pn' }, ''), src, name, el('label', { class: 'chk', title: 'Use the reverse complement of this region' }, rc, 'rev-comp'),
      el('button', { class: 'tx', title: 'Remove this part', onclick: () => { if (rows.length > 1) { rows.splice(rows.indexOf(row), 1); row.el.remove(); renumber(); } } }, '✕')), ta);
    rows.push(row); partsHost.append(row.el); sync(); renumber();
  };
  const renumber = () => rows.forEach((r, i) => { r.el.querySelector('.pn').textContent = rows.length > 1 ? i + 1 : ''; });
  const specs = () => rows.map(r => ({ src: r.src.value, name: r.name.value.trim() || 'insert', rc: r.rc.checked, text: r.ta.value }));
  const addBtn = el('button', { class: 'btn sm', onclick: addPart }, '+ Add insert');

  /* ----- method panels ----- */
  const panel = el('div', { class: 'cpanel' });
  const ctl = {};
  /* choices = unique cutters + the user's selected enzymes; an enzyme with several sites gets one choice per site (value "name|index") */
  const enzChoices = d => {
    const res = analyzeEnzymes(d), mine = new Set(App.settings.enzShow), out = [];
    for (const [name, r] of res) {
      const k = r.cuts.length; if (!k || (k > 1 && !mine.has(name))) continue;
      const cuts = r.cuts.map((c, i) => ({ c, i })).sort((x, y) => x.c.top - y.c.top);
      for (const { c, i } of cuts) out.push({ value: `${name}|${i}`, text: k === 1 ? name : `${name} – site at ${c.top} (${k} sites)`, name, k, sort: name + '\0' + String(c.top).padStart(9, '0') });
    }
    return out.sort((x, y) => x.sort.localeCompare(y.sort));
  };
  const fillEnz = () => {
    const d = getDoc(ctl.vec), list = enzChoices(d);
    for (const s of [ctl.enz1, ctl.enz2]) { const v = s.value; s.innerHTML = ''; s.append(el('option', { value: '' }, s === ctl.enz2 ? '— none (single cut) —' : '— choose —')); for (const o of list) s.append(el('option', { value: o.value }, o.text)); if (list.some(o => o.value === v)) s.value = v; }
    const uq = list.filter(o => o.k === 1).length;
    ctl.enzInfo.textContent = `${uq} unique cutter${uq === 1 ? '' : 's'}${list.length > uq ? ` + ${list.length - uq} site(s) of your selected enzymes` : ''} in ${d.name}`;
  };
  const fillGG = () => {
    const d = getDoc(ctl.vec), res = analyzeEnzymes(d), v = ctl.ggEnz.value; ctl.ggEnz.innerHTML = '';
    const list = ENZYMES.filter(goldenGateEligible).sort((a, b) => (GOLDEN_GATE.has(b.name) - GOLDEN_GATE.has(a.name)) || a.name.localeCompare(b.name));
    for (const e of list) { const k = (res.get(e.name) || { cuts: [] }).cuts.length; ctl.ggEnz.append(el('option', { value: e.name }, `${e.name}  (${k} site${k === 1 ? '' : 's'} in vector)`)); }
    const best = list.find(e => (res.get(e.name) || { cuts: [] }).cuts.length === 2 && GOLDEN_GATE.has(e.name)) || list.find(e => (res.get(e.name) || { cuts: [] }).cuts.length === 2);
    ctl.ggEnz.value = v && list.some(e => e.name === v) && v !== '' ? v : (best ? best.name : 'BsaI');
  };
  const buildPanel = () => {
    panel.innerHTML = ''; const f = (l, c, h) => panel.append(field(l, c, h));
    if (method === 'homology') {
      ctl.vec = sel(docOpts, cur.id); f('Vector (backbone)', ctl.vec);
      ctl.mode = sel([['enz', 'Linearize with restriction enzyme(s)'], ['sel', 'Insert at cursor / replace the selection']], 'enz'); f('Vector opening', ctl.mode);
      ctl.enz1 = sel([['', '']]); ctl.enz2 = sel([['', '']]); ctl.swap = el('input', { type: 'checkbox' }); ctl.keepSites = el('input', { type: 'checkbox' }); ctl.enzInfo = el('small', {});
      ctl.enzBox = el('div', {}, el('div', { class: 'two' }, field('Enzyme 1', ctl.enz1), field('Enzyme 2 (optional)', ctl.enz2)), el('label', { class: 'chk' }, ctl.swap, 'Remove the other fragment (the one that spans the origin)'),
        el('label', { class: 'chk' }, ctl.keepSites, 'Keep the restriction sites in the final plasmid (regenerate them at both junctions)'), ctl.enzInfo);
      panel.append(ctl.enzBox);
      ctl.selInfo = el('div', { class: 'small-note' }); panel.append(ctl.selInfo);
      ctl.H = num(15, 8, 60); ctl.tm = num(60, 45, 72);
      panel.append(el('div', { class: 'two' }, field('Homology arm length (nt)', ctl.H, 'In-Fusion: 15 · Gibson / HiFi: 20–40'), field('Primer Tm target (°C)', ctl.tm)));
      const upd = () => { const sm = ctl.mode.value === 'sel'; ctl.enzBox.style.display = sm ? 'none' : 'block'; ctl.selInfo.style.display = sm ? 'block' : 'none'; const d = getDoc(ctl.vec), r = selRange(d); ctl.selInfo.textContent = r ? `Replaces ${r[0] + 1}..${r[1]} (${r[1] - r[0]} bp) of ${d.name}.` : `Inserts at the cursor (after base ${d.caret}) of ${d.name}.`; };
      ctl.mode.addEventListener('change', upd); ctl.vec.addEventListener('change', () => { fillEnz(); upd(); }); fillEnz(); upd();
    } else if (method === 'gg') {
      ctl.vec = sel(docOpts, cur.id); f('Destination vector', ctl.vec);
      ctl.ggEnz = sel([['BsaI', 'BsaI']]); f('Type IIS enzyme', ctl.ggEnz, 'Eligible enzymes only (cut downstream of the site, 5′ overhang). The vector needs exactly 2 sites.');
      ctl.keep = sel([['auto', 'Automatic (keep the fragment without enzyme sites)'], ['A', 'Keep fragment A (between the sites, no origin crossing)'], ['B', 'Keep fragment B (spans the origin)']]); f('Backbone fragment', ctl.keep);
      ctl.pad = txt('TTTT', '', 120); ctl.tm = num(60, 45, 72); ctl.internal = txt('', 'auto', 260);
      panel.append(el('div', { class: 'two' }, field('5′ pad (helps the enzyme cut)', ctl.pad), field('Primer Tm target (°C)', ctl.tm)), field('Internal overhangs (comma-separated, blank = automatic)', ctl.internal, 'Only needed with 2+ inserts: one overhang per junction between inserts.'));
      ctl.vec.addEventListener('change', fillGG); fillGG();
    } else {
      ctl.rxn = sel([['BP', 'BP reaction – PCR product (attB) + donor vector (attP)'], ['LR', 'LR reaction – entry clone (attL) + destination vector (attR)']], 'BP'); f('Reaction', ctl.rxn);
      const bp = el('div', {}), lr = el('div', {});
      ctl.vec = sel(docOpts, cur.id); bp.append(field('Donor vector (attP1 … attP2)', ctl.vec));
      ctl.extraF = txt('', 'e.g. GCCACC (Kozak)', 160); ctl.extraR = txt('', "e.g. TTA (stop, written 5′→3′ in the primer)", 220); ctl.tm = num(60, 45, 72);
      bp.append(el('div', { class: 'two' }, field('Extra bases after attB1 (forward primer)', ctl.extraF), field('Extra bases after attB2 (reverse primer)', ctl.extraR)), field('Primer Tm target (°C)', ctl.tm));
      ctl.entry = sel(docOpts, cur.id); ctl.dest = sel(docOpts, (docs.find(d => d !== cur) || cur).id);
      lr.append(el('div', { class: 'two' }, field('Entry clone (attL1 … attL2)', ctl.entry), field('Destination vector (attR1 … attR2)', ctl.dest)), el('div', { class: 'small-note' }, 'No primers are needed for an LR reaction.'));
      panel.append(bp, lr);
      const upd = () => { const isBP = ctl.rxn.value === 'BP'; bp.style.display = isBP ? 'block' : 'none'; lr.style.display = isBP ? 'none' : 'block'; partsWrap.style.display = isBP ? 'block' : 'none'; };
      ctl.rxn.addEventListener('change', upd); setTimeout(upd, 0);
    }
    addBtn.style.display = method === 'gateway' ? 'none' : '';
  };

  const partsWrap = el('div', { class: 'partswrap' }, el('h4', {}, 'Insert(s)'), partsHost, addBtn);
  const tabs = el('div', { class: 'ctabs' });
  const setMethod = m => { method = m; $$('button', tabs).forEach(b => b.classList.toggle('on', b.dataset.m === m)); partsWrap.style.display = 'block'; if (m === 'gateway') while (rows.length > 1) { rows.pop().el.remove(); } buildPanel(); out.innerHTML = ''; design = null; };
  for (const [m, t] of [['homology', 'In-Fusion / Gibson'], ['gg', 'Golden Gate'], ['gateway', 'Gateway']]) tabs.append(el('button', { 'data-m': m, class: 'ctab', onclick: () => setMethod(m) }, t));

  /* ----- design + results ----- */
  const out = el('div', { class: 'cresult' });
  const run = () => {
    out.innerHTML = ''; design = null;
    try {
      const T = () => clamp(+ctl.tm.value || 60, 45, 72);
      if (method === 'gateway' && ctl.rxn.value === 'LR') design = designGateway({ reaction: 'LR', entry: getDoc(ctl.entry), dest: getDoc(ctl.dest) });
      else {
        const parts = specs().map(s => resolvePart(s, docs)), vec = getDoc(ctl.vec);
        if (method === 'homology') {
          const H = clamp(+ctl.H.value || 15, 8, 60); let leftEnd, rightStart, addL = '', addR = '', siteNote = '';
          if (ctl.mode.value === 'sel') { const r = selRange(vec); [leftEnd, rightStart] = r || [vec.caret, vec.caret]; }
          else {
            if (!ctl.enz1.value) throw new Error('Choose at least one restriction enzyme for the vector (or use the selection mode).');
            const cutObj = v => { const [n, i] = v.split('|'); return { name: n, ...analyzeEnzymes(vec).get(n).cuts[+i] }; };
            const ca = cutObj(ctl.enz1.value), cb = ctl.enz2.value ? cutObj(ctl.enz2.value) : ca;
            if (ctl.enz2.value && ca.top === cb.top) throw new Error('The two enzymes cut at the same position.');
            const [cx, cy] = ca.top <= cb.top ? [ca, cb] : [cb, ca];
            const [cl, cr] = ctl.swap.checked ? [cy, cx] : [cx, cy];     // enzyme at the left / right junction of the insert
            leftEnd = cl.top; rightStart = cr.top;
            if (ctl.keepSites.checked) {
              const ok = c => c.s < c.top && c.top < c.e;                  // the cut must lie inside the recognition site
              if (!ok(cl) || !ok(cr)) throw new Error(`Restriction sites can only be regenerated for enzymes that cut inside their recognition site (not ${!ok(cl) ? cl.name : cr.name}).`);
              addL = vec.seq.slice(cl.top, cl.e); addR = vec.seq.slice(cr.s, cr.top);
              siteNote = `Restriction sites kept: ${cl.name} (${vec.seq.slice(cl.s, cl.e).toUpperCase()}) is regenerated at the left junction${cl === cr ? ' and again' : ', '} ${cr.name} (${vec.seq.slice(cr.s, cr.e).toUpperCase()}) at the right junction. The primers carry the missing half of each site after the vector homology, so the final plasmid can be cut with the same enzymes.`;
            }
          }
          design = designHomology({ vec, leftEnd, rightStart, parts, H, tm: T(), addL, addR, siteNote });
        } else if (method === 'gg') {
          design = designGoldenGate({ vec, parts, enzyme: ctl.ggEnz.value, tm: T(), pad: ctl.pad.value, keep: ctl.keep.value, internal: ctl.internal.value.split(/[,\s]+/).filter(Boolean) });
        } else {
          design = designGateway({ reaction: 'BP', donor: vec, part: parts[0], extraF: cleanSeq(ctl.extraF.value).toUpperCase(), extraR: cleanSeq(ctl.extraR.value).toUpperCase(), tm: T() });
        }
      }
      design.name = design.name || (getDoc(ctl.vec) ? getDoc(ctl.vec).name : 'cloned') + '_' + (rows.length ? rows.map(r => r.name.value.trim() || 'insert').join('_') : 'insert');
    } catch (e) { out.append(el('div', { class: 'warnline', style: 'display:block' }, '✖ ' + e.message)); return; }
    renderResult();
  };
  const renderResult = () => {
    const d = design; out.innerHTML = '';
    out.append(el('div', { class: 'okline' }, `✔ ${d.method}: new ${d.circular ? 'circular' : 'linear'} plasmid “${d.name}”, ${d.seq.length.toLocaleString('en-US')} bp, ${d.features.length} features.`));
    for (const w of d.warnings) out.append(el('div', { class: 'warnline', style: 'display:block' }, '⚠ ' + w));
    for (const n of d.notes) out.append(el('div', { class: 'small-note' }, n));
    if (d.primers.length) {
      const tbl = el('table', { class: 'ptable' }, el('tr', {}, ['Primer', "Sequence (5′→3′)", 'Length', 'Tm of binding part'].map(h => el('th', {}, h))));
      for (const p of d.primers) {
        const tail = p.seq.slice(0, p.seq.length - p.gs.length);
        tbl.append(el('tr', {}, el('td', {}, el('b', {}, p.name)), el('td', { class: 'mono pseq' }, tail ? el('span', { class: 'tail' }, tail) : '', p.gs), el('td', {}, p.seq.length + ' nt'), el('td', {}, fmtTm(tmPrimer3(p.gs)))));
      }
      out.append(el('div', { class: 'ptablewrap' }, tbl), el('div', { class: 'small-note' }, 'Coloured part = 5′ extension (homology / enzyme site / att site); black = gene-specific binding part. Tm: ' + tmConfText() + '.'));
    }
    const csv = () => 'Name,Sequence,Length,Tm\n' + d.primers.map(p => `${p.name},${p.seq},${p.seq.length},${fmtTm(tmPrimer3(p.gs)).replace(' °C', '')}`).join('\n');
    out.append(el('div', { class: 'btnrow' },
      el('button', { class: 'btn primary', onclick: () => { const doc = addDoc({ name: d.name, seq: d.seq, circular: d.circular, features: d.features, meta: { definition: `Designed with Plasmid Viewer – ${d.method}` } }); doc.dirty = true; toast(`Created “${doc.name}”`); if (modal) modal.close(); } }, 'Create plasmid in a new tab'),
      d.primers.length ? el('button', { class: 'btn', onclick: async () => { await copyToSystem(d.primers.map(p => `${p.name}\t${p.seq}`).join('\n')); toast('Primers copied'); } }, 'Copy primers') : null,
      d.primers.length ? el('button', { class: 'btn', onclick: () => downloadBlob(d.name + '_primers.csv', csv(), 'text/csv') }, 'Download primers (CSV)') : null));
  };

  const tmLine = el('div', { class: 'small-note' }, 'Primer Tm: ' + tmConfText() + '  ', el('a', { href: '#', onclick: e => { e.preventDefault(); openTmDialog(); } }, 'change conditions'), ' (applies to the next design).');
  const body = el('div', { class: 'cloning' }, tabs, panel, tmLine, partsWrap, el('div', { class: 'btnrow' }, el('button', { class: 'btn primary', onclick: run }, 'Design')), out);
  addPart(); setMethod('homology');
  modal = openModal('Cloning tools', body, [{ label: 'Close' }], { wide: true });
}
actions.cloning = openCloningDialog;
