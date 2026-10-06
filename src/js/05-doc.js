'use strict';
/* ---------- document model, editing operations, undo, clipboard ---------- */

const App = {
  docs: [], cur: null,
  settings: Object.assign({
    theme: 'auto', view: 'split', sideOpen: true, mapW: 0, sideW: 0, sideTab: 'features', enzMode: 'unique', enzShow: [],
    showFeatures: true, showOrfs: false, orfMin: 100, showTrans: true, autoDetect: true, detectThr: 0.96, linZoom: 1,
  }, store.get('settings', {})),
  clip: null,
  editing: false,   // sequence editing is locked until the user switches it on (guards against accidental typing)
};
let _lockToast = 0;
function lockedToast() {
  if (Date.now() - _lockToast > 1500) toast('Sequence editing is disabled – click the 🔒 button in the toolbar to allow it', 3200);
  _lockToast = Date.now();
  return false;
}
const saveSettings = () => store.set('settings', App.settings);

let _docUid = 0, _featUid = 0;
function normFeature(f) {
  f.id = ++_featUid;
  f.type = f.type || 'misc_feature';
  f.name = f.name || f.type;
  f.strand = f.strand === -1 ? -1 : f.strand === 0 ? 0 : 1;
  f.color = f.color || typeColor(f.type);
  f.quals = f.quals || {};
  f.locs = f.locs.map(([a, b]) => [a, b]);
  return f;
}
function makeDoc(p) {
  const doc = Object.assign({ name: 'Untitled', seq: '', circular: true, features: [], meta: {}, undo: [], redo: [], dirty: false, rev: 0, anchor: 0, caret: 0, selFid: null, handle: null, _cache: {} }, p);
  doc.id = ++_docUid;
  const n = doc.seq.length;
  doc.features = doc.features.filter(f => f.locs.every(([a, b]) => a >= 0 && b <= n && b > a)).map(normFeature);
  return doc;
}

const featBounds = f => [Math.min(...f.locs.map(l => l[0])), Math.max(...f.locs.map(l => l[1]))];
function selRange(doc) { const a = Math.min(doc.anchor, doc.caret), b = Math.max(doc.anchor, doc.caret); return a === b ? null : [a, b]; }

function setSel(doc, anchor, caret, fid = null, reveal = true) {
  const n = doc.seq.length;
  doc.anchor = clamp(anchor, 0, n); doc.caret = clamp(caret, 0, n); doc.selFid = fid;
  doc._reveal = reveal;
  requestUpdate(false);
}

/* ----- history ----- */
function snapshot(doc) {
  return { seq: doc.seq, features: JSON.stringify(doc.features), circular: doc.circular, anchor: doc.anchor, caret: doc.caret };
}
function restore(doc, s) {
  doc.seq = s.seq; doc.features = JSON.parse(s.features); doc.circular = s.circular; doc.anchor = s.anchor; doc.caret = s.caret; doc.selFid = null;
  doc.rev++; doc._cache = {}; doc.dirty = true; doc._reveal = true;
}
function undo(doc) { if (!doc.undo.length) return; doc.redo.push(snapshot(doc)); restore(doc, doc.undo.pop()); requestUpdate(true); }
function redo(doc) { if (!doc.redo.length) return; doc.undo.push(snapshot(doc)); restore(doc, doc.redo.pop()); requestUpdate(true); }

function mutate(doc, fn, { coalesce = false } = {}) {
  const now = Date.now();
  if (!(coalesce && doc._lastCoalesce && now - doc._lastCoalesce < 1200 && doc.undo.length)) {
    doc.undo.push(snapshot(doc)); if (doc.undo.length > 80) doc.undo.shift();
  }
  doc._lastCoalesce = coalesce ? now : 0;
  doc.redo = [];
  fn();
  doc.rev++; doc._cache = {}; doc.dirty = true; doc._reveal = true;
  requestUpdate(true);
}

/* ----- raw edit primitives (no history) ----- */
function rawInsert(doc, pos, text, feats) {
  const n = text.length;
  doc.seq = doc.seq.slice(0, pos) + text + doc.seq.slice(pos);
  for (const f of doc.features) {
    f.locs = f.locs.map(([a, b]) => {
      if (a >= pos) return [a + n, b + n];
      if (b > pos) return [a, b + n];   // insertion inside the feature → grows
      return [a, b];
    });
  }
  if (feats) for (const pf of feats) {
    const nf = normFeature({ ...JSON.parse(JSON.stringify(pf)), locs: pf.locs.map(([a, b]) => [a + pos, b + pos]) });
    doc.features.push(nf);
  }
}
function rawDelete(doc, s, e) {
  const n = e - s;
  doc.seq = doc.seq.slice(0, s) + doc.seq.slice(e);
  const out = [];
  for (const f of doc.features) {
    const locs = [];
    for (const [a, b] of f.locs) {
      const na = a < s ? a : a >= e ? a - n : s;
      const nb = b <= s ? b : b >= e ? b - n : s;
      if (nb > na) locs.push([na, nb]);
    }
    if (locs.length) { f.locs = locs; out.push(f); }
  }
  doc.features = out;
}

/* ----- public edit operations ----- */
function editReplace(doc, s, e, text, feats, opts) {
  if (!App.editing) return lockedToast();   // changing / inserting / deleting sequence is locked; feature edits are not
  text = cleanSeq(text);
  if (s === e && !text) return;
  mutate(doc, () => {
    if (e > s) rawDelete(doc, s, e);
    if (text) rawInsert(doc, s, text, feats);
    doc.anchor = doc.caret = s + text.length; doc.selFid = null;
  }, opts);
}
function deleteSelection(doc, forward = true) {
  const r = selRange(doc);
  if (r) return editReplace(doc, r[0], r[1], '');
  const n = doc.seq.length;
  if (forward && doc.caret < n) editReplace(doc, doc.caret, doc.caret + 1, '', null, { coalesce: true });
  else if (!forward && doc.caret > 0) editReplace(doc, doc.caret - 1, doc.caret, '', null, { coalesce: true });
}
function typeBases(doc, text) {
  const r = selRange(doc) || [doc.caret, doc.caret];
  editReplace(doc, r[0], r[1], text, null, { coalesce: true });
}

/* features completely inside [s,e), made relative */
function clipFeatures(doc, s, e) {
  const out = [];
  for (const f of doc.features) {
    if (f.locs.every(([a, b]) => a >= s && b <= e)) {
      const c = JSON.parse(JSON.stringify(f)); c.locs = f.locs.map(([a, b]) => [a - s, b - s]); delete c.id; out.push(c);
    }
  }
  return out;
}
function revcompClip(clip) {
  const L = clip.text.length;
  return {
    text: revcomp(clip.text),
    features: clip.features.map(f => ({ ...f, strand: -f.strand, locs: f.locs.map(([a, b]) => [L - b, L - a]).reverse() })),
  };
}

/* ----- clipboard (internal clip keeps features; plain text goes to system clipboard) ----- */
function currentClip(doc) {
  const r = selRange(doc); if (!r) return null;
  return { text: doc.seq.slice(r[0], r[1]), features: clipFeatures(doc, r[0], r[1]) };
}
function storeClip(clip) { App.clip = clip; store.set('clip', clip); }
function lookupClip(text) {
  const c = App.clip || store.get('clip', null);
  const t = cleanSeq(text), T = t.toUpperCase();   // case-insensitive match, but the pasted text keeps its own case
  if (c && c.text.toUpperCase() === T) return { ...c, text: t };
  if (c && c.text && revcomp(c.text).toUpperCase() === T) return { ...revcompClip(c), text: t };
  return { text: t, features: [] };
}
async function copyToSystem(text) {
  try { await navigator.clipboard.writeText(text); return true; } catch (e) { /* fall through */ }
  const ta = el('textarea', { style: 'position:fixed;opacity:0' }); ta.value = text; document.body.append(ta); ta.select();
  let ok = false; try { ok = document.execCommand('copy'); } catch (e) { /* ignore */ }
  ta.remove(); return ok;
}
async function doCopy(doc, cut = false) {
  const clip = currentClip(doc); if (!clip) return toast('Select some sequence first');
  storeClip(clip); await copyToSystem(clip.text);
  toast(`${cut ? 'Cut' : 'Copied'} ${fmtBp(clip.text.length)}${clip.features.length ? ` with ${clip.features.length} feature${clip.features.length > 1 ? 's' : ''}` : ''}`);
  if (cut) deleteSelection(doc);
}
async function doCopyRC(doc) {
  const clip = currentClip(doc); if (!clip) return toast('Select some sequence first');
  const rc = revcompClip(clip); storeClip(rc); await copyToSystem(rc.text); toast('Copied reverse complement');
}
function pasteClip(doc, clip) {
  if (!clip || !clip.text) return toast('Nothing to paste');
  const r = selRange(doc) || [doc.caret, doc.caret];
  if (!App.editing) return lockedToast();
  warnNonATGC(clip.text, 'in the pasted sequence');
  editReplace(doc, r[0], r[1], clip.text, clip.features);
  toast(`Pasted ${fmtBp(clip.text.length)}${clip.features.length ? ` with ${clip.features.length} feature${clip.features.length > 1 ? 's' : ''}` : ''}`);
}
async function doPaste(doc, rc = false) {
  let text = '';
  try { text = await navigator.clipboard.readText(); } catch (e) { /* denied → fall back to the internal clip */ }
  let clip = text ? lookupClip(text) : (App.clip || store.get('clip', null));
  if (!clip) return toast('Clipboard is empty (or access was blocked) – try Ctrl/Cmd+V');
  if (rc) clip = revcompClip(clip);
  pasteClip(doc, clip);
}

/* ----- other transformations ----- */
function reverseComplementRange(doc, s, e) {
  if (!App.editing) return lockedToast();
  mutate(doc, () => {
    const L = e - s, mid = revcomp(doc.seq.slice(s, e));
    doc.seq = doc.seq.slice(0, s) + mid + doc.seq.slice(e);
    for (const f of doc.features) {
      if (f.locs.every(([a, b]) => a >= s && b <= e)) {
        f.locs = f.locs.map(([a, b]) => [s + e - b, s + e - a]).reverse(); f.strand = -f.strand;
      }
    }
    doc.anchor = s; doc.caret = e;
  });
}
function reverseComplementAll(doc) { reverseComplementRange(doc, 0, doc.seq.length); doc.anchor = doc.caret = 0; }

function setOrigin(doc, p) {
  if (!App.editing) return lockedToast();
  const n = doc.seq.length; if (!doc.circular || p <= 0 || p >= n) return;
  mutate(doc, () => {
    doc.seq = doc.seq.slice(p) + doc.seq.slice(0, p);
    for (const f of doc.features) {
      const out = [];
      for (const [a, b] of f.locs) {
        if (a >= p) out.push([a - p, b - p]);
        else if (b <= p) out.push([a + n - p, b + n - p]);
        else { out.push([a + n - p, n]); out.push([0, b - p]); }
      }
      f.locs = out;
    }
    doc.anchor = doc.caret = 0;
  });
}

function addFeature(doc, f, select = true) {
  mutate(doc, () => { const nf = normFeature(f); doc.features.push(nf); if (select) { const [a, b] = featBounds(nf); doc.anchor = a; doc.caret = b; doc.selFid = nf.id; } });
}
function removeFeature(doc, id) { mutate(doc, () => { doc.features = doc.features.filter(f => f.id !== id); doc.selFid = null; }); }
function updateFeature(doc, id, patch) { mutate(doc, () => { const f = doc.features.find(x => x.id === id); if (f) Object.assign(f, patch); }); }

function runDetection(doc, quiet = false) {
  const found = detectFeatures(doc, fullLibrary());
  if (!found.length) { if (!quiet) toast('No new common features found'); return 0; }
  mutate(doc, () => { for (const f of found) doc.features.push(normFeature(f)); });
  if (!quiet) toast(`Added ${found.length} feature${found.length > 1 ? 's' : ''}`);
  return found.length;
}

function addFeatures(doc, feats) {
  mutate(doc, () => { for (const f of feats) doc.features.push(normFeature(f)); });
}

/* primers whose 5' end does not match the sequence (e.g. attB tails): quals.unmatched5 = number of non-matching 5' nt; returns the
   stretch(es) next to the matching part where that overhang is drawn (as a zigzag) */
function hangSegs(f, n) {
  const k = +((f.quals || {}).unmatched5 || [])[0] || 0; if (!k || !f.strand || !f.locs.length || k >= n) return [];
  const x = f.strand === 1 ? f.locs[0][0] - k : f.locs[f.locs.length - 1][1], y = x + k;
  if (x < 0) return [[x + n, n], [0, y]].filter(s => s[1] > s[0]);
  if (y > n) return [[x, n], [0, y - n]].filter(s => s[1] > s[0]);
  return [[x, y]];
}
function visibleFeatures(doc) {
  const list = App.settings.showFeatures ? doc.features.slice() : [];
  if (App.settings.showOrfs) list.push(...findORFs(doc, App.settings.orfMin));
  return list;
}

/* ----- the document list ----- */
function addDoc(doc, { detect = false } = {}) {
  const d = makeDoc(doc);
  if (detect && App.settings.autoDetect && !d.features.length) {
    const found = detectFeatures(d, fullLibrary());
    d.features = found.map(normFeature);
  }
  App.docs.push(d); App.cur = d;
  requestUpdate(true);
  return d;
}
function closeDoc(doc) {
  if (doc.dirty && !confirm(`"${doc.name}" has unsaved changes. Close anyway?`)) return;
  const i = App.docs.indexOf(doc); App.docs.splice(i, 1);
  if (App.cur === doc) App.cur = App.docs[Math.min(i, App.docs.length - 1)] || null;
  requestUpdate(true);
}


/* warn (and later highlight) when input contains anything other than A/C/G/T */
function warnNonATGC(text, where) {
  const st = nonATGCStats(text);
  if (st.count) toast(`⚠ ${st.count.toLocaleString('en-US')} non-ATGC character${st.count > 1 ? 's' : ''} ${where} (${st.detail}) – highlighted in red in the sequence view`, 8000, 'warn');
  return st.count;
}


/* ----- case conversion (length-preserving, so features stay put) ----- */
function changeCase(doc, mode) {
  if (!App.editing) return lockedToast();
  const r = selRange(doc); if (!r) return toast('Select some sequence first');
  const t = doc.seq.slice(r[0], r[1]);
  const nt = mode === 'upper' ? t.toUpperCase() : mode === 'lower' ? t.toLowerCase() : t.replace(/[a-z]/gi, c => (c === c.toUpperCase() ? c.toLowerCase() : c.toUpperCase()));
  if (nt === t) return toast('Nothing to change');
  mutate(doc, () => { doc.seq = doc.seq.slice(0, r[0]) + nt + doc.seq.slice(r[1]); doc.anchor = r[0]; doc.caret = r[1]; });
}

/* ----- translations to the clipboard ----- */
async function copyTranslation(doc, reverse) {
  const r = selRange(doc); if (!r) return toast('Select some sequence first');
  let s = seqU(doc).slice(r[0], r[1]); if (reverse) s = revcomp(s);
  const p = translate(s);
  if (!p) return toast('The selection is shorter than one codon');
  await copyToSystem(p);
  toast(`Copied ${p.length} aa (${reverse ? 'reverse strand' : 'forward strand'}, starting at the first selected base)`);
}
async function copyFeatureTranslation(doc, f) {
  const p = cdsCodons(doc, f).map(c => c.aa).join('').replace(/\*$/, '');
  if (!p) return toast('Nothing to translate');
  await copyToSystem(p);
  toast(`Copied the ${p.length}-aa translation of “${f.name}”`);
}
