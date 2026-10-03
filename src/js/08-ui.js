'use strict';
/* ---------- UI: update loop, panels, dialogs, menus, file I/O, events ---------- */

const UI = { fFilter: '', eFilter: '', eMode: 'cutters', find: { q: '', res: [], idx: -1 } };

/* ===== update loop ===== */
let _pf = false, _pending = false;
function requestUpdate(full) {
  _pf = _pf || full;
  if (_pending) return;
  _pending = true;
  requestAnimationFrame(() => { _pending = false; const f = _pf; _pf = false; update(f); });
}

function renderMap(doc) {
  const pane = $('#mapPane'), host = $('#mapHost');
  pane.classList.toggle('lin', !doc.circular);
  if (doc.circular) renderCircular(doc, host); else renderLinear(doc, host);
}

function update(full) {
  const doc = App.cur;
  document.body.classList.toggle('empty', !doc);
  renderTabs();
  if (!doc) { document.title = 'Plasmid Viewer'; return; }
  document.title = (doc.dirty ? '• ' : '') + doc.name + ' – Plasmid Viewer';
  $('#viewport').className = 'v-' + App.settings.view;
  if (App.settings.view !== 'seq') renderMap(doc);
  if (App.settings.view !== 'map') renderSeq(doc, false);
  renderMini(doc);
  if (full) renderSide(doc); else updateSideSel(doc);
  renderToolbarState(doc);
  renderStatus(doc);
}

function renderTabs() {
  const t = $('#tabs'); t.innerHTML = '';
  for (const d of App.docs) {
    t.append(el('div', { class: 'tab' + (d === App.cur ? ' on' : ''), onclick: () => { App.cur = d; requestUpdate(true); } },
      el('span', { class: 'tn' }, (d.dirty ? '● ' : '') + d.name),
      el('button', { class: 'tx', title: 'Close', onclick: e => { e.stopPropagation(); closeDoc(d); } }, '×')));
  }
}

function renderToolbarState(doc) {
  $$('[data-act^="view-"]').forEach(b => b.classList.toggle('on', b.dataset.act === 'view-' + App.settings.view));
  const eb = $('#editBtn'); eb.classList.toggle('on', App.editing); eb.querySelector('b').textContent = App.editing ? '✎' : '🔒'; eb.querySelector('span').textContent = App.editing ? 'Seq edit on' : 'Seq edit off';
  document.body.classList.toggle('editing', App.editing);
  $('[data-act="undo"]').disabled = !doc.undo.length;
  $('[data-act="redo"]').disabled = !doc.redo.length;
  $('#enzMode').value = App.settings.enzMode;
}

function renderStatus(doc) {
  const sel = selRange(doc), n = doc.seq.length; let left;
  if (sel) {
    const s = doc.seq.slice(sel[0], sel[1]);
    left = `Selected ${sel[0] + 1}..${sel[1]} · ${fmtBp(s.length)} · GC ${gcPercent(s).toFixed(1)}% · Tm ${meltingTemp(s).toFixed(1)}°C`;
  } else left = `Cursor after base ${doc.caret.toLocaleString('en-US')}`;
  if (UI.find.q) { const h = UI.find.res[UI.find.idx]; left += ` · Find: ${UI.find.res.length ? `${UI.find.idx >= 0 ? UI.find.idx + 1 : '–'} of ${UI.find.res.length}` : 'no matches'}${h ? (h.strand === 1 ? ' · forward strand (top)' : ' · reverse strand (bottom)') : ''}`; }
  $('#stLeft').textContent = left;
  $('#stRight').textContent = `${doc.circular ? 'Circular' : 'Linear'} · ${fmtBp(n)} · ${doc.features.length} features · GC ${gcPercent(doc.seq).toFixed(1)}%`;
}

/* ===== side panel ===== */
function renderSide(doc) {
  $$('#sideTabs button').forEach(b => b.classList.toggle('on', b.dataset.tab === App.settings.sideTab));
  const body = $('#sideBody'); body.innerHTML = '';
  const tab = App.settings.sideTab;
  if (tab === 'features') {
    body.append(el('div', { class: 'sidehead' },
      el('input', { type: 'search', placeholder: 'Filter features…', value: UI.fFilter, oninput: e => { UI.fFilter = e.target.value; renderFeatList(App.cur); } }),
      el('button', { class: 'btn sm', onclick: () => openFeatureDialog(App.cur), title: 'Add a feature from the current selection' }, '+ Add')),
      el('div', { id: 'featList', class: 'list' }));
    renderFeatList(doc);
  } else if (tab === 'enzymes') {
    body.append(el('div', { class: 'sidehead' },
      el('input', { type: 'search', placeholder: 'Filter enzymes…', value: UI.eFilter, oninput: e => { UI.eFilter = e.target.value; renderEnzList(App.cur); } }),
      el('select', { onchange: e => { UI.eMode = e.target.value; renderEnzList(App.cur); }, value: UI.eMode },
        el('option', { value: 'cutters' }, 'Cutters'), el('option', { value: 'unique' }, 'Unique'), el('option', { value: 'gg' }, 'Golden Gate'), el('option', { value: 'selected' }, 'My selected'), el('option', { value: 'none' }, 'Non-cutters'), el('option', { value: 'all' }, 'All')),
      el('button', { class: 'btn sm', title: 'Add a custom enzyme', onclick: openEnzymeDialog }, '+')),
      el('div', { id: 'enzList', class: 'list' }));
    renderEnzList(doc);
  } else renderInfo(doc, body);
}

function renderFeatList(doc) {
  const host = $('#featList'); if (!host) return;
  const q = UI.fFilter.toLowerCase();
  const list = doc.features.map(f => ({ f, b: featBounds(f) })).filter(o => !q || o.f.name.toLowerCase().includes(q) || o.f.type.toLowerCase().includes(q)).sort((a, b) => a.b[0] - b.b[0]);
  host.innerHTML = '';
  if (!list.length) host.append(el('div', { class: 'empty-note' }, doc.features.length ? 'No matches.' : 'No features yet. Use “Detect features”, or select a region and click “+ Add”.'));
  for (const { f, b } of list) {
    const len = f.locs.reduce((s, l) => s + l[1] - l[0], 0);
    host.append(el('div', { class: 'item' + (doc.selFid === f.id ? ' on' : ''), 'data-fid': f.id, onclick: () => selectFeature(doc, f), ondblclick: () => openFeatureDialog(doc, f) },
      el('span', { class: 'sw', style: `background:${f.color}` }),
      el('div', { class: 'fn' }, el('b', {}, f.name), el('small', {}, `${f.type} · ${b[0] + 1}..${b[1]} · ${len.toLocaleString('en-US')} bp ${f.strand === 1 ? '→' : f.strand === -1 ? '←' : ''}`))));
  }
}

function renderEnzList(doc) {
  const host = $('#enzList'); if (!host) return;
  const q = UI.eFilter.toLowerCase(), res = analyzeEnzymes(doc), shown = new Set(App.settings.enzShow);
  host.innerHTML = '';
  let count = 0;
  for (const [name, r] of res) {
    const k = r.cuts.length;
    if (q && !name.toLowerCase().includes(q)) continue;
    if (UI.eMode === 'cutters' && !k) continue;
    if (UI.eMode === 'unique' && k !== 1) continue;
    if (UI.eMode === 'gg' && !GOLDEN_GATE.has(name)) continue;
    if (UI.eMode === 'selected' && !shown.has(name)) continue;
    if (UI.eMode === 'none' && k) continue;
    if (++count > 400) break;
    const e = r.enzyme;
    host.append(el('div', { class: 'item enz', 'data-enz': name, onclick: ev => { if (ev.target.tagName === 'INPUT') return; if (r.cuts[0]) setSel(doc, r.cuts[0].s, Math.min(r.cuts[0].e, doc.seq.length)); } },
      el('input', { type: 'checkbox', checked: shown.has(name) || undefined, title: 'Add to “My selected enzymes”', onchange: ev => { const s = new Set(App.settings.enzShow); ev.target.checked ? s.add(name) : s.delete(name); App.settings.enzShow = [...s]; if (ev.target.checked && App.settings.enzMode !== 'custom') { App.settings.enzMode = 'custom'; toast('Showing “My selected enzymes”'); } saveSettings(); requestUpdate(true); } }),
      el('div', { class: 'fn' }, el('b', {}, name), el('small', {}, (GOLDEN_GATE.has(name) ? 'Golden Gate · ' : '') + e.def.replace(/_/g, '') + (k ? ' · ' + r.cuts.slice(0, 6).map(c => c.top).join(', ') + (k > 6 ? '…' : '') : ''))),
      el('span', { class: 'badge' + (k === 1 ? ' u' : '') }, k)));
  }
  if (!count) host.append(el('div', { class: 'empty-note' }, 'Nothing to show.'));
}

function updateSideSel(doc) {
  $$('#featList .item').forEach(i => i.classList.toggle('on', +i.dataset.fid === doc.selFid));
}

function renderInfo(doc, body) {
  const sel = selRange(doc), s = sel ? doc.seq.slice(sel[0], sel[1]) : '';
  const row = (k, v) => el('div', { class: 'kv' }, el('span', {}, k), el('b', {}, v));
  body.append(el('div', { class: 'info' },
    el('label', {}, 'Name', el('input', { value: doc.name, onchange: e => { doc.name = e.target.value.trim() || 'Untitled'; doc.dirty = true; requestUpdate(true); } })),
    el('label', {}, 'Topology', el('select', { value: doc.circular ? 'c' : 'l', onchange: e => { mutate(doc, () => { doc.circular = e.target.value === 'c'; }); } }, el('option', { value: 'c' }, 'Circular'), el('option', { value: 'l' }, 'Linear'))),
    row('Length', fmtBp(doc.seq.length)), row('GC content', gcPercent(doc.seq).toFixed(1) + '%'), row('Features', String(doc.features.length)),
    row('Unique cutters', String([...analyzeEnzymes(doc).values()].filter(r => r.cuts.length === 1).length)),
    el('h4', {}, 'Selection'),
    sel ? el('div', {}, row('Range', `${sel[0] + 1}..${sel[1]}`), row('Length', fmtBp(s.length)), row('GC', gcPercent(s).toFixed(1) + '%'), row('Tm (basic)', meltingTemp(s).toFixed(1) + ' °C'),
      el('div', { class: 'mono wrap' }, el('small', {}, 'Translation (frame 1)'), el('div', {}, translate(s).slice(0, 400) || '—'))) : el('div', { class: 'empty-note' }, 'Nothing selected.'),
    el('h4', {}, 'Tools'),
    el('div', { class: 'btnrow' },
      el('button', { class: 'btn', onclick: () => openInsertDialog(doc) }, 'Insert sequence…'),
      el('button', { class: 'btn', onclick: () => openGotoDialog(doc) }, 'Go to…'),
      el('button', { class: 'btn', onclick: () => reverseComplementAll(doc) }, 'Reverse-complement all'),
      el('button', { class: 'btn', disabled: !doc.circular || undefined, onclick: () => setOrigin(doc, doc.caret) }, 'Set origin at cursor')),
    el('h4', {}, 'Display'),
    el('label', { class: 'chk' }, el('input', { type: 'checkbox', checked: App.settings.showOrfs || undefined, onchange: e => { App.settings.showOrfs = e.target.checked; saveSettings(); requestUpdate(true); } }), 'Show open reading frames (ORFs)'),
    el('label', {}, 'Minimum ORF length (aa)', el('input', { type: 'number', min: 20, max: 2000, value: App.settings.orfMin, onchange: e => { App.settings.orfMin = clamp(+e.target.value || 100, 20, 2000); saveSettings(); requestUpdate(true); } })),
    el('label', { class: 'chk' }, el('input', { type: 'checkbox', checked: App.settings.showTrans || undefined, onchange: e => { App.settings.showTrans = e.target.checked; saveSettings(); requestUpdate(true); } }), 'Show amino acids inside CDS bars'),
    el('label', { class: 'chk' }, el('input', { type: 'checkbox', checked: App.settings.autoDetect || undefined, onchange: e => { App.settings.autoDetect = e.target.checked; saveSettings(); } }), 'Auto-detect common features when opening un-annotated files')));
}

/* ===== selection helpers ===== */
function selectFeature(doc, f) {
  const segs = f.locs, mono = segs.every((l, i) => !i || l[0] >= segs[i - 1][1]);
  const [a, b] = mono ? featBounds(f) : segs[0];
  setSel(doc, a, b, f.id);
}

/* ===== modal / menu ===== */
function openModal(title, content, buttons = [], { wide = false } = {}) {
  const back = el('div', { class: 'backdrop' });
  const close = () => { back.remove(); document.removeEventListener('keydown', onKey, true); };
  const foot = el('div', { class: 'mfoot' }, buttons.map(b => el('button', { class: 'btn' + (b.primary ? ' primary' : '') + (b.danger ? ' danger' : ''), onclick: () => { if (b.action && b.action() === false) return; close(); } }, b.label)));
  const box = el('div', { class: 'modal' + (wide ? ' wide' : ''), role: 'dialog' }, el('div', { class: 'mhead' }, el('h3', {}, title), el('button', { class: 'tx', onclick: close }, '×')), el('div', { class: 'mbody' }, content), buttons.length ? foot : null);
  back.append(box); $('#modalRoot').append(back);
  const onKey = e => {
    if (e.key === 'Escape') { e.stopPropagation(); close(); }
    else if (e.key === 'Enter' && e.target.tagName !== 'TEXTAREA' && buttons.some(b => b.primary)) { e.preventDefault(); e.stopPropagation(); foot.querySelector('.primary').click(); }
  };
  document.addEventListener('keydown', onKey, true);
  back.addEventListener('mousedown', e => { if (e.target === back) close(); });
  const first = box.querySelector('input,textarea,select'); if (first) setTimeout(() => first.focus(), 30);
  return { close };
}
const field = (label, input, hint) => el('label', { class: 'field' }, el('span', {}, label), input, hint ? el('small', {}, hint) : null);

function showMenu(x, y, items) {
  closeMenu();
  const m = el('div', { class: 'menu', id: 'ctxMenu' });
  for (const it of items) {
    if (it === '-') { m.append(el('div', { class: 'sep' })); continue; }
    m.append(el('button', { class: 'mi', disabled: it.disabled || undefined, onclick: () => { closeMenu(); it.action && it.action(); } }, el('span', { class: 'ck' }, it.checked ? '✓' : ''), it.label, it.hint ? el('em', {}, it.hint) : null));
  }
  $('#menuRoot').append(m);
  const r = m.getBoundingClientRect();
  m.style.left = Math.min(x, innerWidth - r.width - 6) + 'px'; m.style.top = Math.min(y, innerHeight - r.height - 6) + 'px';
}
function closeMenu() { const m = $('#ctxMenu'); if (m) m.remove(); }
document.addEventListener('mousedown', e => { if (!e.target.closest('#ctxMenu')) closeMenu(); }, true);
window.addEventListener('blur', closeMenu);

const isMac = /Mac|iPhone|iPad/.test(navigator.platform);
const MOD = isMac ? '⌘' : 'Ctrl+';

function contextMenu(e) {
  const doc = App.cur; if (!doc) return;
  e.preventDefault();
  const fe = e.target.closest('[data-fid]');
  if (fe) { const f = doc.features.find(x => x.id === +fe.dataset.fid); if (f) selectFeature(doc, f); }
  const sel = selRange(doc), f = fe && doc.features.find(x => x.id === +fe.dataset.fid), ro = !App.editing;
  const items = [
    { label: 'Cut', hint: MOD + 'X', disabled: !sel || ro, action: () => doCopy(doc, true) },
    { label: 'Copy', hint: MOD + 'C', disabled: !sel, action: () => doCopy(doc) },
    { label: 'Copy reverse complement', disabled: !sel, action: () => doCopyRC(doc) },
    { label: 'Paste', hint: MOD + 'V', disabled: ro, action: () => doPaste(doc) },
    { label: 'Paste reverse complement', disabled: ro, action: () => doPaste(doc, true) },
    { label: 'Delete', hint: '⌫', disabled: !sel || ro, action: () => deleteSelection(doc) }, '-',
    { label: 'Add feature from selection…', disabled: !sel, action: () => openFeatureDialog(doc) },
  ];
  if (f) items.push({ label: `Edit “${f.name}”…`, action: () => openFeatureDialog(doc, f) }, { label: `Delete feature “${f.name}”`, action: () => removeFeature(doc, f.id) });
  items.push('-', { label: 'Reverse-complement selection', disabled: !sel || ro, action: () => reverseComplementRange(doc, ...sel) },
    { label: 'Set origin here', disabled: !doc.circular || ro, action: () => setOrigin(doc, doc.caret) },
    { label: 'Insert sequence here…', disabled: ro, action: () => openInsertDialog(doc) },
    { label: 'Select all', hint: MOD + 'A', action: () => setSel(doc, 0, doc.seq.length) });
  showMenu(e.clientX, e.clientY, items);
}

/* ===== dialogs ===== */
function parseRanges(text, n) {
  const locs = [];
  for (const part of text.split(/[,;]+/).map(s => s.trim()).filter(Boolean)) {
    const m = part.match(/^(\d+)\s*(?:\.\.|-|–|:)\s*(\d+)$/) || part.match(/^(\d+)$/);
    if (!m) return null;
    const a = +m[1], b = m[2] ? +m[2] : a;
    if (a < 1 || b < a || b > n) return null;
    locs.push([a - 1, b]);
  }
  return locs.length ? locs : null;
}

function openFeatureDialog(doc, f) {
  const sel = selRange(doc);
  if (!f && !sel) return toast('Select a region of the sequence first');
  const isNew = !f;
  const rangesDefault = f ? f.locs.map(l => `${l[0] + 1}..${l[1]}`).join(', ') : `${sel[0] + 1}..${sel[1]}`;
  const name = el('input', { value: f ? f.name : '', placeholder: 'e.g. My promoter' });
  const type = el('select', {}, FEATURE_TYPES.map(t => el('option', { value: t }, t)));
  type.value = f ? f.type : 'misc_feature';
  const strand = el('select', { value: String(f ? f.strand : 1) }, el('option', { value: '1' }, 'Forward →'), el('option', { value: '-1' }, 'Reverse ←'), el('option', { value: '0' }, 'None'));
  const color = el('input', { type: 'color', value: f ? f.color : typeColor('misc_feature') });
  let colorTouched = !!f;
  color.addEventListener('input', () => { colorTouched = true; });
  type.addEventListener('change', () => { if (!colorTouched) color.value = typeColor(type.value); });
  const ranges = el('input', { value: rangesDefault });
  const note = el('textarea', { rows: 2, placeholder: 'Optional note' }); note.value = f ? ((f.quals.note || []).join('\n')) : '';
  const lib = el('input', { type: 'checkbox' });
  const body = el('div', { class: 'form' }, field('Name', name), el('div', { class: 'two' }, field('Type', type), field('Direction', strand)), el('div', { class: 'two' }, field('Colour', color), field('Range(s), 1-based', ranges, 'e.g. 120..845 or 900..1000, 1..40')), field('Note', note),
    el('label', { class: 'chk' }, lib, 'Add to feature library (so it is auto-detected in future)'));
  const buttons = [{ label: 'Cancel' }];
  if (f) buttons.unshift({ label: 'Delete', danger: true, action: () => removeFeature(doc, f.id) });
  buttons.push({
    label: isNew ? 'Add feature' : 'Save', primary: true, action: () => {
      const locs = parseRanges(ranges.value, doc.seq.length);
      if (!locs) { toast('Invalid range – use 1-based start..end'); return false; }
      if (!name.value.trim()) { toast('Please give the feature a name'); return false; }
      const quals = { ...(f ? f.quals : {}) }; const nt = note.value.trim(); if (nt) quals.note = nt.split('\n'); else delete quals.note;
      const data = { name: name.value.trim(), type: type.value, strand: +strand.value, color: color.value, locs, quals };
      if (lib.checked) {
        let s = locs.map(([a, b]) => doc.seq.slice(a, b)).join(''); if (data.strand === -1) s = revcomp(s);
        if (s.length < 6) toast('Too short for the library (min 6 bp)');
        else { const u = userLibrary().filter(x => x.name.toLowerCase() !== data.name.toLowerCase()); u.push({ name: data.name, type: data.type, color: data.color, seq: s, dir: data.strand !== 0 }); saveUserLibrary(u); toast(`“${data.name}” added to the feature library`); }
      }
      if (isNew) addFeature(doc, data); else updateFeature(doc, f.id, data);
    },
  });
  openModal(isNew ? 'Add feature' : 'Edit feature', body, buttons);
}

function openInsertDialog(doc) {
  if (!App.editing) return lockedToast();
  const pos = el('input', { type: 'number', min: 0, max: doc.seq.length, value: selRange(doc) ? selRange(doc)[0] : doc.caret });
  const ta = el('textarea', { rows: 6, class: 'mono', placeholder: 'Paste or type DNA (FASTA allowed)…' });
  const rc = el('input', { type: 'checkbox' });
  openModal('Insert sequence', el('div', { class: 'form' }, field('Insert after base', pos, '0 = at the very start'), field('Sequence', ta), el('label', { class: 'chk' }, rc, 'Insert reverse complement')),
    [{ label: 'Cancel' }, { label: 'Insert', primary: true, action: () => { let s = cleanSeq(ta.value); if (!s) { toast('No valid DNA found'); return false; } if (rc.checked) s = revcomp(s); const p = clamp(+pos.value || 0, 0, doc.seq.length); editReplace(doc, p, p, s); } }]);
}

function openGotoDialog(doc) {
  const inp = el('input', { placeholder: 'e.g. 1500 or 1500..1800' });
  openModal('Go to', field('Position or range (1-based)', inp), [{ label: 'Cancel' }, { label: 'Go', primary: true, action: () => {
    const l = parseRanges(inp.value, doc.seq.length); if (!l) { toast('Invalid position'); return false; }
    if (inp.value.match(/\.\.|-/)) setSel(doc, l[0][0], l[0][1]); else setSel(doc, l[0][1], l[0][1]);
  } }]);
}

function openNewSeqDialog() {
  const name = el('input', { value: 'New plasmid' }), ta = el('textarea', { rows: 8, class: 'mono', placeholder: 'Paste DNA sequence (raw, FASTA or GenBank text)…' });
  const topo = el('select', {}, el('option', { value: 'c' }, 'Circular'), el('option', { value: 'l' }, 'Linear'));
  openModal('New from sequence', el('div', { class: 'form' }, field('Name', name), field('Topology', topo), field('Sequence', ta)), [{ label: 'Cancel' }, { label: 'Create', primary: true, action: () => {
    try {
      let d;
      if (/^LOCUS/m.test(ta.value)) d = parseGenBank(ta.value); else d = { name: name.value.trim() || 'New plasmid', seq: cleanSeq(ta.value), circular: topo.value === 'c', features: [], meta: {} };
      if (!d.seq) { toast('No valid DNA found'); return false; }
      addDoc(d, { detect: true });
    } catch (e) { toast(e.message); return false; }
  } }]);
}

function openEnzymeDialog() {
  const name = el('input', { placeholder: 'e.g. MyEnzI' }), def = el('input', { placeholder: 'G^AATTC   (use ^ for the cut; _ for bottom-strand cut if non-palindromic)' });
  openModal('Add restriction enzyme', el('div', { class: 'form' }, field('Name', name), field('Site', def, 'IUPAC codes allowed (R, Y, N, …)')), [{ label: 'Cancel' }, { label: 'Add', primary: true, action: () => {
    const e = compileEnzyme(name.value.trim(), def.value.trim().toUpperCase());
    if (!e || !name.value.trim()) { toast('Invalid enzyme definition'); return false; }
    const u = store.get('userEnzymes', []).filter(x => x.name !== e.name); u.push({ name: e.name, def: e.def }); store.set('userEnzymes', u);
    buildEnzymes(); if (App.cur) App.cur._cache = {}; requestUpdate(true); toast(`Added ${e.name}`);
  } }]);
}

function openLibraryDialog() {
  const box = el('div', { class: 'libdlg' });
  const render = () => {
    box.innerHTML = '';
    const lib = fullLibrary();
    box.append(el('div', { class: 'small-note' }, `${lib.length} entries · ${lib.filter(l => l.user).length} yours. Sequences are matched on both strands (long ones tolerate ~4% mismatches); protein motifs are matched in all six frames.`));
    const list = el('div', { class: 'list tall' });
    for (const l of lib) {
      list.append(el('div', { class: 'item' }, el('span', { class: 'sw', style: `background:${l.color || typeColor(l.type)}` }),
        el('div', { class: 'fn' }, el('b', {}, l.name), el('small', {}, `${l.type} · ${l.seq ? l.seq.length + ' bp' : l.aa.length + ' aa (protein)'}${l.user ? ' · yours' : ' · built-in'}`)),
        l.user ? el('button', { class: 'tx', title: 'Remove', onclick: () => { saveUserLibrary(userLibrary().filter(x => x.name !== l.name)); render(); } }, '🗑') : null));
    }
    box.append(list);
    const name = el('input', { placeholder: 'Name' }), type = el('select', {}, FEATURE_TYPES.map(t => el('option', { value: t }, t))), color = el('input', { type: 'color', value: typeColor('misc_feature') });
    const kind = el('select', {}, el('option', { value: 'dna' }, 'DNA'), el('option', { value: 'aa' }, 'Protein'));
    const seq = el('textarea', { rows: 3, class: 'mono', placeholder: 'DNA (ACGT…) or protein motif' });
    type.value = 'misc_feature';
    box.append(el('h4', {}, 'Add a new feature'), el('div', { class: 'form' }, el('div', { class: 'two' }, field('Name', name), field('Type', type)), el('div', { class: 'two' }, field('Colour', color), field('Sequence kind', kind)), field('Sequence', seq),
      el('div', { class: 'btnrow' }, el('button', { class: 'btn primary', onclick: () => {
        if (!name.value.trim()) return toast('Name required');
        const e = { name: name.value.trim(), type: type.value, color: color.value };
        if (kind.value === 'dna') { e.seq = cleanSeq(seq.value); if (e.seq.length < 6) return toast('Need at least 6 bp'); } else { e.aa = seq.value.toUpperCase().replace(/[^A-Z]/g, ''); if (e.aa.length < 3) return toast('Need at least 3 residues'); }
        saveUserLibrary(userLibrary().filter(x => x.name.toLowerCase() !== e.name.toLowerCase()).concat([e])); render(); toast('Added to library');
      } }, 'Add to library'),
      el('button', { class: 'btn', onclick: () => downloadBlob('feature-library.json', JSON.stringify(userLibrary(), null, 2), 'application/json') }, 'Export mine…'),
      el('button', { class: 'btn', onclick: () => { const i = el('input', { type: 'file', accept: '.json' }); i.onchange = async () => { try { const arr = JSON.parse(await i.files[0].text()); const cur = userLibrary(); for (const x of arr) if (x.name && (x.seq || x.aa)) { const k = cur.findIndex(c => c.name.toLowerCase() === x.name.toLowerCase()); if (k >= 0) cur[k] = x; else cur.push(x); } saveUserLibrary(cur); render(); toast('Library imported'); } catch (e) { toast('Could not read that file'); } }; i.click(); } }, 'Import…'))));
  };
  render();
  openModal('Feature library', box, [{ label: 'Close', primary: true }], { wide: true });
}

function openHelp() {
  const rows = [[MOD + 'O', 'Open file(s)'], [MOD + 'S', 'Save as GenBank'], [MOD + 'Z / ⇧' + MOD + 'Z', 'Undo / redo'], [MOD + 'C / X / V', 'Copy / cut / paste (features travel with the sequence)'], [MOD + 'A', 'Select all'], [MOD + 'F', 'Find (both strands, IUPAC ok)'], [MOD + 'G', 'Go to position'],
    ['A C G T N …', 'Type bases at the cursor (sequence view)'], ['⌫ / Del', 'Delete selection or base'], ['← → ↑ ↓ (⇧)', 'Move cursor / extend selection'], ['Double-click feature', 'Edit feature'], ['Right-click', 'Context menu']];
  openModal('Shortcuts & tips', el('div', {}, el('table', { class: 'kbd' }, rows.map(r => el('tr', {}, el('td', {}, el('code', {}, r[0])), el('td', {}, r[1])))),
    el('p', { class: 'small-note' }, 'Supported files: GenBank (.gb/.gbk), FASTA, SnapGene .dna (read-only). Save writes GenBank, which SnapGene, Benchling, ApE and others can open. Copying a region copies fully-contained features too; pasting adds them at the new position.')), [{ label: 'Close', primary: true }]);
}

/* ===== find ===== */
function runFind(dir) {
  const doc = App.cur; if (!doc) return;
  const q = $('#find').value;
  if (q !== UI.find.q) { UI.find = { q, res: findAll(doc, q), idx: -1 }; }
  const f = UI.find; if (!f.res.length) { renderStatus(doc); return; }
  const base = selRange(doc) ? (dir > 0 ? selRange(doc)[1] - 1 : selRange(doc)[0]) : doc.caret;
  let idx;
  if (dir > 0) { idx = f.res.findIndex(r => r.s > base); if (idx < 0) idx = 0; }
  else { idx = -1; f.res.forEach((r, i) => { if (r.s < base) idx = i; }); if (idx < 0) idx = f.res.length - 1; }
  f.idx = idx; setSel(doc, f.res[idx].s, f.res[idx].e);
}

/* ===== file I/O ===== */
async function openFiles(files) {
  for (const file of files) {
    try {
      const buf = await file.arrayBuffer();
      const d = parseAny(buf, file.name);
      const doc = addDoc(d, { detect: true });
      if (doc.features.length && !d.features.length) doc.dirty = false;
    } catch (e) { toast(`${file.name}: ${e.message}`, 5000); }
  }
}
function pickFiles() {
  const i = el('input', { type: 'file', multiple: true, accept: '.gb,.gbk,.genbank,.gbff,.fasta,.fa,.fna,.fas,.dna,.txt,.seq' });
  i.onchange = () => openFiles(Array.from(i.files)); i.click();
}

async function saveDoc(doc, saveAs = false) {
  const text = writeGenBank(doc), fname = doc.name.replace(/[^\w.-]+/g, '_') + '.gb';
  try {
    if (window.plasmidNative) {
      const p = await window.plasmidNative.saveFile(fname, text, saveAs ? null : doc.path);
      if (!p) return; doc.path = p;
    } else if (window.showSaveFilePicker) {
      if (!doc.handle || saveAs) doc.handle = await window.showSaveFilePicker({ suggestedName: fname, types: [{ description: 'GenBank', accept: { 'text/plain': ['.gb', '.gbk'] } }] });
      const w = await doc.handle.createWritable(); await w.write(text); await w.close();
      doc.name = doc.name;
    } else downloadBlob(fname, text);
    doc.dirty = false; toast('Saved ' + fname); requestUpdate(true);
  } catch (e) { if (e.name !== 'AbortError') { doc.handle = null; downloadBlob(fname, text); doc.dirty = false; requestUpdate(true); } }
}
function exportDoc(doc, fmt) {
  if (fmt === 'fasta') downloadBlob(doc.name + '.fasta', writeFasta(doc));
  else downloadBlob(doc.name + '.gb', writeGenBank(doc));
}

/* ===== demo plasmid ===== */
function buildDemo() {
  let seed = 12345; const rnd = () => (seed = (seed * 1664525 + 1013904223) >>> 0) / 4294967296;
  const spacer = k => Array.from({ length: k }, () => 'ACGT'[Math.floor(rnd() * 4)]).join('');
  const L = name => BUILTIN_FEATURES.find(f => f.name === name).seq;
  const parts = [spacer(80), L('CMV enhancer'), L('CMV promoter'), spacer(25), 'GAATTCGAGCTCGGTACCCGGGGATCCTCTAGAGTCGACCTGCAGGCATGCAAGCTT', 'GCCACC', L('EGFP').slice(0, -3), 'GGATCCGGC' + 'GACTACAAGGACGACGATGACAAG' + 'CACCATCACCATCACCAC' + 'TAA',
    spacer(60), L('bGH poly(A) signal'), spacer(120), L('f1 ori'), spacer(90), L('SV40 promoter'), spacer(18), L('KanR'), spacer(40), L('SV40 poly(A) signal'), spacer(100),
    revcomp(L('ori')), spacer(70), revcomp(L('AmpR')), spacer(30), revcomp(L('AmpR promoter')), spacer(50)];
  return { name: 'demo_plasmid', seq: parts.join(''), circular: true, features: [], meta: { definition: 'Synthetic demo plasmid (random spacers)' } };
}

/* ===== theme ===== */
function applyTheme() {
  const t = App.settings.theme;
  if (t === 'auto') document.documentElement.removeAttribute('data-theme'); else document.documentElement.setAttribute('data-theme', t);
}

/* ===== actions (toolbar) ===== */
const actions = {
  open: pickFiles, new: openNewSeqDialog, demo: () => addDoc(buildDemo(), { detect: true }),
  save: () => App.cur && saveDoc(App.cur), saveas: () => App.cur && saveDoc(App.cur, true),
  undo: () => App.cur && undo(App.cur), redo: () => App.cur && redo(App.cur),
  'view-map': () => setView('map'), 'view-seq': () => setView('seq'), 'view-split': () => setView('split'),
  detect: () => App.cur && runDetection(App.cur), addfeature: () => App.cur && openFeatureDialog(App.cur), library: openLibraryDialog, help: openHelp,
  'find-next': () => runFind(1), 'find-prev': () => runFind(-1),
  toggleedit: () => { App.editing = !App.editing; toast(App.editing ? 'Sequence editing enabled' : 'Sequence editing disabled (features can still be edited)'); requestUpdate(true); },
  toggleside: () => { App.settings.sideOpen = !App.settings.sideOpen; saveSettings(); applyLayout(); requestUpdate(false); },
  theme: () => { const o = ['auto', 'light', 'dark']; App.settings.theme = o[(o.indexOf(App.settings.theme) + 1) % 3]; saveSettings(); applyTheme(); toast('Theme: ' + App.settings.theme); },
  export: e => { const r = e.target.closest('button').getBoundingClientRect(); showMenu(r.left, r.bottom + 4, [{ label: 'Export GenBank (.gb)', action: () => App.cur && exportDoc(App.cur, 'gb') }, { label: 'Export FASTA (.fasta)', action: () => App.cur && exportDoc(App.cur, 'fasta') }, { label: 'Save as…', action: () => App.cur && saveDoc(App.cur, true) }]); },
  display: e => {
    const r = e.target.closest('button').getBoundingClientRect(), S = App.settings;
    const tog = (k, label) => ({ label, checked: S[k], action: () => { S[k] = !S[k]; saveSettings(); requestUpdate(true); } });
    showMenu(r.left, r.bottom + 4, [tog('showFeatures', 'Features'), tog('showOrfs', 'Open reading frames'), tog('showTrans', 'Amino acids in CDS bars'), '-',
      { label: 'Zoom linear map in', action: () => { S.linZoom = Math.min(S.linZoom * 1.5, 12); saveSettings(); requestUpdate(false); } },
      { label: 'Zoom linear map out', action: () => { S.linZoom = Math.max(S.linZoom / 1.5, 1); saveSettings(); requestUpdate(false); } }]);
  },
};
function setView(v) { App.settings.view = v; saveSettings(); SV.key = ''; requestUpdate(true); }

/* ===== wiring ===== */
function wireMap() {
  const pane = $('#mapPane'); let drag = null;
  pane.addEventListener('mousedown', e => {
    const doc = App.cur; if (!doc || e.button !== 0) return;
    const svg = pane.querySelector('svg'); if (!svg) return;
    const n = doc.seq.length, circ = doc.circular;
    const enz = e.target.closest('[data-enz]');
    if (enz) { setSel(doc, +enz.dataset.s % n, Math.min(+enz.dataset.e, n)); return; }
    const fe = e.target.closest('[data-fid]');
    if (fe) { const f = doc.features.find(x => x.id === +fe.dataset.fid) || findORFs(doc, App.settings.orfMin).find(x => x.id === fe.dataset.fid); if (f) selectFeature(doc, f); return; }
    let pos;
    if (circ) { const c = circularPos(svg, e, n); if (c.r < 120 || c.r > 380) { return; } pos = c.pos; } else pos = linearPos(svg, e, n);
    drag = { start: pos }; setSel(doc, pos, pos, null, false);
    e.preventDefault();
  });
  window.addEventListener('mousemove', e => {
    if (!drag) return; const doc = App.cur, svg = pane.querySelector('svg'); if (!svg) return;
    const n = doc.seq.length, pos = doc.circular ? circularPos(svg, e, n).pos : linearPos(svg, e, n);
    setSel(doc, drag.start, pos, null, false);
  });
  window.addEventListener('mouseup', () => { drag = null; });
  pane.addEventListener('dblclick', e => { const fe = e.target.closest('[data-fid]'), doc = App.cur; if (fe && doc) { const f = doc.features.find(x => x.id === +fe.dataset.fid); if (f) openFeatureDialog(doc, f); } });
  pane.addEventListener('contextmenu', contextMenu);
}

function wireSeq() {
  const scroll = $('#seqScroll'); let drag = null;
  scroll.addEventListener('mousedown', e => {
    const doc = App.cur; if (!doc || e.button !== 0) return;
    scroll.focus({ preventScroll: true }); e.preventDefault();
    const n = doc.seq.length;
    const enz = e.target.closest('[data-enz]');
    if (enz) { setSel(doc, +enz.dataset.s % n, Math.min(+enz.dataset.e, n), null, false); return; }
    const fe = e.target.closest('.fbar');
    if (fe) { const id = fe.dataset.fid; const f = doc.features.find(x => x.id === +id) || findORFs(doc, App.settings.orfMin).find(x => x.id === id); if (f) selectFeature(doc, f); doc._reveal = false; return; }
    if (e.target.closest('.gut')) return;
    const pos = seqPosFromEvent(e);
    if (e.shiftKey) setSel(doc, doc.anchor, pos, null, false); else setSel(doc, pos, pos, null, false);
    drag = true;
  });
  window.addEventListener('mousemove', e => {
    if (!drag || !App.cur) return;
    const r = scroll.getBoundingClientRect();
    if (e.clientY < r.top + 20) scroll.scrollTop -= 18; else if (e.clientY > r.bottom - 20) scroll.scrollTop += 18;
    setSel(App.cur, App.cur.anchor, seqPosFromEvent(e), null, false);
  });
  window.addEventListener('mouseup', () => { drag = null; });
  scroll.addEventListener('dblclick', e => { const fe = e.target.closest('.fbar'), doc = App.cur; if (fe && doc) { const f = doc.features.find(x => x.id === +fe.dataset.fid); if (f) openFeatureDialog(doc, f); } });
  scroll.addEventListener('contextmenu', contextMenu);
}

function inTextField(t) { return t && (['INPUT', 'TEXTAREA', 'SELECT'].includes(t.tagName) || t.isContentEditable); }

let clipTimer = 0;
function wireKeys() {
  document.addEventListener('keydown', e => {
    if (inTextField(e.target) || $('.backdrop')) return;
    const doc = App.cur, mod = e.metaKey || e.ctrlKey, k = e.key;
    if (mod) {
      const lk = k.toLowerCase();
      if (lk === 'o') { e.preventDefault(); pickFiles(); }
      else if (!doc) return;
      else if (lk === 's') { e.preventDefault(); saveDoc(doc, e.shiftKey); }
      else if (lk === 'z') { e.preventDefault(); e.shiftKey ? redo(doc) : undo(doc); }
      else if (lk === 'y') { e.preventDefault(); redo(doc); }
      else if (lk === 'a') { e.preventDefault(); setSel(doc, 0, doc.seq.length, null, false); }
      else if (lk === 'f') { e.preventDefault(); $('#find').focus(); $('#find').select(); }
      else if (lk === 'g') { e.preventDefault(); openGotoDialog(doc); }
      else if ('cxv'.includes(lk) && !e.shiftKey && !e.altKey) {
        // normally the native copy/cut/paste events handle this; this fallback covers environments where they never fire
        clearTimeout(clipTimer);
        clipTimer = setTimeout(() => { if (lk === 'v') doPaste(doc); else doCopy(doc, lk === 'x'); }, 80);
      }
      return;
    }
    if (!doc) return;
    const seqFocus = document.activeElement === $('#seqScroll'), bodyFocus = document.activeElement === document.body || $('#mapPane').contains(document.activeElement);
    if (!seqFocus && !bodyFocus) return;
    const n = doc.seq.length, sel = selRange(doc), bpr = SV.bpr;
    const move = (target) => { target = clamp(target, 0, n); if (e.shiftKey) setSel(doc, doc.anchor, target); else setSel(doc, target, target); };
    switch (k) {
      case 'ArrowLeft': e.preventDefault(); if (sel && !e.shiftKey) setSel(doc, sel[0], sel[0]); else move(doc.caret - 1); return;
      case 'ArrowRight': e.preventDefault(); if (sel && !e.shiftKey) setSel(doc, sel[1], sel[1]); else move(doc.caret + 1); return;
      case 'ArrowUp': e.preventDefault(); move(doc.caret - bpr); return;
      case 'ArrowDown': e.preventDefault(); move(doc.caret + bpr); return;
      case 'Home': e.preventDefault(); move(e.ctrlKey ? 0 : Math.floor(doc.caret / bpr) * bpr); return;
      case 'End': e.preventDefault(); move(e.ctrlKey ? n : Math.min(n, (Math.floor(doc.caret / bpr) + 1) * bpr)); return;
      case 'PageUp': e.preventDefault(); move(doc.caret - bpr * 10); return;
      case 'PageDown': e.preventDefault(); move(doc.caret + bpr * 10); return;
      case 'Escape': if (sel) setSel(doc, doc.caret, doc.caret); return;
      case 'Backspace': e.preventDefault(); deleteSelection(doc, false); return;
      case 'Delete': e.preventDefault(); deleteSelection(doc, true); return;
    }
    if (seqFocus && /^[a-zA-Z]$/.test(k) && /[ACGTUNRYKMSWBDHV]/i.test(k)) { e.preventDefault(); if (!App.editing) return lockedToast(); typeBases(doc, k); }
  });
  document.addEventListener('copy', e => { clearTimeout(clipTimer); const doc = App.cur; if (!doc || inTextField(e.target)) return; const clip = currentClip(doc); if (!clip) return; e.preventDefault(); e.clipboardData.setData('text/plain', clip.text); storeClip(clip); toast(`Copied ${fmtBp(clip.text.length)}${clip.features.length ? ` + ${clip.features.length} feature(s)` : ''}`); });
  document.addEventListener('cut', e => { clearTimeout(clipTimer); const doc = App.cur; if (!doc || inTextField(e.target)) return; const clip = currentClip(doc); if (!clip) return; e.preventDefault(); e.clipboardData.setData('text/plain', clip.text); storeClip(clip); deleteSelection(doc); toast(`Cut ${fmtBp(clip.text.length)}`); });
  document.addEventListener('paste', e => { clearTimeout(clipTimer); const doc = App.cur; if (!doc || inTextField(e.target)) return; e.preventDefault(); const t = e.clipboardData.getData('text/plain'); if (t) pasteClip(doc, lookupClip(t)); });
}

/* ===== hover: red outline on shape + name, and an info pop-up ===== */
function featureById(doc, id) {
  return doc.features.find(x => String(x.id) === String(id)) || findORFs(doc, App.settings.orfMin).find(x => x.id === id);
}
const HOVER_DELAY = 1000;   // ms before the info pop-up appears (the red highlight is immediate)
function wireHover() {
  const tip = $('#tip'); let cur = null, timer = 0, last = null;
  const hide = () => { clearTimeout(timer); if (cur !== null) $$('.hov').forEach(n => n.classList.remove('hov')); cur = null; tip.style.display = 'none'; };
  const move = e => {
    const r = tip.getBoundingClientRect(), enz = tip.classList.contains('enzTip');
    tip.style.right = 'auto';
    tip.style.left = Math.max(4, Math.min(e.clientX + 16, innerWidth - r.width - 8)) + 'px';
    // enzyme pop-up sits to the top-right of the pointer so it never covers the site / cut marks being hovered
    tip.style.top = (enz ? Math.max(4, e.clientY - r.height - 14) : Math.max(4, Math.min(e.clientY + 18, innerHeight - r.height - 8))) + 'px';
  };
  const featureTip = (doc, id) => {
    const f = featureById(doc, id); if (!f) return '';
    const len = f.locs.reduce((s, l) => s + l[1] - l[0], 0), q = f.quals || {};
    const rows = [['Type', f.type], ['Position', f.locs.map(l => `${l[0] + 1}..${l[1]}`).join(', ')], ['Length', fmtBp(len) + (f.type === 'CDS' || f.orf ? ` · ${Math.floor(len / 3) - 1} aa` : '')], ['Direction', f.strand === 1 ? 'Forward →' : f.strand === -1 ? 'Reverse ←' : 'None']];
    const product = (q.product || []).join('; '), note = (q.note || []).join(' ');
    return `<div class="th"><i style="background:${f.color}"></i>${esc(f.name)}</div>` + rows.map(r => `<div class="tr"><span>${r[0]}</span><b>${esc(r[1])}</b></div>`).join('')
      + (product ? `<div class="tn"><u>/product</u> ${esc(product.slice(0, 300))}</div>` : '') + (note ? `<div class="tn"><u>/note</u> ${esc(note.slice(0, 400))}${note.length > 400 ? '…' : ''}</div>` : '');
  };
  document.addEventListener('mouseover', e => {
    const doc = App.cur; if (!doc || !e.target.closest) { hide(); return; }
    last = e;
    const fe = e.target.closest('#viewport [data-fid]'), ee = !fe && e.target.closest('[data-enz]');
    if (!fe && !ee) { if (cur !== null) hide(); return; }
    const key = fe ? 'f' + fe.dataset.fid : 'e' + (ee.dataset.ek || ee.dataset.enz);
    if (key === cur) return;
    hide(); cur = key;
    let html;
    if (fe) { $$('#viewport [data-fid]').forEach(n => { if (n.dataset.fid === fe.dataset.fid) n.classList.add('hov'); }); html = featureTip(doc, fe.dataset.fid); }
    else {
      const nm = ee.dataset.enz;
      $$('#viewport [data-ek]').forEach(n => { if (ee.dataset.ek ? n.dataset.ek === ee.dataset.ek : n.dataset.ek.startsWith(nm + ':')) n.classList.add('hov'); });
      $$('#viewport [data-enz]').forEach(n => { if (n.dataset.enz === nm) n.classList.add('hov'); });
      html = enzymeTipHTML(doc, nm);
    }
    if (!html) return;
    timer = setTimeout(() => {
      tip.classList.toggle('enzTip', !fe);
      tip.innerHTML = html; tip.style.display = 'block';
      move(last);
    }, HOVER_DELAY);
  });
  document.addEventListener('mousemove', e => { last = e; if (cur !== null && tip.style.display === 'block') move(e); });
  document.addEventListener('mouseleave', hide);
  document.addEventListener('mousedown', hide, true);
}

/* ===== nucleotide position under the mouse ===== */
function wirePos() {
  const pt = $('#posTip'), pane = $('#mapPane');
  const hide = () => { pt.style.display = 'none'; };
  document.addEventListener('mousemove', e => {
    const doc = App.cur; if (!doc || !e.target.closest || $('#tip').style.display === 'block' || $('.backdrop')) return hide();
    const n = doc.seq.length; let p = -1;
    const rs = e.target.closest('.rseq');
    if (rs && !e.target.closest('.fbar, .elbl')) {
      const col = Math.floor((e.clientX - rs.getBoundingClientRect().left) / SV.cw), q = +rs.parentElement.dataset.r * SV.bpr + col;
      if (col >= 0 && q < n) p = q;
    } else if (pane.contains(e.target) && !e.target.closest('[data-fid], [data-enz]')) {
      const svg = pane.querySelector('svg');
      if (svg) { if (doc.circular) { const c = circularPos(svg, e, n); if (c.r > 120 && c.r < 380) p = c.base; } else p = linearBase(svg, e, n); }
    }
    if (p < 0) return hide();
    pt.textContent = `${(p + 1).toLocaleString('en-US')} · ${doc.seq[p]}`;
    pt.style.display = 'block';
    pt.style.left = Math.min(e.clientX + 14, innerWidth - pt.offsetWidth - 6) + 'px'; pt.style.top = Math.max(4, e.clientY - 30) + 'px';
  });
  document.addEventListener('mouseleave', hide);
  document.addEventListener('mousedown', hide, true);
}

/* ===== layout: show/hide side panel, draggable splitters ===== */
function applyLayout() {
  const S = App.settings, main = $('#main');
  document.body.classList.toggle('noside', !S.sideOpen);
  $('#sideBtn').classList.toggle('on', S.sideOpen);
  if (S.sideW) main.style.setProperty('--sidew', S.sideW + 'px'); else main.style.removeProperty('--sidew');
  if (S.mapW) main.style.setProperty('--mapw', S.mapW + 'px'); else main.style.removeProperty('--mapw');
}
function wireSplit() {
  const main = $('#main'), S = App.settings;
  const drag = (el, onMove, reset) => {
    el.addEventListener('mousedown', e => {
      e.preventDefault(); el.classList.add('drag'); document.body.classList.add('resizing');
      const mv = ev => { onMove(ev); applyLayout(); };
      const up = () => { el.classList.remove('drag'); document.body.classList.remove('resizing'); window.removeEventListener('mousemove', mv); window.removeEventListener('mouseup', up); saveSettings(); requestUpdate(false); };
      window.addEventListener('mousemove', mv); window.addEventListener('mouseup', up);
    });
    el.addEventListener('dblclick', () => { reset(); applyLayout(); saveSettings(); requestUpdate(false); });
  };
  drag($('#splitMS'), ev => { const vp = $('#viewport').getBoundingClientRect(); S.mapW = Math.round(clamp(ev.clientX - vp.left, 200, vp.width - 260)); }, () => { S.mapW = 0; });
  drag($('#splitSide'), ev => { const r = main.getBoundingClientRect(); S.sideW = Math.round(clamp(r.right - ev.clientX, 220, Math.min(720, r.width - 320))); }, () => { S.sideW = 0; });
}
function wireMini() {
  const host = $('#miniMap'); let down = false;
  const go = e => { const svg = host.querySelector('svg'); if (!svg || !App.cur) return; const r = svg.getBoundingClientRect(); const p = (e.clientX - r.left - MINI.M) / (r.width - 2 * MINI.M) * MINI.n; scrollSeqToPos(clamp(Math.round(p), 0, MINI.n)); };
  host.addEventListener('mousedown', e => {
    const fe = e.target.closest('[data-fid]'), doc = App.cur; if (!doc) return;
    if (fe) { const f = featureById(doc, fe.dataset.fid); if (f) { selectFeature(doc, f); } return; }
    down = true; go(e); e.preventDefault();
  });
  window.addEventListener('mousemove', e => { if (down) go(e); });
  window.addEventListener('mouseup', () => { down = false; });
  $('#seqScroll').addEventListener('scroll', () => { if (App.settings.view === 'seq' && !wireMini.raf) wireMini.raf = requestAnimationFrame(() => { wireMini.raf = 0; updateMiniWindow(App.cur); }); });
}

function init() {
  applyTheme();
  document.addEventListener('click', e => { const b = e.target.closest('[data-act]'); if (b && actions[b.dataset.act]) actions[b.dataset.act](e); });
  $('#sideTabs').addEventListener('click', e => { const b = e.target.closest('button'); if (!b) return; App.settings.sideTab = b.dataset.tab; saveSettings(); if (App.cur) renderSide(App.cur); });
  $('#enzMode').addEventListener('change', e => { App.settings.enzMode = e.target.value; saveSettings(); requestUpdate(true); });
  const find = $('#find');
  find.addEventListener('input', debounce(() => { if (App.cur) { UI.find = { q: find.value, res: findAll(App.cur, find.value), idx: -1 }; requestUpdate(false); } }, 150));
  find.addEventListener('keydown', e => { if (e.key === 'Enter') { e.preventDefault(); runFind(e.shiftKey ? -1 : 1); } else if (e.key === 'Escape') { find.blur(); } });
  applyLayout(); wireMap(); wireSeq(); wireKeys(); wireHover(); wirePos(); wireSplit(); wireMini();
  new ResizeObserver(debounce(() => { if (App.cur) requestUpdate(false); }, 60)).observe($('#viewport'));
  { const ro = new ResizeObserver(debounce(() => { if (App.cur) requestUpdate(false); }, 60)); ro.observe($('#mapPane')); ro.observe($('#seqPane')); }
  // drag & drop
  window.addEventListener('dragover', e => { e.preventDefault(); document.body.classList.add('dragging'); });
  window.addEventListener('dragleave', e => { if (!e.relatedTarget) document.body.classList.remove('dragging'); });
  window.addEventListener('drop', e => { e.preventDefault(); document.body.classList.remove('dragging'); if (e.dataTransfer.files.length) openFiles(Array.from(e.dataTransfer.files)); });
  window.addEventListener('beforeunload', e => { if (App.docs.some(d => d.dirty)) { e.preventDefault(); e.returnValue = ''; } });
  if (window.plasmidNative && window.plasmidNative.onOpenFile) window.plasmidNative.onOpenFile(({ name, data, path }) => { try { const bin = atob(data), u8 = Uint8Array.from(bin, c => c.charCodeAt(0)); const d = addDoc(parseAny(u8.buffer, name), { detect: true }); d.path = path; } catch (err) { toast(err.message); } });
  requestUpdate(true);
}
document.addEventListener('DOMContentLoaded', init);
