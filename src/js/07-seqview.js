'use strict';
/* ---------- sequence editor view (DOM based, monospace grid) ---------- */

const SV = { lineH: 19, laneH: 21, enzH: 15, gutter: 84, cw: 0, bpr: 60, rows: [], tops: [], key: '', touched: [] };

function measureChar() {
  const s = el('span', { class: 'seqfont', style: 'position:absolute;visibility:hidden;white-space:pre' }); s.textContent = 'ACGTACGTAC'.repeat(10);
  document.body.append(s); const w = s.getBoundingClientRect().width / 100; s.remove(); return w || 8.4;
}

function seqLayoutKey(doc, bpr) {
  const s = App.settings;
  return [doc.id, doc.rev, bpr, s.showFeatures, s.showOrfs, s.orfMin, s.showTrans, s.enzMode, s.enzShow.join(','), enzVersion].join('|');
}

function renderSeq(doc, force) {
  const scroll = $('#seqScroll'), host = $('#seqRows');
  if (!doc) { host.innerHTML = ''; return; }
  if (!SV.cw) SV.cw = measureChar();
  const cw = SV.cw, n = doc.seq.length;
  const avail = scroll.clientWidth - SV.gutter - 28;
  const bpr = clamp(Math.floor(avail / cw / 10) * 10, 10, 400);
  const key = seqLayoutKey(doc, bpr);
  if (!force && key === SV.key) { updateSeqSelection(doc); return; }
  const prevScroll = scroll.scrollTop;
  SV.key = key; SV.bpr = bpr;
  const rows = Math.max(1, Math.ceil(n / bpr));

  const showTrans = App.settings.showTrans;
  const feats = visibleFeatures(doc);
  if (showTrans) for (const f of visibleFeatures(doc)) { if (f.type === 'CDS' && !f.orf) { const x = cdsExtension(doc, f); if (x) feats.push(x); } }
  const { lane } = assignLanes(feats.map(f => ({ id: f.id, ivs: f.locs })), 0);
  const isTrans = f => showTrans && (f.type === 'CDS' || f.orf);
  const rowFeat = Array.from({ length: rows }, () => []);
  for (const f of feats) {
    const k = lane.get(f.id), headSeg = f.strand === -1 ? 0 : f.locs.length - 1;
    f.locs.forEach(([a, b], i) => {
      for (let r = Math.floor(a / bpr); r <= Math.floor((b - 1) / bpr); r++) {
        const rs = r * bpr, x = Math.max(a, rs), y = Math.min(b, rs + bpr);
        const head = f.strand !== 0 && i === headSeg;
        rowFeat[r].push({ f, lane: k, x, y, arrR: head && f.strand === 1 && y === b, arrL: head && f.strand === -1 && x === a });
      }
    });
  }
  // enzymes per row
  const enz = enzymeView(doc), rowEnz = Array.from({ length: rows }, () => []), rowSite = Array.from({ length: rows }, () => []), rowBot = Array.from({ length: rows }, () => []);
  const rowOf = p => Math.min(Math.floor(p / bpr), rows - 1);
  for (const e of enz) {
    const r = rowOf(e.top);
    rowEnz[r].push(e);
    rowBot[rowOf(e.bot)].push(e);
    const pieces = e.e > n ? [[e.s, n], [0, e.e - n]] : [[e.s, e.e]];
    for (const [a, b] of pieces) for (let rr = Math.floor(a / bpr); rr <= Math.floor((b - 1) / bpr) && rr < rows; rr++) rowSite[rr].push({ e, x: Math.max(a, rr * bpr), y: Math.min(b, (rr + 1) * bpr) });
  }

  const html = [], tops = []; let y = 0;
  for (let r = 0; r < rows; r++) {
    const rs = r * bpr, re = Math.min(n, rs + bpr), len = re - rs;
    // enzyme label lanes
    const els = rowEnz[r].map(e => ({ e, col: e.top - rs, w: (e.name.length * 6.9 + 8) })).sort((a, b) => a.col - b.col);
    const ends = [];
    for (const it of els) { let k = 0; while (k < ends.length && ends[k] > it.col * cw - 2) k++; ends[k] = it.col * cw + it.w; it.k = k; }
    const eLanes = ends.length, eH = eLanes ? eLanes * SV.enzH + 4 : 6;
    // lane heights: plain bar 21px; translated CDS 31px (residue numbers above); in-frame continuation 44px (+ dashed arrow line)
    const fLanes = rowFeat[r].length ? Math.max(...rowFeat[r].map(o => o.lane)) + 1 : 0;
    const laneH = Array(fLanes).fill(SV.laneH);
    for (const o of rowFeat[r]) laneH[o.lane] = Math.max(laneH[o.lane], o.f.ext ? 44 : isTrans(o.f) ? 31 : SV.laneH);
    const laneTop = []; let acc = 0; for (let k = 0; k < fLanes; k++) { laneTop.push(acc); acc += laneH[k]; }
    const H = eH + 2 * SV.lineH + (fLanes ? acc + 8 : 6) + 8;
    tops.push(y); y += H;

    const b = [];
    b.push(`<div class="row" data-r="${r}" style="height:${H}px">`);
    b.push(`<div class="gut" style="top:${eH}px">${(rs + 1).toLocaleString('en-US')}</div>`);
    b.push(`<div class="rseq" style="left:${SV.gutter}px;width:${Math.max(len, 1) * cw}px">`);
    for (const s of rowSite[r]) b.push(`<div class="eshade" data-ek="${esc(s.e.name)}:${s.e.top}" style="left:${f2((s.x - rs) * cw)}px;width:${f2((s.y - s.x) * cw)}px;top:${eH}px;height:${2 * SV.lineH}px"></div>`);
    for (const it of els) {
      const x = it.col * cw, top = it.k * SV.enzH;
      b.push(`<div class="elbl" data-enz="${esc(it.e.name)}" data-ek="${esc(it.e.name)}:${it.e.top}" data-s="${it.e.s}" data-e="${it.e.e}" style="left:${f2(x)}px;top:${top}px">${esc(it.e.name)}</div>`);
      b.push(`<div class="estem" style="left:${f2(x)}px;top:${top + SV.enzH - 1}px;height:${eH - (top + SV.enzH - 1)}px"></div>`);
    }
    b.push(`<div class="ln top" style="top:${eH}px">${doc.seq.slice(rs, re)}</div>`);
    b.push(`<div class="ln bot" style="top:${eH + SV.lineH}px">${complement(doc.seq.slice(rs, re))}</div>`);
    // cut marks: red zig-zag joining the top-strand and bottom-strand cut positions
    const Hh = SV.lineH, zig = (xt, xb, kind, e) => {
      const x0 = Math.min(xt, xb) - 3, w = Math.abs(xt - xb) + 6, a = xt - x0, c = xb - x0;
      const pts = kind === 'full' ? `${f2(a)},0 ${f2(a)},${Hh} ${f2(c)},${Hh} ${f2(c)},${2 * Hh}` : kind === 'top' ? `${f2(a)},0 ${f2(a)},${Hh}` : `${f2(c)},${Hh} ${f2(c)},${2 * Hh}`;
      return `<svg class="ezig" data-ek="${esc(e.name)}:${e.top}" style="left:${f2(x0)}px;top:${eH}px" width="${f2(w)}" height="${2 * Hh}"><polyline points="${pts}"/></svg>`;
    };
    for (const e of rowEnz[r]) {
      const xt = (e.top - rs) * cw;
      if (rowOf(e.bot) === r) b.push(zig(xt, (e.bot - rs) * cw, 'full', e)); else b.push(zig(xt, xt, 'top', e));
    }
    for (const e of rowBot[r]) if (rowOf(e.top) !== r) { const xb = Math.min(e.bot - rs, len) * cw; b.push(zig(xb, xb, 'bot', e)); }
    b.push(`<div class="selLayer" style="top:${eH}px;height:${2 * SV.lineH}px"></div>`);
    const fy0 = eH + 2 * SV.lineH + 6;
    for (const o of rowFeat[r]) {
      const f = o.f, left = (o.x - rs) * cw, w = (o.y - o.x) * cw, top = fy0 + laneTop[o.lane], tr = isTrans(f);
      const cls = 'fbar' + (o.arrR && !f.ext ? ' arR' : '') + (o.arrL && !f.ext ? ' arL' : '') + (f.orf ? ' orf' : '') + (f.ext ? ' ext' : '') + (doc.selFid === (f.ext ? f.parent.id : f.id) ? ' on' : '');
      const fid = f.ext ? f.parent.id : f.id;
      let inner, nums = '';
      if (tr) {
        const cod = cdsCodons(doc, f), off = f.aaOffset || 0, parts = [];
        cod.forEach((c, i) => {
          const m = c.pos[1]; if (m < o.x || m >= o.y) return;
          parts.push(`<i style="left:${f2((m - o.x) * cw)}px;width:${f2(cw)}px">${c.aa}</i>`);
          const num = off + i + 1; if (num === 1 || num % 10 === 0) nums += `<b style="left:${f2((m - o.x) * cw - 14)}px;width:${f2(cw + 28)}px">${num}</b>`;
        });
        inner = parts.join('');
        b.push(`<div class="fnum" style="left:${f2(left)}px;top:${top}px;width:${f2(w)}px">${nums}</div>`);
      } else inner = `<span>${esc(f.name)}</span>`;
      const by = top + (tr ? 10 : 0), bg = f.ext ? tint(f.color, 0.22) : f.color, fg = f.ext ? 'var(--text)' : textOn(f.color);
      b.push(`<div class="${cls}" data-fid="${fid}" style="left:${f2(left)}px;top:${by}px;width:${f2(w)}px;background:${bg};color:${fg}${f.ext ? ';outline-color:' + f.color : ''}">${inner}</div>`);
      if (f.ext) {
        const label = f.name + (f.stopped ? '' : ' (no stop codon found)');
        b.push(`<div class="extline${o.arrR ? ' hr' : ''}${o.arrL ? ' hl' : ''}" style="left:${f2(left)}px;top:${by + 21}px;width:${f2(w)}px;--ecol:${f.color}">${w > textWidth(label, '600 9.5px system-ui') + 34 ? `<span class="extlbl">${esc(label)}</span>` : ''}</div>`);
      }
    }
    b.push('</div></div>');
    html.push(b.join(''));
  }
  host.innerHTML = html.join('');
  SV.rows = $$('.row', host); SV.tops = tops; SV.heights = null; SV.touched = [];
  host.style.setProperty('--cw', cw + 'px');
  scroll.scrollTop = prevScroll;
  updateSeqSelection(doc);
}

function tint(hex, a) {
  const m = /^#?([0-9a-f]{6})$/i.exec(hex || ''); if (!m) return 'rgba(127,127,127,' + a + ')';
  const v = parseInt(m[1], 16); return `rgba(${v >> 16},${(v >> 8) & 255},${v & 255},${a})`;
}
function seqRowOfPos(p) { return Math.min(Math.floor(p / SV.bpr), Math.max(0, SV.rows.length - 1)); }

function updateSeqSelection(doc) {
  if (!doc || !SV.rows.length) return;
  for (const r of SV.touched) { const l = r.querySelector('.selLayer'); if (l) l.innerHTML = ''; }
  SV.touched = [];
  const n = doc.seq.length, bpr = SV.bpr, cw = SV.cw, sel = selRange(doc), LH = SV.lineH;
  const box = (a0, b0, cls, strand) => {   // paints [a0,b0) over the whole row pair, or only the matched strand's line
    for (let r = Math.floor(a0 / bpr); r <= Math.floor((b0 - 1) / bpr) && r < SV.rows.length; r++) {
      const rs = r * bpr, a = Math.max(a0, rs), b = Math.min(b0, rs + bpr);
      const st = strand === undefined ? '' : `top:${strand === 1 ? 0 : LH}px;bottom:auto;height:${LH}px;`;
      SV.rows[r].querySelector('.selLayer').insertAdjacentHTML('beforeend', `<div class="${cls}" style="${st}left:${f2((a - rs) * cw)}px;width:${f2((b - a) * cw)}px"></div>`);
      SV.touched.push(SV.rows[r]);
    }
  };
  const F = (typeof UI !== 'undefined' && UI.find && UI.find.q) ? UI.find : null;
  const curHit = F && F.idx >= 0 ? F.res[F.idx] : null;
  const curIsSel = !!(curHit && sel && sel[0] === curHit.s && sel[1] === curHit.e);
  if (F && F.res.length <= 3000) for (const m of F.res) if (!(curIsSel && m === curHit)) box(m.s, m.e, 'fm', m.strand);
  if (curIsSel) box(curHit.s, curHit.e, 'fm cur', curHit.strand);
  else if (sel) box(sel[0], sel[1], 'selbox');
  if (!sel) {
    const r = seqRowOfPos(doc.caret), col = doc.caret - r * bpr;
    SV.rows[r].querySelector('.selLayer').insertAdjacentHTML('beforeend', `<div class="caret" style="left:${f2(col * cw)}px"></div>`);
    SV.touched.push(SV.rows[r]);
  }
  if (doc._reveal) {
    doc._reveal = false;
    const scroll = $('#seqScroll'), p = sel ? (doc.caret === sel[1] ? sel[1] - 1 : sel[0]) : doc.caret;
    const r = seqRowOfPos(clamp(p, 0, Math.max(0, n - 1))), top = SV.tops[r], h = SV.rows[r].offsetHeight;
    if (top < scroll.scrollTop + 4) scroll.scrollTop = Math.max(0, top - 30);
    else if (top + h > scroll.scrollTop + scroll.clientHeight - 4) scroll.scrollTop = top + h - scroll.clientHeight + 30;
  }
}

/* mouse → caret position */
function seqPosFromEvent(e) {
  const scroll = $('#seqScroll'), r0 = scroll.getBoundingClientRect();
  const y = e.clientY - r0.top + scroll.scrollTop;
  let lo = 0, hi = SV.tops.length - 1;
  while (lo < hi) { const mid = (lo + hi + 1) >> 1; if (SV.tops[mid] <= y) lo = mid; else hi = mid - 1; }
  const row = SV.rows[lo]; if (!row) return 0;
  const rs = row.querySelector('.rseq').getBoundingClientRect();
  const len = Math.min(SV.bpr, App.cur.seq.length - lo * SV.bpr);
  const col = clamp(Math.round((e.clientX - rs.left) / SV.cw), 0, len);
  return lo * SV.bpr + col;
}


/* ---------- mini linear map under the sequence (sequence-only view) ---------- */
const MINI = { M: 26, W: 0, n: 0 };
function renderMini(doc) {
  const host = $('#miniMap');
  if (!doc || App.settings.view !== 'seq') { host.innerHTML = ''; return; }
  const n = doc.seq.length, W = Math.max(300, host.clientWidth || 600), M = MINI.M, X = p => M + p / n * (W - 2 * M);
  MINI.W = W; MINI.n = n;
  const feats = visibleFeatures(doc), out = [];
  const { lane } = assignLanes(feats.map(f => ({ id: f.id, ivs: f.locs.map(([a, b]) => [X(a), X(b)]) })), 2);
  const baseY = 26, FY = baseY + 12, FH = 9, LHm = 11, H = FY + 3 * LHm + 8;
  const step = niceStep(n, Math.max(4, Math.round(W / 100)));
  out.push(`<line class="mp-ring" x1="${M}" y1="${baseY}" x2="${W - M}" y2="${baseY}" stroke-width="3"/>`);
  for (let p = 0; p <= n; p += step) {
    out.push(`<line class="mp-tick" x1="${f2(X(p))}" y1="${baseY - 8}" x2="${f2(X(p))}" y2="${baseY - 3}"/><text class="mini-t" x="${f2(X(p))}" y="${baseY - 11}" text-anchor="middle">${Math.round(p).toLocaleString('en-US')}</text>`);
  }
  for (const f of feats) {
    const k = Math.min(lane.get(f.id), 2), y = FY + k * LHm;
    for (const [a, b] of f.locs) out.push(`<rect class="mp-feat${f.orf ? ' orf' : ''}${doc.selFid === f.id ? ' on' : ''}" data-fid="${f.id}" x="${f2(X(a))}" y="${y}" width="${f2(Math.max(2, X(b) - X(a)))}" height="${FH}" rx="2" fill="${f.color}"/>`);
  }
  for (const e of enzymeView(doc)) out.push(`<line class="mp-cut" x1="${f2(X(e.top))}" y1="${baseY - 2}" x2="${f2(X(e.top))}" y2="${baseY + 8}"/>`);
  const sel = selRange(doc);
  if (sel) out.push(`<rect class="mp-sel" x="${f2(X(sel[0]))}" y="2" width="${f2(Math.max(2, X(sel[1]) - X(sel[0])))}" height="${H - 4}"/>`);
  else out.push(`<line class="mp-caret" x1="${f2(X(doc.caret))}" y1="2" x2="${f2(X(doc.caret))}" y2="${H - 2}"/>`);
  out.push(`<rect id="miniWin" class="mini-win" x="0" y="2" width="0" height="${H - 4}" rx="3"/>`);
  host.innerHTML = `<svg viewBox="0 0 ${W} ${H}" width="${W}" height="${H}" data-n="${n}">${out.join('')}</svg>`;
  updateMiniWindow(doc);
}
function visibleRange() {
  const sc = $('#seqScroll'), st = sc.scrollTop, vb = st + sc.clientHeight, T = SV.tops;
  if (!T.length) return [0, 0];
  const idx = y => { let lo = 0, hi = T.length - 1; while (lo < hi) { const m = (lo + hi + 1) >> 1; if (T[m] <= y) lo = m; else hi = m - 1; } return lo; };
  return [idx(st) * SV.bpr, Math.min(App.cur.seq.length, (idx(Math.max(st, vb - 1)) + 1) * SV.bpr)];
}
function updateMiniWindow(doc) {
  const w = $('#miniWin'); if (!w || !doc || !SV.tops.length) return;
  const [a, b] = visibleRange(), X = p => MINI.M + p / MINI.n * (MINI.W - 2 * MINI.M);
  w.setAttribute('x', f2(X(a))); w.setAttribute('width', f2(Math.max(4, X(b) - X(a))));
}
function scrollSeqToPos(p) {
  const sc = $('#seqScroll'), r = seqRowOfPos(clamp(p, 0, Math.max(0, App.cur.seq.length - 1)));
  sc.scrollTop = Math.max(0, SV.tops[r] - sc.clientHeight / 2 + (SV.rows[r] ? SV.rows[r].offsetHeight / 2 : 0));
}
