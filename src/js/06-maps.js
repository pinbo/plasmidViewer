'use strict';
/* ---------- circular & linear map renderers (SVG) ---------- */

const CW = 1000, CH = 880, CCX = 500, CCY = 440;
const polar = (r, a) => [CCX + r * Math.sin(a), CCY - r * Math.cos(a)];
const f2 = v => v.toFixed(2);

function arcSeg(r1, r2, a1, a2, head, headPx) {
  // head: 1 → tip at a2, -1 → tip at a1, 0 → none. angles clockwise from top.
  const rm = (r1 + r2) / 2, P = (r, a) => polar(r, a).map(f2).join(' ');
  const large = (x, y) => (y - x > Math.PI ? 1 : 0);
  if (!head) {
    return `M${P(r2, a1)} A${r2} ${r2} 0 ${large(a1, a2)} 1 ${P(r2, a2)} L${P(r1, a2)} A${r1} ${r1} 0 ${large(a1, a2)} 0 ${P(r1, a1)}Z`;
  }
  let da = Math.min(headPx / rm, a2 - a1);
  if (head === 1) {
    const ab = a2 - da, flare = 3;
    if (ab <= a1 + 1e-6) return `M${P(r2 + flare, a1)} L${P(rm, a2)} L${P(r1 - flare, a1)}Z`;
    return `M${P(r2, a1)} A${r2} ${r2} 0 ${large(a1, ab)} 1 ${P(r2, ab)} L${P(r2 + flare, ab)} L${P(rm, a2)} L${P(r1 - flare, ab)} L${P(r1, ab)} A${r1} ${r1} 0 ${large(a1, ab)} 0 ${P(r1, a1)}Z`;
  }
  const ab = a1 + da, flare = 3;
  if (ab >= a2 - 1e-6) return `M${P(r2 + flare, a2)} L${P(rm, a1)} L${P(r1 - flare, a2)}Z`;
  return `M${P(rm, a1)} L${P(r2 + flare, ab)} L${P(r2, ab)} A${r2} ${r2} 0 ${large(ab, a2)} 1 ${P(r2, a2)} L${P(r1, a2)} A${r1} ${r1} 0 ${large(ab, a2)} 0 ${P(r1, ab)} L${P(r1 - flare, ab)}Z`;
}

/* labels of enzymes that cut at the same position share one label: "BanII, SacI (406)" */
function groupEnzymes(list) {
  const g = new Map();
  for (const e of list) { const k = e.top; if (!g.has(k)) g.set(k, { top: e.top, items: [] }); g.get(k).items.push(e); }
  return [...g.values()].map(x => { x.items.sort((p, q) => p.name.localeCompare(q.name)); x.names = x.items.map(i => i.name); return x; });
}
const enzGroupSvg = (g, attrs = '') => `<text ${attrs}>${g.items.map((e, i) => `<tspan data-enz="${esc(e.name)}" data-s="${e.s}" data-e="${e.e}">${esc(e.name)}${i < g.items.length - 1 ? ', ' : ` (${g.top})`}</tspan>`).join('')}</text>`;   // separators and the position belong to a name's tspan so the whole label is hoverable

/* Least-squares placement of labels along a vertical axis: keeps the order, keeps each label as near as possible to its wish (d) with room for the
   labels' heights h between neighbours (pool-adjacent-violators on d_i - offset_i), then fits the whole stack into [lo, hi]. */
function placeStack(d, h, lo, hi) {
  const m = d.length; if (!m) return [];
  const off = [0]; for (let i = 1; i < m; i++) off.push(off[i - 1] + (h[i - 1] + h[i]) / 2);
  const k = off[m - 1] > hi - lo ? (hi - lo) / off[m - 1] : 1; for (let i = 0; i < m; i++) off[i] *= k;
  const blocks = [];
  d.forEach((v, i) => { blocks.push({ s: v - off[i], w: 1, n: 1 }); while (blocks.length > 1 && blocks[blocks.length - 2].s / blocks[blocks.length - 2].w > blocks[blocks.length - 1].s / blocks[blocks.length - 1].w) { const b = blocks.pop(), a = blocks[blocks.length - 1]; a.s += b.s; a.w += b.w; a.n += b.n; } });
  const y = []; for (const b of blocks) for (let q = 0; q < b.n; q++) y.push(b.s / b.w);
  return y.map((v, i) => Math.min(Math.max(v + off[i], lo + off[i]), hi - (off[m - 1] - off[i])));
}
/* Split ring labels (sorted by angle) into a right and a left column. The right column is a contiguous arc starting at the top; its ends are moved
   (at most ~55° away from the vertical axis) so both columns carry a similar number of labels – contiguous arcs keep the leader lines from crossing. */
function splitRingLabels(L) {
  const N = L.length; let start = 0, cnt = L.filter(l => l.theta <= Math.PI).length;
  const dTop = t => Math.min(t, 2 * Math.PI - t), dBot = t => Math.abs(t - Math.PI), lim = 55 * Math.PI / 180, thr = Math.max(2, N * 0.2), cap = 50;   // only rebalance when one column would be too full to keep the names near their sites
  for (let guard = 0; guard < N; guard++) {
    const diff = cnt - (N - cnt); if (Math.abs(diff) <= thr || Math.max(cnt, N - cnt) <= cap) break;
    if (diff > 0) {   // right column too full: give away its first (near the top) or last (near the bottom) label
      const a = L[start % N], z = L[(start + cnt - 1) % N], da = dTop(a.theta), dz = dBot(z.theta);
      if (Math.min(da, dz) > lim) break; if (da <= dz) { start = (start + 1) % N; cnt--; } else cnt--;
    } else {          // left column too full: take the label before the right arc's start (near the top) or the one after its end (near the bottom)
      const a = L[(start - 1 + N) % N], z = L[(start + cnt) % N], da = dTop(a.theta), dz = dBot(z.theta);
      if (Math.min(da, dz) > lim) break; if (da <= dz) { start = (start - 1 + N) % N; cnt++; } else cnt++;
    }
  }
  const right = [], left = [];
  for (let i = 0; i < N; i++) { const l = L[(start + i) % N]; (i < cnt ? right : left).push(l); }
  return { right, left: left.reverse() };   // both lists run from top to bottom
}

function renderCircular(doc, host) {
  const n = doc.seq.length, feats = visibleFeatures(doc);
  const R = 160, BH = 16.5, LH = 19.5;
  const { lane, count } = assignLanes(feats.map(f => ({ id: f.id, ivs: f.locs.concat(hangSegs(f, n)) })), n * 0.002);
  const lanes = Math.max(1, count), rTop = R + (lanes - 1) * LH + BH / 2 + 4;
  const ang = p => 2 * Math.PI * p / n;
  const sel = selRange(doc), segs = selSegs(doc);
  const out = [];
  out.push(`<circle class="mp-ring" cx="${CCX}" cy="${CCY}" r="${R}" fill="none" stroke-width="3"/>`);

  // selection band
  for (const [x, y] of segs) out.push(`<path class="mp-sel" d="${arcSeg(R - 34, rTop + 6, ang(x), Math.min(ang(y), 2 * Math.PI - 0.002), 0, 0)}"/>`);
  // ticks
  const step = niceStep(n, 12), minor = step / 5;
  for (let p = 0; p < n; p += minor) {
    const major = Math.abs(p % step) < 1e-9 || Math.abs((p % step) - step) < 1e-9;
    const a = ang(p), [x1, y1] = polar(R - 14, a), [x2, y2] = polar(R - (major ? 22 : 18), a);
    out.push(`<line class="mp-tick" x1="${f2(x1)}" y1="${f2(y1)}" x2="${f2(x2)}" y2="${f2(y2)}"/>`);
    if (major && p > 0 && (n - p) > step * 0.35) { const [tx, ty] = polar(R - 36, a); out.push(`<text class="mp-ticklbl" x="${f2(tx)}" y="${f2(ty)}" text-anchor="middle" dominant-baseline="central">${Math.round(p).toLocaleString('en-US')}</text>`); }
  }
  // features (leader lines are inserted below them afterwards, so they never run over a shape)
  const labels = [], underAt = out.length;
  for (const f of feats) {
    const k = lane.get(f.id), rc = R + k * LH, r1 = rc - BH / 2, r2 = rc + BH / 2;
    const headSeg = f.strand === -1 ? 0 : f.locs.length - 1;
    f.locs.forEach(([a, b], i) => {
      const head = f.strand === 0 || i !== headSeg ? 0 : f.strand;
      const a1 = ang(a), a2 = Math.min(ang(b), 2 * Math.PI - 0.002);
      out.push(`<path class="mp-feat${doc.selFid === f.id ? ' on' : ''}${f.orf ? ' orf' : ''}" data-fid="${f.id}" d="${arcSeg(r1, r2, a1, a2, head, 14)}" fill="${f.color}"></path>`);
    });
    for (const [a, b] of hangSegs(f, n)) {
      const a1 = ang(a), a2 = ang(b), cnt = Math.max(2, Math.round((a2 - a1) * rc / 7)), pts = [];
      for (let i = 0; i <= cnt; i++) pts.push(polar(i % 2 ? r1 + 2 : r2 - 2, a1 + (a2 - a1) * i / cnt).map(f2).join(','));
      out.push(`<polyline class="mp-hang${doc.selFid === f.id ? ' on' : ''}" data-fid="${f.id}" points="${pts.join(' ')}" stroke="${f.color}"/>`);
    }
    // label anchor: middle of the longest segment
    const seg = f.locs.reduce((m, l) => (l[1] - l[0] > m[1] - m[0] ? l : m));
    labels.push({ kind: 'f', f, theta: ang((seg[0] + seg[1]) / 2), rr: r2 + 2, text: f.name, color: f.color });
  }
  // enzymes (sites cutting at the same position share a label and a line; the location tick is short, the leader continues from its end)
  const isSel = e => !!sel && sel[0] === e.s && sel[1] === e.e, cutR0 = R - BH * 0.55, cutR1 = R + BH * 0.55;   // tick = 10% longer than the feature shapes are wide
  for (const g of groupEnzymes(enzymeView(doc))) {
    const a = ang(g.top), [x1, y1] = polar(cutR0, a), [x2, y2] = polar(cutR1, a);
    out.push(`<line class="mp-cut${g.items.some(isSel) ? ' on' : ''}" data-enz="${esc(g.names.join('|'))}" x1="${f2(x1)}" y1="${f2(y1)}" x2="${f2(x2)}" y2="${f2(y2)}"/>`);
    labels.push({ kind: 'e', g, theta: a, rr: cutR1 });
  }
  // label layout: names stay around the ring in two columns that follow the circle (right / left of the vertical axis); each column is filled
  // evenly, every name as near to its position as the stacking allows, and the order matches the order around the ring so leaders do not cross
  const Rl = rTop + 22, Re = rTop + 8, minY = 22, maxY = CH - 22, under = [], text = [];
  labels.sort((p, q) => p.theta - q.theta);
  const cols = splitRingLabels(labels);
  const segMin = (ax, ay, bx, by) => { const dx = bx - ax, dy = by - ay, t = Math.max(0, Math.min(1, -((ax - CCX) * dx + (ay - CCY) * dy) / (dx * dx + dy * dy || 1))); return Math.hypot(ax + t * dx - CCX, ay + t * dy - CCY); };
  for (const [s, list] of [[1, cols.right], [-1, cols.left]]) {
    const ys = placeStack(list.map(l => CCY - Rl * Math.cos(l.theta)), list.map(l => (l.kind === 'e' ? 14 : 18)), minY, maxY);
    list.forEach((l, i) => {
      const y = ys[i], dy = y - CCY, x = CCX + s * (Math.sqrt(Math.max(Rl * Rl - dy * dy, (0.35 * Rl) ** 2)) + 10);
      const [ax, ay] = polar(l.rr, l.theta), lx = x - s * 3;
      // straight from the start of the leader to the near end of the name (its start on the right, its end on the left); if that would dip back
      // towards the ring, leave the shape radially first
      let pts = `${f2(ax)},${f2(ay)} ${f2(lx)},${f2(y)}`;
      if (segMin(ax, ay, lx, y) < l.rr - 0.5) { const [ex, ey] = polar(Math.max(Re, l.rr), l.theta); pts = `${f2(ax)},${f2(ay)} ${f2(ex)},${f2(ey)} ${f2(lx)},${f2(y)}`; }
      const tattr = `x="${f2(x)}" y="${f2(y)}" text-anchor="${s === 1 ? 'start' : 'end'}" dominant-baseline="central"`;
      if (l.kind === 'e') {
        under.push(`<polyline class="mp-lead${l.g.items.some(isSel) ? ' on' : ''}" data-enz="${esc(l.g.names.join('|'))}" points="${pts}"/>`);
        text.push(enzGroupSvg(l.g, `class="mp-lbl enz" ${tattr}`));
      } else {
        under.push(`<polyline class="mp-lead${doc.selFid === l.f.id ? ' on' : ''}" data-fid="${l.f.id}" points="${pts}"/>`);
        text.push(`<text class="mp-lbl" data-fid="${l.f.id}" ${tattr}>${esc(l.text)}</text>`);
      }
    });
  }
  out.splice(underAt, 0, ...under); out.push(...text);
  // caret
  if (!segs.length) {
    const a = ang(doc.caret), [x1, y1] = polar(R - 30, a), [x2, y2] = polar(rTop + 8, a);
    out.push(`<line class="mp-caret" x1="${f2(x1)}" y1="${f2(y1)}" x2="${f2(x2)}" y2="${f2(y2)}"/>`);
  }
  // centre text
  const nm = doc.name.length > 26 ? doc.name.slice(0, 25) + '…' : doc.name;
  out.push(`<text class="mp-title" x="${CCX}" y="${CCY - 16}" text-anchor="middle">${esc(nm)}</text>`);
  out.push(`<text class="mp-sub" x="${CCX}" y="${CCY + 12}" text-anchor="middle">${n.toLocaleString('en-US')} bp</text>`);
  if (segs.length) out.push(`<text class="mp-sub small" x="${CCX}" y="${CCY + 36}" text-anchor="middle">${segLabel(segs)} · ${segs.reduce((t, s) => t + s[1] - s[0], 0).toLocaleString('en-US')} bp selected</text>`);
  host.innerHTML = `<svg class="map circ" viewBox="0 0 ${CW} ${CH}" data-n="${n}" preserveAspectRatio="xMidYMid meet">${out.join('')}</svg>`;
}

/* converts a mouse event on the circular map to a sequence position + whether it is on the ring area */
function circularPos(svg, e, n) {
  const pt = svg.createSVGPoint(); pt.x = e.clientX; pt.y = e.clientY;
  const p = pt.matrixTransform(svg.getScreenCTM().inverse());
  const dx = p.x - CCX, dy = p.y - CCY;
  let a = Math.atan2(dx, -dy); if (a < 0) a += 2 * Math.PI;
  return { pos: Math.round(a / (2 * Math.PI) * n) % n, base: Math.min(n - 1, Math.floor(a / (2 * Math.PI) * n)), r: Math.hypot(dx, dy) };
}

/* ---------------- linear map ---------------- */
const LIN_M = 40;
function renderLinear(doc, host) {
  const n = doc.seq.length, feats = visibleFeatures(doc);
  const wrap = host.clientWidth || 900, zoom = App.settings.linZoom;
  const W = Math.max(500, Math.round((wrap - 2) * zoom));
  const X = p => LIN_M + p / n * (W - 2 * LIN_M);
  const font = '12px system-ui, sans-serif';
  const out = [];
  const sel = selRange(doc);

  // enzyme labels (above)
  const enz = groupEnzymes(enzymeView(doc)).map(g => ({ ...g, text: `${g.names.join(', ')} (${g.top})`, x: X(g.top) }));
  enz.forEach(e => { e.w = textWidth(e.text, font) + 8; });
  enz.sort((a, b) => a.x - b.x);
  const rowEnd = [];
  for (const e of enz) {
    let k = 0; while (k < rowEnd.length && rowEnd[k] > e.x - 2) k++;
    rowEnd[k] = e.x + e.w; e.row = k;
  }
  const enzRows = rowEnd.length, EH = 15;
  const baseY = 30 + enzRows * EH + 14;

  // features lanes
  const items = feats.map(f => {
    const nameW = textWidth(f.name, '700 16px system-ui, sans-serif') + 8;
    const [a, b] = featBounds(f); const bw = (X(b) - X(a));
    const inside = nameW <= bw - 14;
    const ivs = f.locs.concat(hangSegs(f, n)).map(([s, e]) => [X(s), X(e)]);
    if (!inside) { const last = f.strand === -1 ? ivs[0] : ivs[ivs.length - 1]; ivs.push(f.strand === -1 ? [last[0] - nameW, last[0]] : [last[1], last[1] + nameW]); }
    return { id: f.id, ivs, inside, nameW };
  });
  const { lane, count } = assignLanes(items.map(i => ({ id: i.id, ivs: i.ivs })), 3);
  const info = new Map(items.map(i => [i.id, i]));
  const FH = 16.5, LH = 22.5, featY0 = baseY + 46;
  const H = featY0 + Math.max(1, count) * LH + 24;

  if (sel) out.push(`<rect class="mp-sel" x="${f2(X(sel[0]))}" y="${enzRows * EH + 10}" width="${f2(Math.max(1.5, X(sel[1]) - X(sel[0])))}" height="${H - enzRows * EH - 14}"/>`);
  out.push(`<line class="mp-ring" x1="${LIN_M}" y1="${baseY}" x2="${W - LIN_M}" y2="${baseY}" stroke-width="3"/>`);
  const step = niceStep(n, Math.max(4, Math.round(W / 110)));
  for (let p = 0; p <= n; p += step / 5) {
    const major = Math.abs(p % step) < 1e-9 || Math.abs(p % step - step) < 1e-9;
    out.push(`<line class="mp-tick" x1="${f2(X(p))}" y1="${baseY + 4}" x2="${f2(X(p))}" y2="${baseY + (major ? 12 : 8)}"/>`);
    if (major) out.push(`<text class="mp-ticklbl" x="${f2(X(p))}" y="${baseY + 24}" text-anchor="middle">${Math.round(p).toLocaleString('en-US')}</text>`);
  }
  // enzymes
  for (const e of enz) {
    const y = 24 + e.row * EH + 6;
    out.push(`<line class="mp-cut${sel && e.items.some(i => sel[0] === i.s && sel[1] === i.e) ? ' on' : ''}" data-enz="${esc(e.names.join('|'))}" x1="${f2(e.x)}" y1="${y + 4}" x2="${f2(e.x)}" y2="${baseY}"/>`);
    out.push(enzGroupSvg(e, `class="mp-lbl enz" x="${f2(e.x + 3)}" y="${y}" dominant-baseline="central"`));
  }
  // features
  for (const f of feats) {
    const k = lane.get(f.id), y = featY0 + k * LH, it = info.get(f.id);
    const headSeg = f.strand === -1 ? 0 : f.locs.length - 1;
    f.locs.forEach(([a, b], i) => {
      const x1 = X(a), x2 = X(b), h = FH, hl = Math.min(11, x2 - x1);
      const head = f.strand === 0 || i !== headSeg ? 0 : f.strand;
      let d;
      if (head === 1) d = `M${f2(x1)} ${y} H${f2(x2 - hl)} L${f2(x2)} ${y + h / 2} L${f2(x2 - hl)} ${y + h} H${f2(x1)}Z`;
      else if (head === -1) d = `M${f2(x2)} ${y} H${f2(x1 + hl)} L${f2(x1)} ${y + h / 2} L${f2(x1 + hl)} ${y + h} H${f2(x2)}Z`;
      else d = `M${f2(x1)} ${y} H${f2(x2)} V${y + h} H${f2(x1)}Z`;
      out.push(`<path class="mp-feat${doc.selFid === f.id ? ' on' : ''}${f.orf ? ' orf' : ''}" data-fid="${f.id}" d="${d}" fill="${f.color}"></path>`);
    });
    for (const [a, b] of hangSegs(f, n)) {
      const x1 = X(a), x2 = X(b), cnt = Math.max(2, Math.round((x2 - x1) / 6)), pts = [];
      for (let i = 0; i <= cnt; i++) pts.push(`${f2(x1 + (x2 - x1) * i / cnt)},${f2(y + (i % 2 ? FH - 3 : 3))}`);
      out.push(`<polyline class="mp-hang${doc.selFid === f.id ? ' on' : ''}" data-fid="${f.id}" points="${pts.join(' ')}" stroke="${f.color}"/>`);
    }
    const [a, b] = featBounds(f);
    if (it.inside) out.push(`<text class="mp-flbl" data-fid="${f.id}" x="${f2((X(a) + X(b)) / 2)}" y="${y + FH / 2}" text-anchor="middle" dominant-baseline="central" fill="${textOn(f.color)}">${esc(f.name)}</text>`);
    else if (f.strand === -1) out.push(`<text class="mp-lbl" data-fid="${f.id}" x="${f2(X(a) - 4)}" y="${y + FH / 2}" text-anchor="end" dominant-baseline="central">${esc(f.name)}</text>`);
    else out.push(`<text class="mp-lbl" data-fid="${f.id}" x="${f2(X(b) + 4)}" y="${y + FH / 2}" dominant-baseline="central">${esc(f.name)}</text>`);
  }
  if (!sel) out.push(`<line class="mp-caret" x1="${f2(X(doc.caret))}" y1="${enzRows * EH + 12}" x2="${f2(X(doc.caret))}" y2="${H - 8}"/>`);
  out.push(`<text class="mp-title lin" x="${LIN_M}" y="16">${esc(doc.name)} <tspan class="mp-sub small">· ${n.toLocaleString('en-US')} bp · linear</tspan></text>`);
  host.innerHTML = `<svg class="map lin" width="${W}" height="${H}" data-n="${n}" data-w="${W}">${out.join('')}</svg>`;
}
function linearPos(svg, e, n) {
  const r = svg.getBoundingClientRect(), W = +svg.dataset.w;
  return clamp(Math.round((e.clientX - r.left - LIN_M) / (W - 2 * LIN_M) * n), 0, n);
}

function linearBase(svg, e, n) {
  const r = svg.getBoundingClientRect(), W = +svg.dataset.w, x = e.clientX - r.left - LIN_M;
  return x < 0 || x > W - 2 * LIN_M ? -1 : Math.min(n - 1, Math.floor(x / (W - 2 * LIN_M) * n));
}


/* ---------- export the map as a standalone SVG / PNG ---------- */
const EXPORT_CSS = `
svg { font-family: "Helvetica Neue", Helvetica, Arial, sans-serif; }
.mp-ring { stroke: #6a7686; opacity: .55; fill: none; }
.mp-tick { stroke: #6a7686; stroke-width: 1; opacity: .75; }
.mp-ticklbl { fill: #6a7686; font-size: 11px; }
.mp-feat { stroke: rgba(0,0,0,.28); stroke-width: .8; }
.mp-feat.orf { opacity: .55; }
.mp-lbl { fill: #1d2530; font-size: 16px; font-weight: 700; }
.mp-lbl.enz { fill: #1c55c7; font-weight: 600; font-size: 12px; }
.mp-flbl { font-size: 15px; font-weight: 600; }
.mp-hang { fill: none; stroke-width: 2; stroke-linejoin: round; cursor: pointer; }
.mp-hang.hov, .mp-hang.on { stroke: #e5322d; }
.mp-lead { fill: none; stroke: #6a7686; stroke-width: .8; opacity: .55; pointer-events: none; }
.mp-lead.on, .mp-lead.hov { stroke: #e5322d; stroke-width: 1.4; opacity: 1; }
.mp-cut { stroke: #1c55c7; stroke-width: 1; opacity: .45; pointer-events: none; }
.mp-cut.hov, .mp-cut.on { stroke: #e5322d; stroke-width: 1.8; opacity: 1; }
.mp-title { fill: #1d2530; font-size: 24px; font-weight: 700; }
.mp-title.lin { font-size: 14px; }
.mp-sub { fill: #6a7686; font-size: 17px; }
.mp-sub.small { font-size: 13px; font-weight: 400; }
`;

/* Renders the current plasmid's map offscreen (always the full map, without selection/cursor) and returns {svg, width, height} */
function buildMapSVG(doc) {
  const tmp = el('div', { style: 'position:fixed;left:-10000px;top:0;width:1400px;height:900px;visibility:hidden' });
  document.body.append(tmp);
  try {
    if (doc.circular) renderCircular(doc, tmp); else renderLinear(doc, tmp);
    const svg = tmp.querySelector('svg');
    svg.querySelectorAll('.mp-sel, .mp-caret').forEach(n => n.remove());
    if (doc.circular) svg.querySelectorAll('.mp-sub.small').forEach(n => n.remove());   // "… selected" caption
    svg.querySelectorAll('.on').forEach(n => n.classList.remove('on'));
    [svg, ...svg.querySelectorAll('*')].forEach(n => { for (const a of Array.from(n.attributes)) if (a.name.startsWith('data-')) n.removeAttribute(a.name); });
    let w, h;
    if (doc.circular) { w = CW; h = CH; svg.setAttribute('width', w); svg.setAttribute('height', h); }
    else { w = +svg.getAttribute('width'); h = +svg.getAttribute('height'); svg.setAttribute('viewBox', `0 0 ${w} ${h}`); }
    svg.removeAttribute('class'); svg.removeAttribute('preserveAspectRatio');
    svg.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
    svg.insertAdjacentHTML('afterbegin', `<style>${EXPORT_CSS}</style><rect width="100%" height="100%" fill="#ffffff"/>`);
    return { svg: new XMLSerializer().serializeToString(svg), width: w, height: h };
  } finally { tmp.remove(); }
}
function svgToPngBlob(svgText, w, h, scale = 2) {
  return new Promise((resolve, reject) => {
    const url = URL.createObjectURL(new Blob([svgText], { type: 'image/svg+xml;charset=utf-8' })), img = new Image();
    img.onload = () => {
      const c = el('canvas'); c.width = Math.round(w * scale); c.height = Math.round(h * scale);
      const g = c.getContext('2d'); g.fillStyle = '#fff'; g.fillRect(0, 0, c.width, c.height); g.drawImage(img, 0, 0, c.width, c.height);
      URL.revokeObjectURL(url); c.toBlob(b => (b ? resolve(b) : reject(new Error('PNG encoding failed'))), 'image/png');
    };
    img.onerror = () => { URL.revokeObjectURL(url); reject(new Error('Could not rasterise the map')); };
    img.src = url;
  });
}
