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

/* spread labels vertically so that they do not overlap */
function spreadLabels(items, gap, minY, maxY) {
  items.sort((a, b) => a.y - b.y);
  for (let i = 0; i < items.length; i++) { items[i].y = Math.max(items[i].y, i ? items[i - 1].y + gap : minY); }
  if (items.length && items[items.length - 1].y > maxY) {
    items[items.length - 1].y = maxY;
    for (let i = items.length - 2; i >= 0; i--) items[i].y = Math.min(items[i].y, items[i + 1].y - gap);
  }
}

function renderCircular(doc, host) {
  const n = doc.seq.length, feats = visibleFeatures(doc);
  const R = 215, BH = 22, LH = 25;
  const { lane, count } = assignLanes(feats.map(f => ({ id: f.id, ivs: f.locs })), n * 0.002);
  const lanes = Math.max(1, count), rTop = R + (lanes - 1) * LH + BH / 2 + 4;
  const ang = p => 2 * Math.PI * p / n;
  const sel = selRange(doc);
  const out = [];
  out.push(`<circle class="mp-ring" cx="${CCX}" cy="${CCY}" r="${R}" fill="none" stroke-width="3"/>`);

  // selection band
  if (sel) {
    const a1 = ang(sel[0]), a2 = ang(sel[1]);
    out.push(`<path class="mp-sel" d="${arcSeg(R - 34, rTop + 6, a1, a2, 0, 0)}"/>`);
  }
  // ticks
  const step = niceStep(n, 12), minor = step / 5;
  for (let p = 0; p < n; p += minor) {
    const major = Math.abs(p % step) < 1e-9 || Math.abs((p % step) - step) < 1e-9;
    const a = ang(p), [x1, y1] = polar(R - 14, a), [x2, y2] = polar(R - (major ? 22 : 18), a);
    out.push(`<line class="mp-tick" x1="${f2(x1)}" y1="${f2(y1)}" x2="${f2(x2)}" y2="${f2(y2)}"/>`);
    if (major && p > 0 && (n - p) > step * 0.35) { const [tx, ty] = polar(R - 36, a); out.push(`<text class="mp-ticklbl" x="${f2(tx)}" y="${f2(ty)}" text-anchor="middle" dominant-baseline="central">${Math.round(p).toLocaleString('en-US')}</text>`); }
  }
  // features
  const labels = [];
  for (const f of feats) {
    const k = lane.get(f.id), rc = R + k * LH, r1 = rc - BH / 2, r2 = rc + BH / 2;
    const headSeg = f.strand === -1 ? 0 : f.locs.length - 1;
    f.locs.forEach(([a, b], i) => {
      const head = f.strand === 0 || i !== headSeg ? 0 : f.strand;
      const a1 = ang(a), a2 = Math.min(ang(b), 2 * Math.PI - 0.002);
      out.push(`<path class="mp-feat${doc.selFid === f.id ? ' on' : ''}${f.orf ? ' orf' : ''}" data-fid="${f.id}" d="${arcSeg(r1, r2, a1, a2, head, 14)}" fill="${f.color}"></path>`);
    });
    // label anchor: middle of the longest segment
    const seg = f.locs.reduce((m, l) => (l[1] - l[0] > m[1] - m[0] ? l : m));
    labels.push({ kind: 'f', f, theta: ang((seg[0] + seg[1]) / 2), rr: r2 + 2, text: f.name, color: f.color });
  }
  // enzymes
  const enz = enzymeView(doc);
  for (const e of enz) {
    const a = ang(e.top), [x1, y1] = polar(R - 13, a), [x2, y2] = polar(rTop + 2, a);
    out.push(`<line class="mp-cut" x1="${f2(x1)}" y1="${f2(y1)}" x2="${f2(x2)}" y2="${f2(y2)}"/>`);
    labels.push({ kind: 'e', e, theta: a, rr: rTop + 2, text: `${e.name} (${e.top})` });
  }
  // label layout
  const Rl = rTop + 22;
  const sides = { 1: [], '-1': [] };
  for (const l of labels) { const s = Math.sin(l.theta) >= 0 ? 1 : -1; l.side = s; l.y = CCY - Rl * Math.cos(l.theta); sides[s].push(l); }
  for (const s of [1, -1]) {
    spreadLabels(sides[s], 15, 22, CH - 22);
    for (const l of sides[s]) {
      const dy = l.y - CCY, x = CCX + s * (Math.sqrt(Math.max(Rl * Rl - dy * dy, 0)) + 10);
      const [ax, ay] = polar(l.rr, l.theta);
      const cls = l.kind === 'e' ? 'mp-lbl enz' : 'mp-lbl';
      const attr = l.kind === 'e' ? `data-enz="${esc(l.e.name)}" data-s="${l.e.s}" data-e="${l.e.e}"` : `data-fid="${l.f.id}"`;
      out.push(`<polyline class="mp-lead" points="${f2(ax)},${f2(ay)} ${f2(x - s * 6)},${f2(l.y)} ${f2(x - s * 3)},${f2(l.y)}"/>`);
      out.push(`<text class="${cls}" ${attr} x="${f2(x)}" y="${f2(l.y)}" text-anchor="${s === 1 ? 'start' : 'end'}" dominant-baseline="central">${esc(l.text)}</text>`);
    }
  }
  // caret
  if (!sel) {
    const a = ang(doc.caret), [x1, y1] = polar(R - 30, a), [x2, y2] = polar(rTop + 8, a);
    out.push(`<line class="mp-caret" x1="${f2(x1)}" y1="${f2(y1)}" x2="${f2(x2)}" y2="${f2(y2)}"/>`);
  }
  // centre text
  const nm = doc.name.length > 26 ? doc.name.slice(0, 25) + '…' : doc.name;
  out.push(`<text class="mp-title" x="${CCX}" y="${CCY - 16}" text-anchor="middle">${esc(nm)}</text>`);
  out.push(`<text class="mp-sub" x="${CCX}" y="${CCY + 12}" text-anchor="middle">${n.toLocaleString('en-US')} bp</text>`);
  if (sel) out.push(`<text class="mp-sub small" x="${CCX}" y="${CCY + 36}" text-anchor="middle">${sel[0] + 1}..${sel[1]} · ${(sel[1] - sel[0]).toLocaleString('en-US')} bp selected</text>`);
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
  const enz = enzymeView(doc).map(e => ({ ...e, text: `${e.name} (${e.top})`, x: X(e.top) }));
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
    const nameW = textWidth(f.name, font) + 8;
    const [a, b] = featBounds(f); const bw = (X(b) - X(a));
    const inside = nameW <= bw - 14;
    const ivs = f.locs.map(([s, e]) => [X(s), X(e)]);
    if (!inside) { const last = f.strand === -1 ? ivs[0] : ivs[ivs.length - 1]; ivs.push(f.strand === -1 ? [last[0] - nameW, last[0]] : [last[1], last[1] + nameW]); }
    return { id: f.id, ivs, inside, nameW };
  });
  const { lane, count } = assignLanes(items.map(i => ({ id: i.id, ivs: i.ivs })), 3);
  const info = new Map(items.map(i => [i.id, i]));
  const FH = 22, LH = 28, featY0 = baseY + 46;
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
    out.push(`<line class="mp-cut" x1="${f2(e.x)}" y1="${y + 4}" x2="${f2(e.x)}" y2="${baseY}"/>`);
    out.push(`<text class="mp-lbl enz" data-enz="${esc(e.name)}" data-s="${e.s}" data-e="${e.e}" x="${f2(e.x + 3)}" y="${y}" dominant-baseline="central">${esc(e.text)}</text>`);
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
