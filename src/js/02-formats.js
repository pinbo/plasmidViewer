'use strict';
/* ---------- file formats: GenBank, FASTA, SnapGene .dna (read) ---------- */

function splitTop(s) {
  const out = []; let depth = 0, cur = '';
  for (const c of s) {
    if (c === '(') depth++; else if (c === ')') depth--;
    if (c === ',' && depth === 0) { out.push(cur); cur = ''; } else cur += c;
  }
  if (cur) out.push(cur);
  return out;
}

function parseLocation(str) {
  str = str.replace(/\s+/g, '');
  const wrapped = /^complement\(/.test(str);
  function rec(s, st) {
    let m;
    if ((m = s.match(/^complement\((.*)\)$/))) return rec(m[1], -st);
    if ((m = s.match(/^(?:join|order)\((.*)\)$/))) { const out = []; for (const p of splitTop(m[1])) out.push(...rec(p, st)); return out; }
    const r = s.replace(/[<>]/g, '').replace(/^[^:]*:/, '');
    if ((m = r.match(/^(\d+)\.\.(\d+)$/))) return [{ s: +m[1] - 1, e: +m[2], st }];
    if ((m = r.match(/^(\d+)\^(\d+)$/))) return [{ s: +m[1] - 1, e: +m[1], st }];
    if ((m = r.match(/^(\d+)\.(\d+)$/))) return [{ s: +m[1] - 1, e: +m[2], st }];
    if ((m = r.match(/^(\d+)$/))) return [{ s: +m[1] - 1, e: +m[1], st }];
    return [];
  }
  let segs = rec(str, 1);
  if (!segs.length) return null;
  if (!wrapped && segs.length > 1 && segs.every(x => x.st === -1)) segs = segs.reverse();
  return { strand: segs[0].st, locs: segs.map(x => [x.s, x.e]) };
}

function parseGenBank(text) {
  const lines = text.replace(/\r/g, '').split('\n');
  let i = 0;
  while (i < lines.length && !/^LOCUS/.test(lines[i])) i++;
  if (i >= lines.length) throw new Error('No LOCUS line found – not a GenBank file.');
  const doc = { name: 'Untitled', circular: false, seq: '', features: [], meta: {} };
  const lm = lines[i].match(/^LOCUS\s+(\S+)/);
  if (lm) doc.name = lm[1];
  doc.circular = /\bcircular\b/i.test(lines[i]);
  let section = '', cur = null, seq = [];

  const finish = () => {
    if (!cur) return;
    const q = {};
    for (const raw of cur.quals) {
      const m = raw.match(/^\/([^=\s]+)(?:=(.*))?$/s);
      if (!m) continue;
      let v = m[2] === undefined ? '' : m[2].trim();
      if (v.startsWith('"')) v = v.replace(/^"/, '').replace(/"$/, '').replace(/""/g, '"');
      (q[m[1]] = q[m[1]] || []).push(v);
    }
    cur = { ...cur, q }; cur.done = true;
    const loc = parseLocation(cur.loc);
    if (loc && cur.key !== 'source') {
      const first = k => (q[k] && q[k][0]) || '';
      const name = first('label') || first('ApEinfo_label') || first('gene') || first('product') || first('standard_name') || (first('note') || '').slice(0, 40) || cur.key;
      let color = first('ApEinfo_fwdcolor') || first('ApEinfo_revcolor') || first('color') || first('colour');
      const keep = {};
      for (const k of Object.keys(q)) if (!/^(label|ApEinfo_|translation|color|colour)/.test(k)) keep[k] = q[k];
      // SnapGene stores display info inside /note, e.g. /note="color: #ccffcc; direction: RIGHT" – pick the colour up and keep the rest of the note
      if (keep.note) {
        keep.note = keep.note.map(t => {
          t = t.replace(/\bcolou?r:\s*#?([0-9a-fA-F]{6}|[0-9a-fA-F]{3})\b;?\s*/gi, (_, h) => { if (!/^#[0-9a-f]{6}$/i.test(color)) color = '#' + (h.length === 3 ? h.replace(/./g, c => c + c) : h); return ''; });
          return t.replace(/\bdirection:\s*\w+;?\s*/gi, '').replace(/^[\s;]+|[\s;]+$/g, '');
        }).filter(Boolean);
        if (!keep.note.length) delete keep.note;
      }
      doc.features.push({ name, type: cur.key, strand: loc.strand, locs: loc.locs, color: /^#[0-9a-f]{6}$/i.test(color) ? color.toLowerCase() : '', quals: keep });
    }
    cur = null;
  };

  for (i++; i < lines.length; i++) {
    const L = lines[i];
    if (L.startsWith('//')) break;
    if (/^ORIGIN/.test(L)) { finish(); section = 'origin'; continue; }
    if (section === 'origin') { seq.push(L.replace(/[^A-Za-z]/g, '')); continue; }
    if (/^FEATURES/.test(L)) { section = 'features'; continue; }
    if (/^\S/.test(L)) {
      finish();
      section = L.split(/\s+/)[0];
      if (section === 'DEFINITION') doc.meta.definition = L.slice(12).trim();
      continue;
    }
    if (section !== 'features') continue;
    if (/^ {5}\S/.test(L)) {
      finish();
      cur = { key: L.slice(5, 21).trim(), loc: L.slice(21).trim(), quals: [], inQ: false };
    } else if (cur && /^ {21}/.test(L)) {
      const t = L.slice(21).trim();
      if (t.startsWith('/')) { cur.quals.push(t); cur.inQ = true; }
      else if (cur.inQ) cur.quals[cur.quals.length - 1] += (cur.quals[cur.quals.length - 1].startsWith('/translation') ? '' : ' ') + t;
      else cur.loc += t;
    }
  }
  finish();
  doc.seq = cleanSeq(seq.join(''));
  if (!doc.seq) throw new Error('GenBank file has no sequence.');
  return doc;
}

function parseFasta(text) {
  const m = text.match(/^\s*>(.*)$/m);
  const name = m ? (m[1].trim().split(/\s+/)[0] || 'Untitled') : 'Untitled';
  const body = text.split(/^\s*>/m).filter(Boolean);
  const seq = cleanSeq(m ? (body.length ? body[0].replace(/^.*\n?/, '') : '') : text);
  if (!seq) throw new Error('No sequence found.');
  return { name, circular: true, seq, features: [], meta: {} };
}

function parseSnapGene(buf, fallbackName) {
  const u8 = new Uint8Array(buf), dv = new DataView(buf);
  let p = 0, seq = '', circular = false, xml = '';
  const dec = new TextDecoder('utf-8');
  while (p + 5 <= u8.length) {
    const type = u8[p], len = dv.getUint32(p + 1, false), off = p + 5;
    if (off + len > u8.length) break;
    if (type === 0) { circular = (u8[off] & 1) === 1; seq = dec.decode(u8.subarray(off + 1, off + len)); }
    else if (type === 10) xml = dec.decode(u8.subarray(off, off + len));
    p = off + len;
  }
  seq = cleanSeq(seq);
  if (!seq) throw new Error('Could not read the sequence of this SnapGene file.');
  const doc = { name: fallbackName, circular, seq, features: [], meta: {} };
  if (xml) {
    const x = new DOMParser().parseFromString(xml, 'text/xml');
    for (const f of $$('Feature', x)) {
      const type = f.getAttribute('type') || 'misc_feature';
      const dir = f.getAttribute('directionality');
      const locs = []; let color = '';
      for (const s of Array.from(f.children).filter(c => c.tagName === 'Segment')) {
        const m = (s.getAttribute('range') || '').match(/(\d+)-(\d+)/);
        if (m) { locs.push([+m[1] - 1, +m[2]]); color = color || s.getAttribute('color') || ''; }
      }
      if (!locs.length) continue;
      const quals = {};
      for (const q of Array.from(f.children).filter(c => c.tagName === 'Q')) {
        const v = q.querySelector('V'); const nm = q.getAttribute('name');
        if (nm && v) { const val = v.getAttribute('text') || v.getAttribute('int') || ''; if (val && nm !== 'label') (quals[nm] = quals[nm] || []).push(String(val).replace(/<[^>]+>/g, '')); }
      }
      doc.features.push({
        name: f.getAttribute('name') || type, type, strand: dir === '2' ? -1 : dir === '0' ? 0 : 1,
        locs, color: /^#[0-9a-f]{6}$/i.test(color) ? color : '', quals,
      });
    }
  }
  return doc;
}

function parseAny(buf, fileName) {
  const base = fileName.replace(/\.[^.]*$/, '') || 'Untitled';
  const u8 = new Uint8Array(buf);
  const magic = String.fromCharCode(...u8.subarray(5, 13));
  let doc;
  if (u8[0] === 9 && magic === 'SnapGene') doc = parseSnapGene(buf, base);
  else {
    const text = new TextDecoder('utf-8').decode(u8);
    if (/^LOCUS/m.test(text)) doc = parseGenBank(text);
    else doc = parseFasta(text);
    if (doc.name === 'Untitled' || !doc.name) doc.name = base;
  }
  doc.fileName = fileName;
  return doc;
}

/* ---------- writers ---------- */
function locString(f) {
  const parts = f.locs.map(([a, b]) => (b - a === 1 ? `${a + 1}` : `${a + 1}..${b}`));
  let s = parts.length > 1 ? `join(${parts.join(',')})` : parts[0];
  if (f.strand === -1) s = `complement(${s})`;
  return s;
}

function writeGenBank(doc) {
  const n = doc.seq.length;
  const d = new Date();
  const mon = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC'][d.getMonth()];
  const date = `${String(d.getDate()).padStart(2, '0')}-${mon}-${d.getFullYear()}`;
  const lines = [];
  lines.push('LOCUS       ' + doc.name.replace(/\s+/g, '_').slice(0, 16).padEnd(16) + ' ' + String(n).padStart(11) + ' bp    DNA     ' + (doc.circular ? 'circular' : 'linear  ') + ' SYN ' + date);
  lines.push('DEFINITION  ' + (doc.meta.definition || '.'));
  lines.push('FEATURES             Location/Qualifiers');
  lines.push('     source          1..' + n);
  lines.push('                     /mol_type="other DNA"');
  const q = (k, v) => {
    const s = v === null ? `/${k}` : `/${k}="${String(v).replace(/"/g, "'").replace(/\s*\n\s*/g, ' ')}"`;
    for (let i = 0; i < s.length; i += 58) lines.push(' '.repeat(21) + s.slice(i, i + 58));
  };
  for (const f of doc.features) {
    const key = f.type || 'misc_feature';
    const loc = locString(f);
    lines.push('     ' + key.padEnd(16) + loc);
    q('label', f.name);
    for (const [k, vals] of Object.entries(f.quals || {})) for (const v of vals) q(k, v);
    const c = f.color || typeColor(f.type);
    q('ApEinfo_fwdcolor', c); q('ApEinfo_revcolor', c);
  }
  lines.push('ORIGIN');
  for (let i = 0; i < n; i += 60) {
    const chunk = doc.seq.slice(i, i + 60).match(/.{1,10}/g).join(' ');
    lines.push(String(i + 1).padStart(9) + ' ' + chunk);
  }
  lines.push('//');
  return lines.join('\n') + '\n';
}

function writeFasta(doc) {
  return '>' + doc.name.replace(/\s+/g, '_') + ' ' + doc.seq.length + ' bp ' + (doc.circular ? 'circular' : 'linear') + '\n' + (doc.seq.match(/.{1,70}/g) || []).join('\n') + '\n';
}
