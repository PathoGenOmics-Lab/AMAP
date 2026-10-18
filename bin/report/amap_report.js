/* AMAP's interactive report. Hand-written, no libraries: the data is the JSON
   amap_report.py embeds, and every figure is drawn here, as SVG or, where
   there are too many marks for SVG, on a canvas. */
(function () {
  'use strict';

  // ── Small helpers ──────────────────────────────────────────────────────

  const SVGNS = 'http://www.w3.org/2000/svg';

  function h(tag, props, ...kids) {
    const e = document.createElement(tag);
    setProps(e, props);
    append(e, kids);
    return e;
  }

  function s(tag, props, ...kids) {
    const e = document.createElementNS(SVGNS, tag);
    setProps(e, props);
    append(e, kids);
    return e;
  }

  function setProps(e, props) {
    if (!props) return;
    for (const k of Object.keys(props)) {
      const v = props[k];
      if (v == null || v === false) continue;
      if (k === 'class') e.setAttribute('class', v);
      else if (k === 'style' && typeof v === 'object') {
        for (const p of Object.keys(v)) e.style.setProperty(p, v[p]);
      } else if (k.startsWith('on') && typeof v === 'function') e.addEventListener(k.slice(2), v);
      else if (k === 'text') e.textContent = v;
      else e.setAttribute(k, v === true ? '' : v);
    }
  }

  function append(e, kids) {
    for (const k of kids.flat(Infinity)) {
      if (k == null || k === false) continue;
      e.append(k.nodeType ? k : document.createTextNode(String(k)));
    }
  }

  function clear(e) {
    while (e.firstChild) e.removeChild(e.firstChild);
    return e;
  }

  const nf = new Intl.NumberFormat('en-US');
  const fmt = (n) => (n == null || Number.isNaN(n) ? '' : nf.format(n));
  const pct = (x, d = 1) => (x == null || Number.isNaN(x) ? '' : x.toFixed(d) + '%');
  const kb = (bp) => (bp >= 1e6 ? (bp / 1e6).toFixed(2) + ' Mb' : bp >= 1e3 ? (bp / 1e3).toFixed(1) + ' kb' : bp + ' bp');
  const plural = (n, one, many) => fmt(n) + ' ' + (n === 1 ? one : (many || (/[^aeiou]y$/.test(one) ? one.slice(0, -1) + 'ies' : one + 's')));
  const median = (xs) => {
    const v = xs.filter((x) => x != null && !Number.isNaN(x)).sort((a, b) => a - b);
    if (!v.length) return null;
    const m = v.length >> 1;
    return v.length % 2 ? v[m] : (v[m - 1] + v[m]) / 2;
  };
  const clamp = (x, a, b) => Math.max(a, Math.min(b, x));
  const css = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();

  function unrle(text) {
    const out = [];
    if (!text) return out;
    for (const tok of text.split(',')) {
      const star = tok.lastIndexOf('*');
      if (star > 0) {
        const v = tok.slice(0, star);
        const n = +tok.slice(star + 1);
        for (let i = 0; i < n; i++) out.push(v);
      } else out.push(tok);
    }
    return out;
  }

  function download(name, text, type) {
    const blob = new Blob([text], { type: type || 'text/plain' });
    const a = h('a', { href: URL.createObjectURL(blob), download: name });
    document.body.append(a);
    a.click();
    setTimeout(() => {
      URL.revokeObjectURL(a.href);
      a.remove();
    }, 1000);
  }

  function csvCell(v) {
    const t = v == null ? '' : String(v);
    return /[",\n\r]/.test(t) ? '"' + t.replace(/"/g, '""') + '"' : t;
  }

  // An SVG as a file: the theme's colours written in, so it opens anywhere.
  function exportSvg(svg, name) {
    const clone = svg.cloneNode(true);
    const src = svg.querySelectorAll('*');
    const dst = clone.querySelectorAll('*');
    for (let i = 0; i < src.length; i++) {
      const cs = getComputedStyle(src[i]);
      for (const p of ['fill', 'stroke', 'stroke-width', 'opacity', 'font-size', 'font-weight', 'font-family']) {
        const v = cs.getPropertyValue(p);
        if (v) dst[i].setAttribute(p, v);
      }
    }
    const vb = (svg.getAttribute('viewBox') || '').split(/\s+/).map(Number);
    if (vb.length === 4) {
      clone.insertBefore(s('rect', { x: vb[0], y: vb[1], width: vb[2], height: vb[3], fill: css('--bg-elev') }), clone.firstChild);
    }
    clone.setAttribute('xmlns', SVGNS);
    download(name, '<?xml version="1.0" encoding="UTF-8"?>\n' + new XMLSerializer().serializeToString(clone), 'image/svg+xml');
  }

  // ── Tooltip ────────────────────────────────────────────────────────────

  const tip = {
    el: null,
    show(x, y, content) {
      const t = this.el;
      clear(t);
      append(t, [content]);
      t.hidden = false;
      const r = t.getBoundingClientRect();
      const vw = window.innerWidth;
      const vh = window.innerHeight;
      let left = x + 14;
      let top = y + 14;
      if (left + r.width > vw - 8) left = Math.max(8, x - r.width - 14);
      if (top + r.height > vh - 8) top = Math.max(8, y - r.height - 14);
      t.style.left = left + 'px';
      t.style.top = top + 'px';
    },
    hide() {
      if (this.el) this.el.hidden = true;
    },
  };

  // ── Data ───────────────────────────────────────────────────────────────

  async function loadData() {
    const node = document.getElementById('amap-data');
    const text = node.textContent.trim();
    if (node.dataset.encoding === 'json') return JSON.parse(text);
    if (typeof DecompressionStream === 'undefined') {
      throw new Error('This browser cannot read the report\'s compressed data. Open it in a browser from 2023 or later (Chrome or Edge 80, Firefox 113, Safari 16.4).');
    }
    const bin = atob(text);
    const bytes = new Uint8Array(bin.length);
    for (let i = 0; i < bin.length; i++) bytes[i] = bin.charCodeAt(i);
    const stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream('deflate'));
    return JSON.parse(await new Response(stream).text());
  }

  // Gene states, as amap_report.py writes them.
  const ST = {
    M: { code: 1, key: 'missing', label: 'Missing', color: '--st-missing', what: 'no gene names it' },
    P: { code: 2, key: 'pseudo', label: 'Pseudogene', color: '--st-pseudo', what: 'flagged as a pseudogene' },
    T: { code: 3, key: 'partial', label: 'Partial', color: '--st-partial', what: 'cut short by a contig end' },
    C: { code: 4, key: 'copies', label: 'More copies', color: '--st-copies', what: 'found more than once' },
    L: { code: 5, key: 'length', label: 'Length differs', color: '--st-length', what: 'a different length from the reference' },
    S: { code: 6, key: 'shared', label: 'Shared copy', color: '--st-shared', what: 'only in a gene that names several reference genes' },
  };
  const BY_CODE = [null, 'M', 'P', 'T', 'C', 'L', 'S'];
  // Where in its row the map draws each state, top to bottom as fractions
  const BAND = { M: [0, 1], T: [0, 1], P: [0.22, 0.78], C: [0, 0.3], L: [0.7, 1] };
  // Intact genes: the map's background, and AMAP's teal wherever they are counted.
  const INTACT = { label: 'Intact', color: '--amap' };
  const BAKTA_ORIGINS = new Set(['1', '2', '3', '4', 'b']);
  const ORIGIN_LABEL = {
    '2': 'case 2, perfect match',
    '3': 'case 3, imperfect match',
    '4': 'case 4, Bakta\'s call kept',
    '1': 'case 1, no reference gene',
    R: 'reference, over a Bakta call',
    A: 'reference, where Bakta called nothing',
    b: 'Bakta call',
    r: 'from the reference',
    '?': 'unknown',
  };

  // MTBC lineages in the colours the lab's figures use (mycolorsTB), with the
  // darkest lifted for the dark theme. Other groups take the categorical
  // palette in a fixed order.
  const LINEAGE_LIGHT = {
    L1: '#ff3091', L2: '#001aff', L3: '#8a0bd2', L4: '#ff0000', L5: '#995200', L6: '#1eb040',
    L7: '#c9c400', L8: '#ff9d00', L9: '#2fcf28', L10: '#6f9d81', A1: '#d1ae00', A2: '#3fbf8f',
    A3: '#3d9be0', A4: '#e86fbd', La1: '#e86fbd', La2: '#3fbf8f', La3: '#3d9be0', Untyped: '#8a8a8a', NA: '#808080',
  };
  const LINEAGE_DARK = Object.assign({}, LINEAGE_LIGHT, {
    L2: '#5b72ff', L3: '#b05be8', L4: '#ff5050', L5: '#c27a2c', L7: '#fbff00', L10: '#8fbda1',
    A2: '#8ef5c8', A3: '#73c2ff', A4: '#ff9cdb', La1: '#ff9cdb', La2: '#8ef5c8', La3: '#73c2ff',
    Untyped: '#a3a3a3', NA: '#999999',
  });
  const CATEGORICAL = ['#3987e5', '#d95926', '#199e70', '#c98500', '#d55181', '#008300', '#9085e9', '#e66767'];

  // ── Model ──────────────────────────────────────────────────────────────

  let D, REF, NREF, OK, FAILED, LINEAGES, STATUS, DETAIL, REF_PLUS, REF_INDEX, REFLEN;
  const geneCache = new Map();
  const refRowsCache = new Map();

  function buildModel(data) {
    D = data;
    REF = D.ref;
    NREF = REF.names.length;
    REF_PLUS = Array.from(REF.strand || '', (c) => c === '+');
    REF_INDEX = new Map(REF.names.map((n, i) => [n, i]));
    REFLEN = D.meta.reference.length || Math.max(0, ...REF.end.filter((x) => x != null));
    OK = D.samples.filter((x) => !x.failed);
    FAILED = D.samples.filter((x) => x.failed);
    OK.forEach((x, i) => (x.i = i));
    const seen = [];
    for (const x of OK) if (!seen.includes(x.lineage)) seen.push(x.lineage);
    LINEAGES = seen;
    let k = 0;
    const extra = {};
    for (const l of LINEAGES) if (!(l in LINEAGE_LIGHT)) extra[l] = CATEGORICAL[k++ % CATEGORICAL.length];
    Object.assign(LINEAGE_LIGHT, extra);
    Object.assign(LINEAGE_DARK, extra);

    STATUS = new Uint8Array(OK.length * NREF);
    DETAIL = new Map();
    const re = /^(\d+)([A-Z])(.*)$/;
    OK.forEach((smp, si) => {
      const text = D.status[smp.name] || '';
      let idx = 0;
      if (!text) return;
      for (const tok of text.split(',')) {
        const m = re.exec(tok);
        if (!m) continue;
        idx += +m[1];
        STATUS[si * NREF + idx] = ST[m[2]].code;
        if (m[3] !== '') DETAIL.set(si * NREF + idx, m[3]);
      }
    });
  }

  const statusAt = (si, gi) => STATUS[si * NREF + gi];
  const detailAt = (si, gi) => DETAIL.get(si * NREF + gi);
  const lineageColor = (l) => (theme() === 'light' ? LINEAGE_LIGHT : LINEAGE_DARK)[l] || '#888';
  const refName = (v) => (typeof v === 'number' ? REF.names[v] : v);
  const geneLabel = (gi) => (REF.name[gi] ? REF.name[gi] + ' (' + REF.names[gi] + ')' : REF.names[gi]);

  // The genes of one genome, decoded from amap_report.py's columns
  // (encode_genes); the two must agree.
  function genesOf(si) {
    if (geneCache.has(si)) return geneCache.get(si);
    const smp = OK[si];
    const enc = D.tables[smp.name];
    if (!enc) return null;
    const g = decodeGenes(smp.name, enc, REF, REF_PLUS);
    geneCache.set(si, g);
    return g;
  }

  // One genome's genes from amap_report.py's columns (encode_genes there; the
  // two must agree, which tests/test_amap_report.py checks with node). Each
  // gene's start is predicted from the last gene on its contig that names a
  // single reference gene, along the reference or against it, and the
  // columns hold what is left over.
  function decodeGenes(name, enc, ref, refPlus) {
    const n = enc.n;
    const C = unrle(enc.c), R = unrle(enc.r), S = unrle(enc.s), L = unrle(enc.l);
    const T = unrle(enc.t), O = unrle(enc.o), B = unrle(enc.b);
    const start = new Int32Array(n), end = new Int32Array(n), ri = new Int32Array(n);
    const contig = new Int32Array(n), plus = new Uint8Array(n), ids = new Array(n);
    const last = new Map(), prevStart = new Map();
    const nameOf = (v) => (typeof v === 'number' ? ref.names[v] : v);
    let prevNum = 0, bi = 0, replaced = 0;
    for (let row = 0; row < n; row++) {
      const key = String(row);
      const c = +C[row];
      contig[row] = c;
      const tok = R[row];
      let r = -1, pred = 0, reverse = false;
      if (tok === 'x' || tok[0] === '=') {
        if (tok !== 'x') r = +tok.slice(1);
        pred = prevStart.has(c) ? prevStart.get(c) : 0;
      } else {
        const lastC = last.get(c);
        r = lastC[0] + +tok;
        if (r > lastC[0]) pred = lastC[1] + ref.start[r] - ref.start[lastC[0]];
        else if (r < lastC[0]) {
          pred = lastC[1] + ref.end[lastC[0]] - ref.end[r];
          reverse = true;
        } else pred = lastC[1];
      }
      const st = pred + +S[row];
      const len = +L[row] + (r >= 0 ? ref.end[r] - ref.start[r] + 1 : 0);
      start[row] = st;
      // A gene across the origin of a circular sequence: its end written out
      end[row] = enc.w && key in enc.w ? +enc.w[key] : st + len - 1;
      ri[row] = r;
      const bit = +T[row];
      plus[row] = r >= 0 ? ((refPlus[r] !== reverse) === (bit === 0) ? 1 : 0) : bit === 1 ? 1 : 0;
      if (key in enc.id) {
        const lit = enc.id[key];
        ids[row] = lit[0] === '=' ? lit.slice(1) : name + '_' + lit;
      } else if (BAKTA_ORIGINS.has(O[row])) {
        // One number further on for each call the reference replaced since
        prevNum += 1 + replaced + +B[bi++];
        replaced = 0;
        ids[row] = name + '_' + String(prevNum).padStart(5, '0');
      } else {
        ids[row] = name + '_' + (key in enc.x ? nameOf(enc.x[key][0]) : ref.names[r]);
      }
      if (O[row] === 'R') replaced++;
      if (r >= 0) last.set(c, [r, st]);
      prevStart.set(c, st);
    }
    return { n, start, end, ri, contig, plus, ids, origin: O, enc, contigs: enc.cn, contigLen: enc.cl };
  }

  function rvsOf(g, row) {
    const key = String(row);
    if (key in g.enc.x) return g.enc.x[key].map(refName);
    return g.ri[row] >= 0 ? [REF.names[g.ri[row]]] : [];
  }

  // {reference gene index: [rows]} of one genome, every gene naming it.
  function refRows(si) {
    if (refRowsCache.has(si)) return refRowsCache.get(si);
    const g = genesOf(si);
    const map = new Map();
    if (g) {
      for (let row = 0; row < g.n; row++) {
        const key = String(row);
        const list = key in g.enc.x ? g.enc.x[key].filter((v) => typeof v === 'number') : g.ri[row] >= 0 ? [g.ri[row]] : [];
        for (const gi of list) {
          if (!map.has(gi)) map.set(gi, []);
          map.get(gi).push(row);
        }
      }
    }
    refRowsCache.set(si, map);
    return map;
  }

  function geneRow(si, row) {
    const g = genesOf(si);
    const key = String(row);
    const note = key in g.enc.p ? D.notes[g.enc.p[key]] : '';
    const contigLen = (g.contigLen || [])[g.contig[row]];
    const across = g.start[row] > g.end[row];
    return {
      id: g.ids[row],
      contig: g.contigs[g.contig[row]],
      start: g.start[row],
      end: g.end[row],
      strand: g.plus[row] ? '+' : '-',
      // Across the origin, from its start round to its end
      length: across ? (contigLen ? contigLen - g.start[row] + 1 + g.end[row] : null) : g.end[row] - g.start[row] + 1,
      across,
      rvs: rvsOf(g, row),
      pseudo: key in g.enc.p,
      note,
      partial: g.enc.q[key] || '',
      origin: g.origin[row],
      product: key in g.enc.pr ? D.products[g.enc.pr[key]] : '',
    };
  }

  // The Bakta calls a genome's final annotation does not keep (encode_gone).
  function goneOf(si) {
    const g = genesOf(si);
    return g && g.enc.gone ? decodeGone(OK[si].name, g.enc.gone, g) : null;
  }

  // The Bakta calls a genome's final annotation does not keep (encode_gone in
  // amap_report.py): by locus tag number, each next to the final gene in its
  // place when one overlaps it, as the offset of its ends from that gene's.
  function decodeGone(name, enc, g) {
    const K = unrle(enc.k), N = unrle(enc.num), RW = unrle(enc.row), S = unrle(enc.s), E = unrle(enc.e);
    const out = [];
    let prevNum = 0, prevRow = 0, ni = 0;
    for (let k = 0; k < enc.n; k++) {
      const key = String(k);
      let tag;
      if (key in enc.tag) tag = enc.tag[key];
      else {
        prevNum += +N[ni++];
        tag = name + '_' + String(prevNum).padStart(5, '0');
      }
      if (RW[k] === 'x') {
        const loose = enc.loose[key];
        out.push({ tag, case: K[k], row: -1, start: +S[k], end: +E[k], strand: loose[0], contig: loose[1] });
      } else {
        prevRow += +RW[k];
        out.push({
          tag, case: K[k], row: prevRow,
          start: g.start[prevRow] + +S[k], end: g.end[prevRow] + +E[k],
          strand: g.plus[prevRow] ? '+' : '-', contig: g.contigs[g.contig[prevRow]],
        });
      }
    }
    return out;
  }

  // ── State: which genomes are shown ─────────────────────────────────────

  const state = { lineages: null, text: '' };
  const listeners = [];

  function visible() {
    const q = state.text.toLowerCase();
    return OK.filter((x) => state.lineages.has(x.lineage) && (!q || x.name.toLowerCase().includes(q) || String(x.id).toLowerCase().includes(q)));
  }

  function onFilter(fn) {
    listeners.push(fn);
  }

  function filtersChanged() {
    renderFilterBar();
    for (const fn of listeners) fn();
  }

  function theme() {
    return document.documentElement.dataset.theme === 'light' ? 'light' : 'dark';
  }

  // ── Generic pieces ─────────────────────────────────────────────────────

  function lineageChip(l) {
    return h('span', { class: 'lin', style: { '--c': lineageColor(l) } }, h('span', { class: 'sw' }), l);
  }

  function stateChip(letter, extra) {
    const st = letter ? ST[letter] : INTACT;
    return h('span', { class: 'st', style: { '--c': 'var(' + st.color + ')' } },
      h('span', { class: 'sw' + (letter === 'L' ? ' low' : '') }), st.label + (extra ? ' ' + extra : ''));
  }

  function callout(kind, title, text) {
    return h('div', { class: 'callout ' + (kind || '') }, h('p', null, title ? h('b', null, title + ' ') : null, text));
  }

  function card(title, desc, tools, ...body) {
    return h('div', { class: 'card' },
      h('div', { class: 'card-head' }, h('div', null, h('h3', null, title), desc ? h('p', null, desc) : null),
        tools ? h('div', { class: 'tools' }, tools) : null),
      body);
  }

  // An SVG chart drawn at its box's width, and again when that changes, so
  // its text stays at reading size on a phone instead of shrinking with it.
  function responsive(box, drawAt) {
    let lastW = 0;
    const api = {
      redraw() {
        const w = Math.round(box.clientWidth);
        if (!w) return;
        lastW = w;
        clear(box).append(drawAt(w));
      },
    };
    new ResizeObserver(() => {
      const w = Math.round(box.clientWidth);
      if (w && Math.abs(w - lastW) >= 2) api.redraw();
    }).observe(box);
    return api;
  }

  function svgButton(getSvg, name) {
    return h('button', { class: 'mini-btn', type: 'button', onclick: () => exportSvg(getSvg(), name) }, 'SVG');
  }

  // A sortable, searchable table that pages its rows and saves them as CSV.
  function table(opts) {
    const cols = opts.columns;
    const st = { sort: opts.sort || null, dir: opts.dir || 1, q: '', page: 0 };
    const size = opts.pageSize || 25;
    const root = h('div', { class: 'table-block' });
    const count = h('span', { class: 'count' });
    const search = opts.search === false ? null : h('input', {
      class: 'field', type: 'search', placeholder: opts.placeholder || 'Search', 'aria-label': 'Search this table',
      oninput: (e) => { st.q = e.target.value.toLowerCase(); st.page = 0; draw(); },
    });
    const csvBtn = h('button', { class: 'mini-btn', type: 'button', onclick: () => saveCsv() }, 'CSV');
    root.append(h('div', { class: 'table-tools' }, count, search, opts.tools || null, csvBtn));
    const thead = h('thead');
    const tbody = h('tbody');
    root.append(h('div', { class: 'table-wrap' }, h('table', null, thead, tbody)));
    const pager = h('div', { class: 'pager' });
    root.append(pager);
    const val = (c, r) => (c.value ? c.value(r) : r[c.key]);
    const shown = cols.filter((c) => !c.hidden);

    function rows() {
      let out = opts.rows;
      if (st.q) {
        out = out.filter((r) => cols.some((c) => {
          const v = val(c, r);
          return v != null && String(v).toLowerCase().includes(st.q);
        }));
      }
      if (st.sort != null) {
        const c = cols.find((x) => x.key === st.sort);
        if (c) {
          out = out.slice().sort((a, b) => {
            const x = val(c, a), y = val(c, b);
            if (x == null || x === '') return 1;
            if (y == null || y === '') return -1;
            return (typeof x === 'number' && typeof y === 'number' ? x - y : String(x).localeCompare(String(y), 'en', { numeric: true })) * st.dir;
          });
        }
      }
      return out;
    }

    function head() {
      clear(thead);
      thead.append(h('tr', null, shown.map((c) => {
        const sorted = st.sort === c.key;
        return h('th', { class: c.num ? 'num' : null, scope: 'col', 'aria-sort': sorted ? (st.dir > 0 ? 'ascending' : 'descending') : null, title: c.title || null },
          c.nosort ? c.label : h('button', { type: 'button', onclick: () => {
            if (st.sort === c.key) st.dir = -st.dir;
            else { st.sort = c.key; st.dir = c.num ? -1 : 1; }
            st.page = 0;
            draw();
          } }, c.label, h('span', { class: 'arrow', 'aria-hidden': 'true' }, sorted ? (st.dir > 0 ? '▲' : '▼') : '↕')));
      })));
    }

    function draw() {
      head();
      const all = rows();
      const pages = Math.max(1, Math.ceil(all.length / size));
      st.page = clamp(st.page, 0, pages - 1);
      clear(tbody);
      const slice = all.slice(st.page * size, (st.page + 1) * size);
      if (!slice.length) {
        tbody.append(h('tr', { class: 'empty-row' }, h('td', { colspan: shown.length }, opts.empty || 'Nothing to show.')));
      }
      for (const r of slice) {
        const tr = h('tr', { class: [opts.onRow ? 'clickable' : '', opts.rowClass ? opts.rowClass(r) : ''].join(' ').trim() || null,
          tabindex: opts.onRow ? '0' : null });
        if (opts.onRow) {
          tr.addEventListener('click', () => opts.onRow(r));
          tr.addEventListener('keydown', (e) => { if (e.key === 'Enter') opts.onRow(r); });
        }
        for (const c of shown) {
          const v = c.render ? c.render(r) : val(c, r);
          tr.append(h('td', { class: [c.num ? 'num' : '', c.wrap ? 'wrap-cell' : '', c.cellClass ? c.cellClass(r) : ''].join(' ').trim() || null, title: c.cellTitle ? c.cellTitle(r) : null },
            c.num && typeof v === 'number' ? fmt(v) : v));
        }
        tbody.append(tr);
      }
      count.textContent = all.length === opts.rows.length ? plural(opts.rows.length, opts.unit || 'row') : fmt(all.length) + ' of ' + plural(opts.rows.length, opts.unit || 'row');
      clear(pager);
      if (pages > 1) {
        pager.append(
          h('button', { class: 'mini-btn', type: 'button', disabled: st.page === 0, onclick: () => { st.page--; draw(); } }, 'Previous'),
          h('span', null, 'Page ' + (st.page + 1) + ' of ' + pages),
          h('button', { class: 'mini-btn', type: 'button', disabled: st.page >= pages - 1, onclick: () => { st.page++; draw(); } }, 'Next'));
      }
    }

    function saveCsv() {
      const out = [cols.map((c) => csvCell(c.csvLabel || c.label)).join(',')];
      for (const r of rows()) out.push(cols.map((c) => csvCell(c.csv ? c.csv(r) : val(c, r))).join(','));
      download(opts.csv || 'table.csv', out.join('\n') + '\n', 'text/csv');
    }

    draw();
    root.redraw = (newRows) => {
      if (newRows) opts.rows = newRows;
      draw();
    };
    return root;
  }

  // A canvas sized to its box in CSS pixels, drawn at the screen's density.
  function sizeCanvas(canvas, width, height) {
    const dpr = window.devicePixelRatio || 1;
    canvas.width = Math.max(1, Math.round(width * dpr));
    canvas.height = Math.max(1, Math.round(height * dpr));
    canvas.style.height = height + 'px';
    const ctx = canvas.getContext('2d');
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    return ctx;
  }

  function colorOf(letter) {
    return css(letter ? ST[letter].color : '--st-intact');
  }

  // ── Filter bar ─────────────────────────────────────────────────────────

  function renderFilterBar() {
    const pills = document.getElementById('lineage-pills');
    clear(pills);
    for (const l of LINEAGES) {
      const n = OK.filter((x) => x.lineage === l).length;
      const on = state.lineages.has(l);
      pills.append(h('button', {
        class: 'pill', type: 'button', 'aria-pressed': on ? 'true' : 'false', style: { '--c': lineageColor(l) },
        title: (on ? 'Hide ' : 'Show ') + l,
        onclick: (e) => {
          if (e.altKey || e.metaKey) state.lineages = new Set([l]);
          else if (on) state.lineages.delete(l);
          else state.lineages.add(l);
          filtersChanged();
        },
      }, h('span', { class: 'sw' }), l, h('span', { class: 'n' }, String(n))));
    }
    const v = visible().length;
    document.getElementById('filters-count').textContent = v === OK.length ? plural(OK.length, 'genome') + ', all shown' : fmt(v) + ' of ' + plural(OK.length, 'genome') + ' shown';
    document.getElementById('filters-reset').hidden = v === OK.length;
  }

  // ── Hero ───────────────────────────────────────────────────────────────

  function renderHero() {
    const m = D.meta;
    const ref = m.reference;
    const refLabel = ref.custom === true ? 'a custom reference (' + baseName(ref.annotation) + ')' : ref.custom === false ? 'H37Rv' : 'the reference';
    document.getElementById('hero-eyebrow').textContent = 'AMAP ' + (m.version ? 'v' + String(m.version).replace(/^v/, '') : '') + ' · report written ' + m.generated;
    const found = OK.map((x) => x.found_pct);
    const lo = found.length ? Math.min(...found) : 0, hi = found.length ? Math.max(...found) : 0;
    const flagged = OK.filter((x) => x.flags && x.flags.length).length;
    const lead = [
      plural(OK.length, 'genome') + (LINEAGES.length && m.modes.typing ? ' of ' + plural(LINEAGES.length, 'lineage') : '') +
        ', annotated against ' + refLabel + ' (' + plural(m.reference.genes, 'reference gene') + ').',
      OK.length ? (OK.length === 1 ? ' It has ' + pct(lo) : ' Each has between ' + pct(lo) + ' and ' + pct(hi)) + ' of the reference genes.' : '',
      flagged ? ' ' + plural(flagged, 'genome') + ' stand' + (flagged === 1 ? 's' : '') + ' out from the rest in section 1.' : '',
      FAILED.length ? ' ' + plural(FAILED.length, 'genome') + ' produced no annotation.' : '',
    ];
    document.getElementById('hero-lead').textContent = lead.join('');
    const tiles = document.getElementById('hero-tiles');
    const flows = OK.filter((x) => x.flow);
    // Bakta's calls the final annotation does not keep, as they were called
    const changed = flows.length ? median(flows.map((x) => {
      const f = x.flow;
      const total = f.c1 + f.c2 + f.c3 + f.c4;
      return total ? 100 * (total - f.k1 - f.k2 - f.k3 - f.k4) / total : null;
    })) : null;
    const items = [
      [fmt(OK.length), 'genomes annotated'],
      [pct(median(found), 1), 'of the reference genes, median genome'],
      [fmt(Math.round(median(OK.map((x) => x.genes)))), 'genes per genome, median'],
      changed != null ? [pct(changed, 0), 'of Bakta\'s calls replaced or removed, median'] : null,
      [fmt(Math.round(median(OK.map((x) => x.pseudo)))), 'pseudogenes per genome, median'],
    ].filter(Boolean);
    clear(tiles);
    for (const [b, t] of items) tiles.append(h('div', { class: 'tile' }, h('b', null, b), h('span', null, t)));
    const modes = document.getElementById('hero-modes');
    clear(modes);
    const mode = (on, label) => h('span', { class: 'mode ' + (on ? 'on' : 'off') }, h('i', null, on ? '✓' : '–'), label);
    modes.append(
      mode(ref.custom != null, 'Reference: ' + (ref.custom === true ? 'custom' : ref.custom === false ? 'H37Rv, bundled' : 'not recorded')),
      mode(m.modes.typing, m.modes.typing ? 'Lineages typed' : 'No lineage typing'),
      mode(m.modes.cases, m.modes.cases ? 'Case tables' : 'No case tables'),
      mode(m.modes.functional, m.modes.functional ? 'Functional annotation (eggNOG)' : 'No functional annotation'),
      mode(m.modes.mycobrowser, m.modes.mycobrowser ? 'Mycobrowser categories' : 'No Mycobrowser categories'),
      mode(!!D.pangene, D.pangene ? 'Pangenome' : 'No pangenome'));
    const alerts = document.getElementById('hero-alerts');
    clear(alerts);
    if (FAILED.length || (m.expected_missing && m.expected_missing.length)) {
      const names = FAILED.map((x) => x.id).concat(m.expected_missing || []);
      alerts.append(h('div', { class: 'callout warn', style: { 'margin-top': '18px' } },
        h('p', null, h('b', null, plural(names.length, 'genome') + ' of the sample sheet ' + (names.length === 1 ? 'has' : 'have') + ' no annotation: '),
          names.join(', ') + '. The pipeline ignores a genome whose tasks keep failing and finishes without it; ',
          h('code', null, 'pipeline_info/execution_trace_*.txt'), ' says which task failed.')));
    }
    for (const w of m.warnings || []) alerts.append(h('div', { class: 'callout warn', style: { 'margin-top': '10px' } }, h('p', null, w)));
    document.getElementById('foot-generated').textContent = 'Written ' + m.generated + (m.version ? ' by AMAP v' + String(m.version).replace(/^v/, '') : '') + '.';
  }

  // ── 1. Genomes ─────────────────────────────────────────────────────────

  const QC_NAMES = {
    size: 'genome size', contigs: 'contigs', n50: 'N50', found_pct: 'reference genes found', missing: 'missing genes',
    pseudo: 'pseudogenes', partial: 'partial genes', noref: 'genes not in the reference',
  };

  function flagText(f) {
    const v = f.key === 'size' || f.key === 'n50' ? kb(f.value) : f.key === 'found_pct' ? pct(f.value, 2) : fmt(f.value);
    const m = f.key === 'size' || f.key === 'n50' ? kb(f.median) : f.key === 'found_pct' ? pct(f.median, 2) : fmt(f.median);
    return QC_NAMES[f.key] + ' ' + v + ', against a median of ' + m + ' in ' + f.group +
      (f.z != null ? ' (' + Math.abs(f.z) + ' robust SD ' + (f.z > 0 ? 'above' : 'below') + ')' : ', at least half of which have exactly that value');
  }

  function renderSamples() {
    const body = clear(document.getElementById('samples-body'));
    const rows = visible();
    const max = (k) => Math.max(1, ...rows.map((x) => x[k] || 0));
    const mx = { missing: max('missing'), pseudo: max('pseudo'), partial: max('partial'), noref: max('noref'), copies: max('copies') };
    const flagged = (r, k) => r.flags && r.flags.some((f) => f.key === k);
    const bar = (r, k, color) => h('span', { class: 'bar', style: { width: Math.round(36 * (r[k] || 0) / mx[k]) + 'px', '--c': 'var(' + color + ')' }, 'aria-hidden': 'true' });
    const cols = [
      { key: 'name', label: 'Genome', value: (r) => r.id, csv: (r) => r.id, title: 'The sample ID; its files are named after the lineage and the ID',
        render: (r) => h('span', null, h('b', null, r.id), (r.flags || []).length ? h('span', { class: 'flag', style: { 'margin-left': '8px' }, title: r.flags.map(flagText).join('\n') }, '⚑ ' + r.flags.length) : null) },
      { key: 'file', label: 'Name in the outputs', value: (r) => r.name, hidden: true },
      { key: 'lineage', label: 'Lineage', render: (r) => lineageChip(r.lineage), value: (r) => r.lineage },
      { key: 'size', label: 'Size', num: true, value: (r) => r.size, render: (r) => (r.size ? (r.size / 1e6).toFixed(2) + ' Mb' : ''), cellClass: (r) => (flagged(r, 'size') ? 'flagged' : ''), title: 'Genome size' },
      { key: 'contigs', label: 'Contigs', num: true, cellClass: (r) => (flagged(r, 'contigs') ? 'flagged' : '') },
      { key: 'n50', label: 'N50', num: true, render: (r) => (r.n50 ? kb(r.n50) : ''), cellClass: (r) => (flagged(r, 'n50') ? 'flagged' : '') },
      { key: 'bakta', label: 'Bakta calls', num: true, title: 'CDS Bakta called' },
      { key: 'genes', label: 'AMAP genes', num: true, title: 'Genes in the final annotation' },
      { key: 'found_pct', label: 'Reference found', num: true, render: (r) => pct(r.found_pct, 2), cellClass: (r) => (flagged(r, 'found_pct') ? 'flagged' : ''), title: 'Reference genes some gene of this genome names' },
      { key: 'missing', label: 'Missing', num: true, render: (r) => h('span', null, fmt(r.missing), bar(r, 'missing', '--st-missing')), cellClass: (r) => (flagged(r, 'missing') ? 'flagged' : ''), title: 'Reference genes no gene names' },
      { key: 'pseudo', label: 'Pseudogenes', num: true, render: (r) => h('span', null, fmt(r.pseudo), bar(r, 'pseudo', '--st-pseudo')), cellClass: (r) => (flagged(r, 'pseudo') ? 'flagged' : '') },
      { key: 'partial', label: 'Partial', num: true, render: (r) => h('span', null, fmt(r.partial), bar(r, 'partial', '--st-partial')), cellClass: (r) => (flagged(r, 'partial') ? 'flagged' : ''), title: 'Genes a contig end cuts short' },
      { key: 'copies', label: 'More copies', num: true, render: (r) => h('span', null, fmt(r.copies), bar(r, 'copies', '--st-copies')), title: 'Reference genes found more than once' },
      { key: 'noref', label: 'Not in reference', num: true, render: (r) => h('span', null, fmt(r.noref), bar(r, 'noref', '--accent')), cellClass: (r) => (flagged(r, 'noref') ? 'flagged' : ''), title: 'Genes that name no reference gene' },
      { key: 'flags', label: 'Stands out', value: (r) => (r.flags || []).length, csv: (r) => (r.flags || []).map(flagText).join('; '),
        render: (r) => (r.flags || []).map((f) => h('span', { class: 'flag', title: flagText(f) }, '⚑ ' + QC_NAMES[f.key])) },
    ];
    if (!D.meta.modes.genomes) cols.splice(cols.findIndex((c) => c.key === 'size'), 3);
    body.append(table({
      columns: cols, rows, sort: 'lineage', dir: 1, pageSize: 50, unit: 'genome', csv: 'amap_genomes.csv',
      placeholder: 'Search genomes', onRow: (r) => openSample(r.i), empty: 'No genome matches the filters.',
    }));
    const notes = [];
    if (OK.length < 4) notes.push(callout('absent', 'Too few genomes to compare.', 'A value is flagged when it lies far from comparable genomes, which needs four genomes or more.'));
    else notes.push(callout('', 'How a value is flagged.', 'Each genome is compared with its own lineage when that lineage has five genomes or more, and with every other genome otherwise. A value is flagged when it lies more than 3.5 robust standard deviations (1.4826 times the median absolute deviation) from their median on the side that is bad, and beyond a minimum gap, so that a difference of a gene or two between identical genomes is not flagged.'));
    if (!D.meta.modes.genomes) notes.push(callout('absent', 'No genome files.', 'The report was written without the genomes, so it has no sizes, contigs or N50.'));
    for (const n of notes) {
      n.style.marginTop = '14px';
      body.append(n);
    }
  }

  // ── The genome drawer ──────────────────────────────────────────────────

  let drawerReturn = null;

  function openSample(si) {
    const smp = OK[si];
    const drawer = document.getElementById('drawer');
    drawerReturn = document.activeElement;
    document.getElementById('drawer-title').textContent = smp.id;
    const body = clear(document.getElementById('drawer-body'));
    const g = genesOf(si);
    const facts = h('div', { class: 'facts' },
      [['Lineage', lineageChip(smp.lineage)], ['Name in the outputs', smp.name],
        ['Genome', smp.size ? (smp.size / 1e6).toFixed(3) + ' Mb, ' + plural(smp.contigs, 'contig') : 'not given'],
        ['N50', smp.n50 ? kb(smp.n50) : ''], ['GC', smp.gc != null ? pct(smp.gc, 2) : ''],
        ['Bakta calls', fmt(smp.bakta)], ['AMAP genes', fmt(smp.genes)],
        ['Reference genes found', fmt(smp.found) + ' of ' + fmt(NREF) + ' (' + pct(smp.found_pct, 2) + ')']]
        .map(([k, v]) => h('div', null, h('span', { class: 'k' }, k), v)));
    body.append(h('div', { class: 'card' }, facts,
      (smp.flags || []).length ? h('div', null, smp.flags.map((f) => h('p', { style: { margin: '6px 0 0' } }, h('span', { class: 'flag' }, '⚑ ' + QC_NAMES[f.key]), ' ', flagText(f) + '.'))) : null,
      h('div', { class: 'tools', style: { 'margin-top': '12px' } },
        h('button', { class: 'mini-btn', type: 'button', onclick: () => { closeDrawer(); mapApi.highlight(si); location.hash = '#map'; } }, 'Show on the map'),
        smp.flow ? h('button', { class: 'mini-btn', type: 'button', onclick: () => { closeDrawer(); correctionApi.select(si); location.hash = '#correction'; } }, 'Show its correction') : null)));

    // Its genes by state
    const lists = [];
    const byState = { M: [], P: [], T: [], C: [], L: [], S: [] };
    for (let gi = 0; gi < NREF; gi++) {
      const c = statusAt(si, gi);
      if (c) byState[BY_CODE[c]].push(gi);
    }
    const rr = refRows(si);
    const refCols = [
      { key: 'locus', label: 'Locus', value: (r) => REF.names[r.gi] },
      { key: 'name', label: 'Gene', value: (r) => REF.name[r.gi] || '' },
      { key: 'product', label: 'Product', value: (r) => REF.product[r.gi] || '', wrap: true },
    ];
    const coordCols = [
      { key: 'id', label: 'Locus tag', value: (r) => (r.row != null ? g.ids[r.row] : '') },
      { key: 'where', label: 'Where', value: (r) => (r.row != null ? g.contigs[g.contig[r.row]] + ':' + fmt(g.start[r.row]) + '-' + fmt(g.end[r.row]) + ' ' + (g.plus[r.row] ? '+' : '-') : '') },
    ];
    const refPos = { key: 'refpos', label: 'In the reference', value: (r) => (REF.start[r.gi] != null ? fmt(REF.start[r.gi]) + '-' + fmt(REF.end[r.gi]) : '') };
    const rowsFor = (letter) => byState[letter].map((gi) => ({ gi, row: (rr.get(gi) || [])[0], detail: detailAt(si, gi) }));
    lists.push(['Missing', byState.M.length, () => table({
      columns: refCols.concat([refPos]), rows: rowsFor('M'), unit: 'gene', csv: smp.name + '_missing.csv', pageSize: 15,
      onRow: (r) => { closeDrawer(); lookupApi.show(r.gi); },
    })]);
    // Pseudogenes and partial genes as the annotation has them, so the counts
    // are the table's: a copy or a gene outside the reference counts too.
    const flagged = (test) => {
      const out = [];
      if (g) for (let row = 0; row < g.n; row++) if (test(String(row))) out.push(geneRow(si, row));
      return out;
    };
    const annCols = [
      { key: 'id', label: 'Locus tag' },
      { key: 'rvs', label: 'Reference gene', value: (r) => r.rvs.map((x) => { const gi = REF_INDEX.get(x); return gi != null && REF.name[gi] ? REF.name[gi] + ' (' + x + ')' : x; }).join(', '), wrap: true },
      { key: 'where', label: 'Where', value: (r) => r.contig + ':' + fmt(r.start) + '-' + fmt(r.end) + ' ' + r.strand },
    ];
    const pseudoRows = flagged((key) => key in g.enc.p);
    lists.push(['Pseudogenes', pseudoRows.length, () => table({
      columns: annCols.concat([{ key: 'note', label: 'Why', wrap: true }]),
      rows: pseudoRows, unit: 'gene', csv: smp.name + '_pseudogenes.csv', pageSize: 15,
      onRow: (r) => { const gi = REF_INDEX.get(r.rvs[0]); if (gi != null) { closeDrawer(); lookupApi.show(gi, si); } },
    })]);
    const partialRows = flagged((key) => key in g.enc.q);
    lists.push(['Partial', partialRows.length, () => table({
      columns: annCols.concat([{ key: 'partial', label: 'Cut at', value: (r) => r.partial.replace('\'', '') + '\' end' }]),
      rows: partialRows, unit: 'gene', csv: smp.name + '_partial.csv', pageSize: 15,
      onRow: (r) => { const gi = REF_INDEX.get(r.rvs[0]); if (gi != null) { closeDrawer(); lookupApi.show(gi, si); } },
    })]);
    lists.push(['More copies', byState.C.length, () => table({
      columns: refCols.concat([{ key: 'copies', label: 'Copies', num: true, value: (r) => +r.detail },
        { key: 'tags', label: 'Locus tags', value: (r) => (rr.get(r.gi) || []).map((row) => g.ids[row]).join(', '), wrap: true }]),
      rows: rowsFor('C'), unit: 'gene', csv: smp.name + '_copies.csv', pageSize: 15,
      onRow: (r) => { closeDrawer(); lookupApi.show(r.gi); },
    })]);
    lists.push(['Length differs', byState.L.length, () => table({
      columns: refCols.concat(coordCols, [{ key: 'delta', label: 'Difference', num: true, value: (r) => +r.detail, render: (r) => (+r.detail > 0 ? '+' : '') + fmt(+r.detail) + ' bp' }]),
      rows: rowsFor('L'), unit: 'gene', csv: smp.name + '_length.csv', pageSize: 15,
      onRow: (r) => { closeDrawer(); lookupApi.show(r.gi); },
    })]);
    const noref = [];
    if (g) for (let row = 0; row < g.n; row++) if (rvsOf(g, row).length === 0) noref.push(geneRow(si, row));
    lists.push(['Not in the reference', noref.length, () => table({
      columns: [{ key: 'id', label: 'Locus tag' }, { key: 'contig', label: 'Contig' },
        { key: 'start', label: 'Start', num: true }, { key: 'end', label: 'End', num: true }, { key: 'strand', label: 'Strand' },
        { key: 'length', label: 'Length', num: true }, { key: 'product', label: 'Product', wrap: true },
        { key: 'origin', label: 'Origin', value: (r) => ORIGIN_LABEL[r.origin] || r.origin }],
      rows: noref, unit: 'gene', csv: smp.name + '_not_in_reference.csv', pageSize: 15,
    })]);
    if (g && g.contigs.length > 1) {
      const counts = new Map();
      for (let row = 0; row < g.n; row++) counts.set(g.contig[row], (counts.get(g.contig[row]) || 0) + 1);
      lists.push(['Contigs', g.contigs.length, () => table({
        columns: [{ key: 'name', label: 'Contig' }, { key: 'length', label: 'Length', num: true }, { key: 'genes', label: 'Genes', num: true }],
        rows: g.contigs.map((c, i) => ({ name: c, length: (g.contigLen || [])[i], genes: counts.get(i) || 0 })),
        unit: 'contig', csv: smp.name + '_contigs.csv', pageSize: 15,
      })]);
    }
    const tabs = h('div', { class: 'tabs', role: 'tablist' });
    const panel = h('div', { role: 'tabpanel' });
    lists.forEach(([label, n, make], k) => {
      const b = h('button', { type: 'button', role: 'tab', 'aria-selected': k === 0 ? 'true' : 'false',
        onclick: () => {
          for (const x of tabs.children) x.setAttribute('aria-selected', 'false');
          b.setAttribute('aria-selected', 'true');
          clear(panel).append(make());
        } }, label, h('span', { class: 'n' }, fmt(n)));
      tabs.append(b);
    });
    clear(panel).append(lists[0][2]());
    body.append(h('div', { class: 'card' }, h('h3', { style: { 'margin-bottom': '10px' } }, 'Its genes, by state'), tabs, panel));
    drawer.classList.add('open');
    drawer.setAttribute('aria-hidden', 'false');
    document.getElementById('drawer-close').focus();
  }

  function closeDrawer() {
    const drawer = document.getElementById('drawer');
    drawer.classList.remove('open');
    drawer.setAttribute('aria-hidden', 'true');
    if (drawerReturn && drawerReturn.focus) drawerReturn.focus();
  }

  // ── 2. Reference map ───────────────────────────────────────────────────

  const mapApi = { highlight() {}, zoomTo() {} };

  function renderMapSection() {
    const body = clear(document.getElementById('map-body'));
    const placed = [];
    for (let gi = 0; gi < NREF; gi++) if (REF.start[gi] != null) placed.push(gi);
    if (!placed.length) {
      body.append(callout('absent', 'No reference coordinates.', 'The reference proteome\'s headers give no coordinates (NAME_START-END), so the genes cannot be placed.'));
      return;
    }
    const enabled = { M: true, P: true, T: true, C: true, L: false };
    const view = { x0: 0, x1: REFLEN, highlight: -1 };
    const where = h('span', { class: 'where' });
    const goto = h('input', { class: 'field', type: 'search', placeholder: 'Go to a gene or 1.2-1.3 Mb', 'aria-label': 'Go to a gene or a region', style: { width: '200px' } });
    const legend = h('div', { class: 'legend' });
    for (const letter of ['M', 'P', 'T', 'C', 'L']) {
      const b = h('button', { class: 'st', type: 'button', 'aria-pressed': enabled[letter] ? 'true' : 'false', style: { '--c': 'var(' + ST[letter].color + ')' },
        title: ST[letter].what,
        onclick: () => { enabled[letter] = !enabled[letter]; b.setAttribute('aria-pressed', enabled[letter] ? 'true' : 'false'); draw(); } },
      h('span', { class: 'sw band-' + letter }), ST[letter].label);
      legend.append(b);
    }
    const zoomBy = (f) => {
      const c = (view.x0 + view.x1) / 2;
      const half = clamp((view.x1 - view.x0) * f / 2, 1500, REFLEN / 2);
      setView(c - half, c + half);
    };
    const pan = (f) => {
      const w = view.x1 - view.x0;
      setView(view.x0 + w * f, view.x1 + w * f);
    };
    const toolbar = h('div', { class: 'map-toolbar' }, where,
      h('div', { class: 'tools' },
        h('button', { class: 'mini-btn', type: 'button', onclick: () => pan(-0.5), 'aria-label': 'Move left' }, '◀'),
        h('button', { class: 'mini-btn', type: 'button', onclick: () => zoomBy(0.5), 'aria-label': 'Zoom in' }, '+'),
        h('button', { class: 'mini-btn', type: 'button', onclick: () => zoomBy(2), 'aria-label': 'Zoom out' }, '−'),
        h('button', { class: 'mini-btn', type: 'button', onclick: () => pan(0.5), 'aria-label': 'Move right' }, '▶'),
        h('button', { class: 'mini-btn', type: 'button', onclick: () => setView(0, REFLEN) }, 'Whole genome'),
        goto,
        h('button', { class: 'mini-btn', type: 'button', onclick: () => exportMap() }, 'SVG')));
    const canvas = h('canvas', { tabindex: '0', 'aria-label': 'Map of the reference genes in every genome. Drag to zoom into a region; arrow keys move, plus and minus zoom.' });
    const brush = h('div', { style: { position: 'absolute', top: '0', bottom: '0', background: 'color-mix(in srgb, var(--accent) 18%, transparent)', border: '1px solid var(--accent)', display: 'none', 'pointer-events': 'none' } });
    const scroller = h('div', { class: 'canvas-scroll' }, canvas);
    const fig = h('div', { class: 'canvas-fig', style: { position: 'relative' } }, scroller, brush);
    const regionBox = h('div');
    const blocksBox = h('div');
    body.append(h('div', { class: 'card' },
      toolbar, h('div', { style: { 'margin-bottom': '10px' } }, legend), fig,
      h('p', { class: 'map-hint' }, 'Hover for a gene, click it for that gene in every genome, drag across to zoom. The strip above the rows counts, per gene, how many of the genomes shown have it in each state.')));
    body.append(regionBox, blocksBox);

    let layout = null;

    function rowsLayout() {
      const vis = visible();
      const width = Math.max(320, fig.clientWidth || 800);
      const narrow = width < 640;
      const n = vis.length;
      const rowH = clamp(Math.floor(480 / Math.max(1, n)), 3, 18);
      const gap = rowH >= 8 ? 2 : rowH >= 5 ? 1 : 0;
      const labels = rowH >= 10;
      const gutter = labels ? (narrow ? 96 : 150) : (narrow ? 34 : 46);
      const axisH = 22, geneH = 0, sumH = 36, top = axisH + sumH + 10;
      let y = top;
      const rows = [];
      let prevL = null;
      for (const smp of vis) {
        if (prevL !== null && smp.lineage !== prevL) y += 6;
        rows.push({ smp, y, h: rowH - gap });
        y += rowH;
        prevL = smp.lineage;
      }
      return { vis, width, gutter, plotW: width - gutter - 10, axisH, geneH, sumH, top, rows, height: Math.max(y + 8, top + 40), labels, rowH };
    }

    const xOf = (pos) => layout.gutter + (pos - view.x0) / (view.x1 - view.x0) * layout.plotW;
    const posOf = (x) => view.x0 + (x - layout.gutter) / layout.plotW * (view.x1 - view.x0);

    // Reference genes in the view, by binary search over the placed genes.
    function inView(x0, x1) {
      let lo = 0, hi = placed.length;
      while (lo < hi) {
        const mid = (lo + hi) >> 1;
        if (REF.end[placed[mid]] < x0) lo = mid + 1;
        else hi = mid;
      }
      const out = [];
      for (let k = Math.max(0, lo - 30); k < placed.length; k++) {
        const gi = placed[k];
        if (REF.start[gi] > x1) break;
        if (REF.end[gi] >= x0) out.push(gi);
      }
      return out;
    }

    function draw() {
      layout = rowsLayout();
      const L = layout;
      const ctx = sizeCanvas(canvas, L.width, L.height);
      ctx.clearRect(0, 0, L.width, L.height);
      ctx.font = '11px ' + css('--font');
      const ink = css('--ink-soft'), muted = css('--muted'), line = css('--line-strong'), intact = css('--st-intact');
      const colors = { M: colorOf('M'), P: colorOf('P'), T: colorOf('T'), C: colorOf('C'), L: colorOf('L'), S: colorOf('S') };
      const genes = inView(view.x0, view.x1);
      const span = view.x1 - view.x0;

      // Axis
      ctx.fillStyle = muted;
      ctx.strokeStyle = line;
      ctx.lineWidth = 1;
      const step = niceStep(span / Math.max(2, L.plotW / 110));
      ctx.beginPath();
      ctx.moveTo(L.gutter, L.axisH - 0.5);
      ctx.lineTo(L.gutter + L.plotW, L.axisH - 0.5);
      ctx.stroke();
      for (let t = Math.ceil(view.x0 / step) * step; t <= view.x1; t += step) {
        const x = Math.round(xOf(t)) + 0.5;
        ctx.beginPath();
        ctx.moveTo(x, L.axisH - 5);
        ctx.lineTo(x, L.axisH);
        ctx.stroke();
        ctx.textAlign = 'center';
        ctx.fillText(axisLabel(t, step), x, L.axisH - 9);
      }

      // Summary strip: per gene, the share of the shown genomes in each state
      const n = L.vis.length;
      const sumTop = L.axisH + 4, sumH = L.sumH;
      ctx.fillStyle = intact;
      ctx.fillRect(L.gutter, sumTop, L.plotW, sumH);
      if (n) {
        const cols = new Map();
        for (const gi of genes) {
          const x0 = Math.floor(xOf(REF.start[gi]));
          const x1 = Math.max(x0 + 1, Math.ceil(xOf(REF.end[gi] + 1)));
          const cnt = { M: 0, P: 0, T: 0, C: 0, L: 0 };
          for (const smp of L.vis) {
            const c = BY_CODE[statusAt(smp.i, gi)];
            if (c && enabled[c]) cnt[c]++;
          }
          const tot = cnt.M + cnt.P + cnt.T + cnt.C + cnt.L;
          if (!tot) continue;
          for (let x = Math.max(x0, L.gutter); x < Math.min(x1, L.gutter + L.plotW); x++) {
            const prev = cols.get(x);
            if (!prev || prev.tot < tot) cols.set(x, { tot, cnt });
          }
        }
        for (const [x, v] of cols) {
          let y = sumTop + sumH;
          for (const letter of ['M', 'P', 'T', 'C', 'L']) {
            const hh = v.cnt[letter] / n * sumH;
            if (!hh) continue;
            ctx.fillStyle = colors[letter];
            ctx.fillRect(x, y - hh, 1, hh);
            y -= hh;
          }
        }
      }
      ctx.fillStyle = muted;
      ctx.textAlign = 'right';
      ctx.fillText('all', L.gutter - 6, sumTop + 10);
      ctx.fillText('shown', L.gutter - 6, sumTop + 22);

      // Rows
      const zoomed = span <= 120000;
      let prevLineage = null;
      for (const row of L.rows) {
        const si = row.smp.i;
        const y = row.y, hh = Math.max(1, row.h);
        ctx.fillStyle = intact;
        ctx.fillRect(L.gutter, y, L.plotW, hh);
        if (zoomed) {
          // Close in: every reference gene, so the gene structure shows
          ctx.fillStyle = css('--bg-subtle');
          for (const gi of genes) {
            const x0 = xOf(REF.start[gi]), x1 = xOf(REF.end[gi] + 1);
            ctx.fillRect(x0 + 0.5, y, Math.max(1, x1 - x0 - 1), hh);
          }
        }
        // A gene that is not there, or cut, fills the row; a pseudogene is a
        // band through its middle, more copies a band on top and a length
        // change one at the bottom. Absent genes are drawn last, over the rest.
        for (const pass of [['P', 'C', 'L'], ['M', 'T']]) {
          for (const gi of genes) {
            const c = BY_CODE[statusAt(si, gi)];
            if (!c || !enabled[c] || !pass.includes(c)) continue;
            const x0 = xOf(REF.start[gi]), x1 = xOf(REF.end[gi] + 1);
            const band = BAND[c];
            ctx.fillStyle = colors[c];
            ctx.fillRect(x0, y + hh * band[0], Math.max(1, x1 - x0 - (zoomed ? 1 : 0)), Math.max(1, hh * (band[1] - band[0])));
          }
        }
        // Lineage band and labels
        ctx.fillStyle = lineageColor(row.smp.lineage);
        ctx.fillRect(L.gutter - 6, y, 3, hh);
        if (L.labels) {
          ctx.fillStyle = si === view.highlight ? css('--accent') : ink;
          ctx.textAlign = 'right';
          ctx.fillText(trimText(ctx, row.smp.id, L.gutter - 44), L.gutter - 12, y + hh / 2 + 4);
          if (row.smp.lineage !== prevLineage) {
            ctx.fillStyle = lineageColor(row.smp.lineage);
            ctx.textAlign = 'left';
            ctx.font = '700 11px ' + css('--font');
            ctx.fillText(trimText(ctx, row.smp.lineage, 34), 2, y + hh / 2 + 4);
            ctx.font = '11px ' + css('--font');
          }
        } else if (row.smp.lineage !== prevLineage) {
          ctx.fillStyle = ink;
          ctx.textAlign = 'right';
          ctx.fillText(row.smp.lineage, L.gutter - 10, y + 9);
        }
        if (si === view.highlight) {
          ctx.strokeStyle = css('--accent');
          ctx.lineWidth = 1.5;
          ctx.strokeRect(L.gutter + 0.5, y - 1, L.plotW - 1, hh + 2);
        }
        prevLineage = row.smp.lineage;
      }
      if (!L.rows.length) {
        ctx.fillStyle = muted;
        ctx.textAlign = 'center';
        ctx.fillText('No genome matches the filters.', L.gutter + L.plotW / 2, L.top + 20);
      }
      where.textContent = (view.x1 - view.x0 >= REFLEN - 1 ? 'Whole reference, ' : '') + fmt(Math.round(view.x0)) + ' to ' + fmt(Math.round(view.x1)) + ' (' + kb(Math.round(view.x1 - view.x0)) + ', ' + plural(genes.length, 'gene') + ')';
      renderRegion(genes);
    }

    function trimText(ctx, text, w) {
      if (ctx.measureText(text).width <= w) return text;
      let t = text;
      while (t.length > 1 && ctx.measureText(t + '…').width > w) t = t.slice(0, -1);
      return t + '…';
    }

    function setView(a, b) {
      let w = Math.max(1500, b - a);
      w = Math.min(w, REFLEN);
      let x0 = a, x1 = a + w;
      if (x0 < 0) { x0 = 0; x1 = w; }
      if (x1 > REFLEN) { x1 = REFLEN; x0 = REFLEN - w; }
      view.x0 = x0;
      view.x1 = x1;
      draw();
    }

    function hit(ev) {
      const r = canvas.getBoundingClientRect();
      const x = ev.clientX - r.left, y = ev.clientY - r.top;
      if (!layout || x < layout.gutter || x > layout.gutter + layout.plotW) return null;
      const row = layout.rows.find((rw) => y >= rw.y - 1 && y <= rw.y + rw.h + 1);
      const pos = posOf(x);
      const tol = (view.x1 - view.x0) / layout.plotW * 2;
      const genes = inView(pos - tol, pos + tol);
      if (!genes.length) return { row, gi: -1, x, y };
      let best = genes[0];
      if (row) {
        const marked = genes.filter((gi) => {
          const c = BY_CODE[statusAt(row.smp.i, gi)];
          return c && enabled[c];
        });
        if (marked.length) best = marked[0];
        else best = genes.find((gi) => REF.start[gi] <= pos && REF.end[gi] >= pos) ?? genes[0];
      }
      return { row, gi: best, x, y, inStrip: y < layout.top };
    }

    let drag = null;
    canvas.addEventListener('pointerdown', (ev) => {
      if (ev.button !== 0) return;
      const r = canvas.getBoundingClientRect();
      const x = ev.clientX - r.left;
      if (!layout || x < layout.gutter) return;
      drag = { x0: x, x1: x, id: ev.pointerId };
      canvas.setPointerCapture(ev.pointerId);
    });
    canvas.addEventListener('pointermove', (ev) => {
      const r = canvas.getBoundingClientRect();
      if (drag) {
        drag.x1 = clamp(ev.clientX - r.left, layout.gutter, layout.gutter + layout.plotW);
        const a = Math.min(drag.x0, drag.x1), b = Math.max(drag.x0, drag.x1);
        if (b - a > 4) {
          brush.style.display = 'block';
          brush.style.left = a + 'px';
          brush.style.width = b - a + 'px';
        }
        tip.hide();
        return;
      }
      const hv = hit(ev);
      if (!hv || hv.gi < 0) { tip.hide(); return; }
      const gi = hv.gi;
      const lines = [h('div', null, h('b', null, geneLabel(gi)))];
      if (REF.product[gi]) lines.push(h('div', null, REF.product[gi]));
      lines.push(h('div', { class: 'muted' }, 'Reference ' + fmt(REF.start[gi]) + '-' + fmt(REF.end[gi])));
      if (hv.row && !hv.inStrip) {
        const c = BY_CODE[statusAt(hv.row.smp.i, gi)];
        lines.push(h('div', { style: { 'margin-top': '4px' } }, h('b', null, hv.row.smp.id), ' (' + hv.row.smp.lineage + '): ', stateChip(c), stateDetail(c, detailAt(hv.row.smp.i, gi))));
      } else {
        const counts = countStates(gi, layout.vis);
        lines.push(h('div', { style: { 'margin-top': '4px' } }, Object.keys(counts).filter((k) => counts[k]).map((k) => h('div', null, (k === 'I' ? 'Intact' : ST[k].label) + ': ' + counts[k]))));
      }
      tip.show(ev.clientX, ev.clientY, lines);
    });
    canvas.addEventListener('pointerleave', () => tip.hide());
    canvas.addEventListener('pointercancel', () => {
      drag = null;
      brush.style.display = 'none';
    });
    canvas.addEventListener('pointerup', (ev) => {
      if (!drag) return;
      const a = Math.min(drag.x0, drag.x1), b = Math.max(drag.x0, drag.x1);
      brush.style.display = 'none';
      drag = null;
      if (b - a > 4) setView(posOf(a), posOf(b));
      else {
        const hv = hit(ev);
        if (hv && hv.gi >= 0) lookupApi.show(hv.gi, hv.row ? hv.row.smp.i : -1);
      }
    });
    canvas.addEventListener('dblclick', () => zoomBy(2));
    canvas.addEventListener('keydown', (ev) => {
      if (ev.key === 'ArrowLeft') pan(-0.25);
      else if (ev.key === 'ArrowRight') pan(0.25);
      else if (ev.key === '+' || ev.key === '=') zoomBy(0.5);
      else if (ev.key === '-') zoomBy(2);
      else if (ev.key === '0') setView(0, REFLEN);
      else return;
      ev.preventDefault();
    });
    goto.addEventListener('keydown', (ev) => {
      if (ev.key !== 'Enter') return;
      const q = goto.value.trim();
      const m = /^([\d.,]+)\s*(mb|kb|bp)?\s*[-–:]\s*([\d.,]+)\s*(mb|kb|bp)?$/i.exec(q);
      if (m) {
        const unit = (u) => (/mb/i.test(u || '') ? 1e6 : /kb/i.test(u || '') ? 1e3 : 1);
        const u2 = m[4] || m[2];
        setView(parseFloat(m[1].replace(/,/g, '')) * unit(m[2] || u2), parseFloat(m[3].replace(/,/g, '')) * unit(u2));
        return;
      }
      const gi = findGene(q);
      if (gi >= 0 && REF.start[gi] != null) {
        setView(REF.start[gi] - 5000, REF.end[gi] + 5000);
        view.highlight = view.highlight;
      } else goto.setCustomValidity('No gene or region called ' + q);
      goto.reportValidity();
      setTimeout(() => goto.setCustomValidity(''), 1500);
    });

    function renderRegion(genes) {
      clear(regionBox);
      const span = view.x1 - view.x0;
      if (span > 400000) return;
      const vis = layout.vis;
      const rows = genes.map((gi) => Object.assign({ gi }, countStates(gi, vis)));
      regionBox.append(card('Genes in this region', 'The reference genes between ' + fmt(Math.round(view.x0)) + ' and ' + fmt(Math.round(view.x1)) + ', and how many of the ' + plural(vis.length, 'genome') + ' shown have each in each state. A click opens the gene.', null,
        table({
          columns: [
            { key: 'locus', label: 'Locus', value: (r) => REF.names[r.gi] },
            { key: 'name', label: 'Gene', value: (r) => REF.name[r.gi] || '' },
            { key: 'start', label: 'Start', num: true, value: (r) => REF.start[r.gi] },
            { key: 'end', label: 'End', num: true, value: (r) => REF.end[r.gi] },
            { key: 'product', label: 'Product', value: (r) => REF.product[r.gi] || '', wrap: true },
            { key: 'I', label: 'Intact', num: true },
            { key: 'M', label: 'Missing', num: true },
            { key: 'P', label: 'Pseudo', num: true },
            { key: 'T', label: 'Partial', num: true },
            { key: 'C', label: 'Copies', num: true },
            { key: 'L', label: 'Length', num: true },
          ],
          rows, unit: 'gene', csv: 'amap_region_' + Math.round(view.x0) + '_' + Math.round(view.x1) + '.csv', pageSize: 15,
          onRow: (r) => lookupApi.show(r.gi),
        })));
      regionBox.firstChild.style.marginTop = '16px';
    }

    function renderBlocks() {
      clear(blocksBox);
      const vis = new Set(visible().map((x) => x.name));
      const blocks = D.blocks.map((b) => {
        const runs = b.runs.filter((r) => vis.has(r[0]));
        return { b, runs };
      }).filter((x) => x.runs.length);
      const lineageOf = new Map(OK.map((x) => [x.name, x.lineage]));
      const total = new Map();
      for (const x of visible()) total.set(x.lineage, (total.get(x.lineage) || 0) + 1);
      const rows = blocks.map(({ b, runs }) => {
        const first = Math.min(...runs.map((r) => r[1])), last = Math.max(...runs.map((r) => r[2]));
        const byL = new Map();
        for (const name of new Set(runs.map((r) => r[0]))) byL.set(lineageOf.get(name), (byL.get(lineageOf.get(name)) || 0) + 1);
        const genes = new Set();
        for (const r of runs) for (let gi = r[1]; gi <= r[2]; gi++) genes.add(gi);
        return {
          first, last, genes: genes.size, span: REF.end[last] - REF.start[first] + 1,
          genomes: new Set(runs.map((r) => r[0])).size,
          lineages: [...byL.entries()].sort((a, b2) => LINEAGES.indexOf(a[0]) - LINEAGES.indexOf(b2[0])),
          whole: [...byL.entries()].filter(([l, n]) => n === total.get(l) && n > 1).map(([l]) => l),
        };
      });
      if (!rows.length) {
        blocksBox.append(h('div', { class: 'callout', style: { 'margin-top': '16px' } }, h('p', null, 'No genome shown lacks two or more consecutive reference genes.')));
        return;
      }
      const el = card('Blocks of missing genes', 'Two or more consecutive reference genes missing from a genome, grouped where genomes overlap: regions of difference and lineage deletions. A lineage in bold lacks the block in every one of its genomes. A click shows the block on the map.', null,
        table({
          columns: [
            { key: 'region', label: 'Region', value: (r) => REF.names[r.first] + ' to ' + REF.names[r.last] },
            { key: 'genes', label: 'Genes', num: true },
            { key: 'span', label: 'Span', num: true, render: (r) => kb(r.span) },
            { key: 'start', label: 'Reference start', num: true, value: (r) => REF.start[r.first] },
            { key: 'genomes', label: 'Genomes lacking it', num: true },
            { key: 'lineages', label: 'Lineages', value: (r) => r.lineages.map(([l, n]) => l + ' ' + n + '/' + total.get(l)).join(', '),
              render: (r) => r.lineages.map(([l, n]) => h('span', { class: 'lin', style: { '--c': lineageColor(l), 'margin-right': '10px', 'font-weight': r.whole.includes(l) ? '750' : '400' } }, h('span', { class: 'sw' }), l + ' ' + n + '/' + total.get(l))) },
          ],
          rows, unit: 'block', csv: 'amap_missing_blocks.csv', pageSize: 10, sort: 'genomes', dir: -1,
          onRow: (r) => { setView(REF.start[r.first] - 3000, REF.end[r.last] + 3000); fig.scrollIntoView({ behavior: 'smooth', block: 'center' }); },
        }));
      el.style.marginTop = '16px';
      blocksBox.append(el);
    }

    function exportMap() {
      const L = layout;
      const svg = s('svg', { viewBox: '0 0 ' + L.width + ' ' + L.height, width: L.width, height: L.height, 'font-family': css('--font'), 'font-size': '11' });
      svg.append(s('rect', { x: 0, y: 0, width: L.width, height: L.height, fill: css('--bg-elev') }));
      const genes = inView(view.x0, view.x1);
      for (const row of L.rows) {
        svg.append(s('rect', { x: L.gutter, y: row.y, width: L.plotW, height: Math.max(1, row.h), fill: css('--st-intact') }));
        svg.append(s('rect', { x: L.gutter - 6, y: row.y, width: 3, height: Math.max(1, row.h), fill: lineageColor(row.smp.lineage) }));
        if (L.labels) svg.append(s('text', { x: L.gutter - 12, y: row.y + row.h / 2 + 4, 'text-anchor': 'end', fill: css('--ink-soft') }, row.smp.id));
        for (const pass of [['P', 'C', 'L'], ['M', 'T']]) {
          for (const gi of genes) {
            const c = BY_CODE[statusAt(row.smp.i, gi)];
            if (!c || !enabled[c] || !pass.includes(c)) continue;
            const x0 = xOf(REF.start[gi]), x1 = xOf(REF.end[gi] + 1);
            const hh = Math.max(1, row.h), band = BAND[c];
            svg.append(s('rect', { x: x0.toFixed(2), y: (row.y + hh * band[0]).toFixed(2), width: Math.max(0.6, x1 - x0).toFixed(2), height: Math.max(0.6, hh * (band[1] - band[0])).toFixed(2), fill: colorOf(c) },
              s('title', null, row.smp.id + ': ' + REF.names[gi] + ', ' + ST[c].label)));
          }
        }
      }
      const step = niceStep((view.x1 - view.x0) / Math.max(2, L.plotW / 110));
      for (let t = Math.ceil(view.x0 / step) * step; t <= view.x1; t += step) {
        svg.append(s('text', { x: xOf(t), y: L.axisH - 9, 'text-anchor': 'middle', fill: css('--muted') }, axisLabel(t, step)));
      }
      download('amap_reference_map.svg', '<?xml version="1.0" encoding="UTF-8"?>\n' + new XMLSerializer().serializeToString(svg), 'image/svg+xml');
    }

    mapApi.highlight = (si) => {
      view.highlight = si;
      draw();
      const row = layout.rows.find((r) => r.smp.i === si);
      if (row) scroller.scrollTop = Math.max(0, row.y - 60);
    };
    mapApi.zoomTo = (gi) => setView(REF.start[gi] - 5000, REF.end[gi] + 5000);
    mapApi.redraw = () => draw();
    onFilter(() => { draw(); renderBlocks(); });
    new ResizeObserver(() => draw()).observe(fig);
    draw();
    renderBlocks();
  }

  function trimText2(ctx, text, w) {
    if (ctx.measureText(text).width <= w) return text;
    let s2 = text;
    while (s2.length > 1 && ctx.measureText(s2 + '…').width > w) s2 = s2.slice(0, -1);
    return s2 + '…';
  }

  function niceStep(x) {
    const p = Math.pow(10, Math.floor(Math.log10(Math.max(1, x))));
    for (const m of [1, 2, 2.5, 5, 10]) if (m * p >= x) return m * p;
    return 10 * p;
  }

  function axisLabel(t, step) {
    if (step >= 1e6) return (t / 1e6).toFixed(0) + ' Mb';
    if (step >= 1e5) return (t / 1e6).toFixed(1) + ' Mb';
    if (step >= 1e3) return (t / 1e3).toFixed(step >= 1e4 ? 0 : 1) + ' kb';
    return fmt(t);
  }

  function countStates(gi, vis) {
    const out = { I: 0, M: 0, P: 0, T: 0, C: 0, L: 0, S: 0 };
    for (const smp of vis) {
      const c = BY_CODE[statusAt(smp.i, gi)];
      out[c || 'I']++;
    }
    return out;
  }

  function stateDetail(letter, detail) {
    if (!letter || detail == null) return '';
    if (letter === 'P') return ': ' + (D.notes[+detail] || 'no reason given');
    if (letter === 'T') return ', cut at its ' + detail + '\' end';
    if (letter === 'C') return ', ' + detail + ' copies';
    if (letter === 'L') return ', ' + (+detail > 0 ? '+' : '') + detail + ' bp';
    if (letter === 'S') return ', in ' + plural(+detail, 'gene') + ' naming several reference genes';
    return '';
  }

  // ── 3. Gene lookup ─────────────────────────────────────────────────────

  const lookupApi = { show() {} };
  let searchIndex = null;

  function findGene(q) {
    const t = q.trim().toLowerCase();
    if (!t) return -1;
    for (let gi = 0; gi < NREF; gi++) if (REF.names[gi].toLowerCase() === t) return gi;
    for (let gi = 0; gi < NREF; gi++) if ((REF.name[gi] || '').toLowerCase() === t) return gi;
    const m = suggest(t, 1);
    return m.length ? m[0] : -1;
  }

  function suggest(q, limit) {
    if (!searchIndex) {
      searchIndex = REF.names.map((n, gi) => [n.toLowerCase(), (REF.name[gi] || '').toLowerCase(), (REF.product[gi] || '').toLowerCase()]);
    }
    const t = q.toLowerCase();
    const scored = [];
    for (let gi = 0; gi < NREF; gi++) {
      const [locus, name, product] = searchIndex[gi];
      let sc = 0;
      if (locus === t || name === t) sc = 100;
      else if (locus.startsWith(t) || name.startsWith(t)) sc = 60 - Math.min(20, locus.length - t.length);
      else if (locus.includes(t) || name.includes(t)) sc = 30;
      else if (t.length >= 3 && product.includes(t)) sc = 10;
      if (sc) scored.push([sc, gi]);
    }
    scored.sort((a, b) => b[0] - a[0] || a[1] - b[1]);
    return scored.slice(0, limit).map((x) => x[1]);
  }

  function renderLookupSection() {
    const body = clear(document.getElementById('lookup-body'));
    const input = h('input', { class: 'field', type: 'search', placeholder: 'katG, Rv1908c, catalase…', 'aria-label': 'Search a reference gene', autocomplete: 'off', spellcheck: 'false', role: 'combobox', 'aria-expanded': 'false', 'aria-controls': 'lookup-suggest' });
    const list = h('ul', { class: 'suggest', id: 'lookup-suggest', role: 'listbox', hidden: true });
    const out = h('div');
    let current = -1, focusSample = -1, active = -1, options = [];
    body.append(h('div', { class: 'lookup-search' }, input, list));

    // Examples: the drug-resistance genes Mycobrowser lists, or else the
    // genes whose state varies most between the genomes.
    const examples = [];
    if (D.meta.modes.mycobrowser) {
      for (const want of ['katG', 'rpoB', 'inhA', 'pncA', 'embB', 'gyrA', 'rpsL', 'ethA', 'Rv0678', 'pepQ']) {
        const gi = REF.name.findIndex((n, i) => n === want || REF.names[i] === want);
        if (gi >= 0) examples.push(gi);
      }
    }
    if (examples.length < 4) {
      const varied = [];
      for (let gi = 0; gi < NREF; gi++) {
        let n = 0;
        for (const smp of OK) if (statusAt(smp.i, gi)) n++;
        if (n && n < OK.length) varied.push([Math.min(n, OK.length - n), gi]);
      }
      varied.sort((a, b) => b[0] - a[0]);
      for (const [, gi] of varied.slice(0, 8 - examples.length)) examples.push(gi);
    }
    body.append(h('div', { class: 'examples' }, 'Try', examples.map((gi) => h('button', { class: 'mini-btn', type: 'button', onclick: () => show(gi) }, REF.name[gi] || REF.names[gi]))));
    body.append(out);

    function openList(items) {
      options = items;
      active = -1;
      clear(list);
      for (const gi of items) {
        list.append(h('li', { role: 'option', id: 'opt-' + gi, onmousedown: (e) => { e.preventDefault(); show(gi); } },
          h('b', null, REF.name[gi] || REF.names[gi]), h('span', null, (REF.name[gi] ? REF.names[gi] + ' · ' : '') + (REF.product[gi] || ''))));
      }
      list.hidden = !items.length;
      input.setAttribute('aria-expanded', items.length ? 'true' : 'false');
    }
    input.addEventListener('input', () => openList(input.value.trim() ? suggest(input.value.trim(), 12) : []));
    input.addEventListener('blur', () => setTimeout(() => openList([]), 120));
    input.addEventListener('keydown', (e) => {
      if (e.key === 'ArrowDown' || e.key === 'ArrowUp') {
        if (!options.length) return;
        active = clamp(active + (e.key === 'ArrowDown' ? 1 : -1), 0, options.length - 1);
        [...list.children].forEach((li, k) => li.setAttribute('aria-selected', k === active ? 'true' : 'false'));
        input.setAttribute('aria-activedescendant', 'opt-' + options[active]);
        e.preventDefault();
      } else if (e.key === 'Enter') {
        const gi = active >= 0 ? options[active] : findGene(input.value);
        if (gi >= 0) show(gi);
      } else if (e.key === 'Escape') openList([]);
    });

    function show(gi, si) {
      current = gi;
      focusSample = si == null ? -1 : si;
      input.value = REF.name[gi] || REF.names[gi];
      openList([]);
      draw();
      document.getElementById('lookup').scrollIntoView({ behavior: 'smooth' });
    }

    function draw() {
      clear(out);
      if (current < 0) return;
      const gi = current;
      const vis = visible();
      const refLen = REF.start[gi] != null ? REF.end[gi] - REF.start[gi] + 1 : null;
      const counts = countStates(gi, vis);
      const drugs = (REF.drugs && REF.drugs[gi]) || '';
      const facts = [
        ['In the reference', REF.start[gi] != null ? fmt(REF.start[gi]) + '-' + fmt(REF.end[gi]) + ' (' + (REF_PLUS[gi] ? '+' : '-') + ')' : 'no coordinates'],
        ['Length', (REF.aa[gi] ? fmt(REF.aa[gi]) + ' aa' : '') + (refLen ? ', ' + fmt(refLen) + ' bp' : '')],
        REF.cat && REF.cat[gi] ? ['Functional category', REF.cat[gi]] : null,
        REF.cog && REF.cog[gi] ? ['COG', REF.cog[gi]] : null,
        REF.refpseudo && REF.refpseudo[gi] >= 0 ? ['Pseudogene in the reference', REF.refpseudo[gi] ? 'yes (Mycobrowser)' : 'no (Mycobrowser)'] : null,
        drugs ? ['Drug resistance', h('span', null, drugs.split(',').map((d) => h('span', { class: 'drug', title: 'Mycobrowser lists resistance mutations in this gene for ' + d.trim() }, d.trim())))] : null,
      ].filter(Boolean);
      const stack = h('div', { class: 'stack', role: 'img', 'aria-label': 'States across the genomes shown' });
      const legend = h('div', { class: 'legend' });
      for (const k of ['I', 'M', 'P', 'T', 'C', 'L', 'S']) {
        if (!counts[k]) continue;
        const color = k === 'I' ? 'var(--amap)' : 'var(' + ST[k].color + ')';
        stack.append(h('span', { style: { flex: String(counts[k]), '--c': color }, title: (k === 'I' ? 'Intact' : ST[k].label) + ': ' + counts[k] }));
        legend.append(h('span', { class: 'st', style: { '--c': color } }, h('span', { class: 'sw' }), (k === 'I' ? 'Intact' : ST[k].label) + ' ' + counts[k]));
      }
      // One row per genome and gene naming the reference gene
      const rows = [];
      for (const smp of vis) {
        const c = BY_CODE[statusAt(smp.i, gi)];
        const det = detailAt(smp.i, gi);
        const list = refRows(smp.i).get(gi) || [];
        if (!list.length) rows.push({ smp, letter: c, det, gene: null });
        for (const row of list) rows.push({ smp, letter: c, det, gene: geneRow(smp.i, row) });
      }
      const product = REF.product[gi] || '';
      out.append(h('div', { class: 'card gene-card' },
        h('div', { class: 'gene-title' }, h('h3', null, REF.name[gi] || REF.names[gi]), REF.name[gi] ? h('span', { class: 'locus' }, REF.names[gi]) : null),
        product ? h('p', { style: { margin: '6px 0 0', color: 'var(--ink-soft)' } }, product, REF.src && REF.src[gi] === 'e' ? h('span', { class: 'muted' }, ' (from the genes\' own annotation)') : null) : null,
        h('div', { class: 'facts' }, facts.map(([k, v]) => h('div', null, h('span', { class: 'k' }, k), v))),
        stack, legend,
        h('div', { class: 'tools', style: { 'margin-top': '10px' } },
          REF.start[gi] != null ? h('button', { class: 'mini-btn', type: 'button', onclick: () => { mapApi.zoomTo(gi); location.hash = '#map'; } }, 'Show on the map') : null)));
      const tbl = table({
        columns: [
          { key: 'genome', label: 'Genome', value: (r) => r.smp.id, render: (r) => h('b', null, r.smp.id) },
          { key: 'lineage', label: 'Lineage', value: (r) => r.smp.lineage, render: (r) => lineageChip(r.smp.lineage) },
          { key: 'state', label: 'State', value: (r) => (r.letter ? ST[r.letter].label : 'Intact'), render: (r) => stateChip(r.letter) },
          { key: 'id', label: 'Locus tag', value: (r) => (r.gene ? r.gene.id : '') },
          { key: 'contig', label: 'Contig', value: (r) => (r.gene ? r.gene.contig : '') },
          { key: 'start', label: 'Start', num: true, value: (r) => (r.gene ? r.gene.start : null) },
          { key: 'end', label: 'End', num: true, value: (r) => (r.gene ? r.gene.end : null) },
          { key: 'strand', label: 'Strand', value: (r) => (r.gene ? r.gene.strand : '') },
          { key: 'aa', label: 'Length', num: true, value: (r) => (r.gene ? r.gene.length : null), render: (r) => (r.gene && r.gene.length != null ? fmt(r.gene.length) + ' bp' : '') },
          { key: 'delta', label: 'vs reference', num: true, value: (r) => (r.gene && refLen && r.gene.length != null ? r.gene.length - refLen : null), render: (r) => (r.gene && refLen && r.gene.length != null ? (r.gene.length - refLen > 0 ? '+' : '') + fmt(r.gene.length - refLen) + ' bp' : '') },
          { key: 'note', label: 'Note', wrap: true, value: (r) => (r.gene ? [r.gene.note, r.gene.partial ? 'cut at its ' + r.gene.partial.replace('\'', '') + '\' end by a contig end' : '', r.gene.across ? 'runs across the origin of its circular sequence' : '', r.gene.rvs.length > 1 ? 'also names ' + r.gene.rvs.filter((x) => x !== REF.names[gi]).join(', ') : ''].filter(Boolean).join('; ') : (r.letter === 'M' ? 'no gene names it' : '')) },
          { key: 'origin', label: 'How AMAP got it', value: (r) => (r.gene ? ORIGIN_LABEL[r.gene.origin] || '' : '') },
        ],
        rows, unit: 'row', csv: 'amap_gene_' + REF.names[gi] + '.csv', pageSize: 50, sort: 'lineage',
        rowClass: (r) => (r.smp.i === focusSample ? 'sel' : ''),
        onRow: (r) => openSample(r.smp.i),
      });
      tbl.style.marginTop = '14px';
      out.append(tbl);
    }

    lookupApi.show = show;
    onFilter(() => draw());
  }

  // ── 4. Correction ──────────────────────────────────────────────────────

  const correctionApi = { select() {} };

  function renderCorrectionSection() {
    const body = clear(document.getElementById('correction-body'));
    if (!D.meta.modes.cases) {
      body.append(callout('absent', 'No case tables.', 'The report was written without TYPE_CASES\' tables, so it cannot say which of Bakta\'s calls AMAP kept, replaced or removed.'));
      return;
    }
    let selected = -1; // -1: every genome shown, added up
    const pick = h('select', { class: 'field', 'aria-label': 'Genome' });
    const figBox = h('div', { class: 'fig' });
    const listBox = h('div');
    const tableBox = h('div');
    body.append(card('From Bakta\'s calls to AMAP\'s genes', 'Each bar starts where the one above ends: what the cases removed, then what the reference added. For one genome, a click on a bar lists its genes.',
      [pick, svgButton(() => chartBox.querySelector('svg'), 'amap_correction.svg')], figBox, listBox));
    body.append(tableBox);

    function fillPick() {
      clear(pick);
      pick.append(h('option', { value: '-1' }, 'Every genome shown, added up'));
      for (const smp of visible()) if (smp.flow) pick.append(h('option', { value: String(smp.i), selected: smp.i === selected }, smp.id + ' (' + smp.lineage + ')'));
      if (selected >= 0 && !visible().some((x) => x.i === selected)) selected = -1;
      pick.value = String(selected);
    }
    pick.addEventListener('change', () => { selected = +pick.value; draw(); });

    function flowOf() {
      const f = { c1: 0, c2: 0, c3: 0, c4: 0, k1: 0, k2: 0, k3: 0, k4: 0, R: 0, A: 0, other: 0, twice: 0 };
      const list = selected >= 0 ? [OK[selected]] : visible().filter((x) => x.flow);
      for (const smp of list) for (const k of Object.keys(f)) f[k] += smp.flow[k] || 0;
      return { f, n: list.length };
    }

    let current = null;
    const chartBox = h('div');
    const chart = responsive(chartBox, (w) => waterfall(w));

    // The waterfall at a given width: labels to the left of the bars, or
    // above them where the width is a phone's.
    function waterfall(W) {
      const { steps, total, final } = current;
      const narrow = W < 600;
      const labelW = narrow ? 0 : Math.min(300, Math.round(W * 0.4)), valW = 72, barH = 18;
      const rowH = narrow ? 44 : 30;
      const H = steps.length * rowH + 12;
      const max = Math.max(total, final, 1);
      const x = (v) => labelW + v / max * (W - labelW - valW);
      const svg = s('svg', { viewBox: '0 0 ' + W + ' ' + H, width: W, height: H, role: 'img', 'aria-label': 'Waterfall from Bakta\'s calls to AMAP\'s genes' });
      let level = 0, prevBar = 0;
      steps.forEach((st, k) => {
        const top = 6 + k * rowH;
        const y = narrow ? top + 18 : top;
        let a, b;
        if (st.kind === 'total') { a = 0; b = st.value; level = st.value; }
        else if (st.kind === 'minus') { a = level - st.value; b = level; level -= st.value; }
        else { a = level; b = level + st.value; level += st.value; }
        const g = s('g', { class: 'hit', tabindex: selected >= 0 ? '0' : null, role: selected >= 0 ? 'button' : null, 'aria-label': st.label + ': ' + fmt(st.value) });
        if (narrow) g.append(s('text', { x: 0, y: top + 12 }, st.label));
        else g.append(s('text', { x: labelW - 10, y: y + barH / 2 + 4, 'text-anchor': 'end' }, st.label));
        if (st.parts) {
          let x0 = a;
          st.parts.forEach(([pid, v, what], pk) => {
            if (!v) return;
            const r = s('rect', { class: 'mark', x: x(x0), y, width: Math.max(1, x(x0 + v) - x(x0) - (pk < st.parts.length - 1 ? 2 : 0)), height: barH, rx: 3,
              fill: 'var(' + (st.id === 'amap' && pid === 'ref' ? '--ref' : st.color) + ')', opacity: st.id === 'bakta' ? String(1 - pk * 0.18) : '1' });
            const pg = s('g', { class: 'hit', tabindex: selected >= 0 ? '0' : null, role: selected >= 0 ? 'button' : null }, r, s('title', null, what + ': ' + fmt(v)));
            pg.addEventListener('click', (ev) => { ev.stopPropagation(); listFor(st.id + ':' + pid, what); });
            pg.addEventListener('keydown', (ev) => { if (ev.key === 'Enter') listFor(st.id + ':' + pid, what); });
            g.append(pg);
            x0 += v;
          });
        } else {
          g.append(s('rect', { class: 'mark', x: x(a), y, width: Math.max(st.value ? 1 : 0, x(b) - x(a)), height: barH, rx: 3, fill: 'var(' + st.color + ')' }));
          if (k > 0) g.append(s('line', { x1: x(st.kind === 'minus' ? b : a), x2: x(st.kind === 'minus' ? b : a), y1: prevBar + barH, y2: y, stroke: 'var(--line-strong)', 'stroke-width': 1 }));
          g.append(s('title', null, st.label + ': ' + fmt(st.value) + (st.note ? '. ' + st.note : '')));
          g.addEventListener('click', () => listFor(st.id, st.label));
          g.addEventListener('keydown', (ev) => { if (ev.key === 'Enter') listFor(st.id, st.label); });
        }
        g.append(s('text', { x: x(b) + 8, y: y + barH / 2 + 4, class: 'num', fill: 'var(--ink)' }, (st.value && st.kind === 'minus' ? '−' : st.value && st.kind === 'plus' ? '+' : '') + fmt(st.value)));
        svg.append(g);
        prevBar = y;
      });
      return svg;
    }

    function draw() {
      const { f, n } = flowOf();
      clear(listBox);
      const total = f.c1 + f.c2 + f.c3 + f.c4;
      const kept = f.k1 + f.k2 + f.k3 + f.k4;
      const final = kept + f.R + f.A + f.other;
      const steps = [
        { id: 'bakta', label: 'Bakta\'s calls', kind: 'total', value: total, color: '--bakta',
          parts: [['2', f.c2, 'case 2, perfect match'], ['3', f.c3, 'case 3, imperfect match'], ['4', f.c4, 'case 4, conflict'], ['1', f.c1, 'case 1, no reference gene']] },
        { id: 'r1', label: 'Case 1 calls removed', kind: 'minus', value: f.c1 - f.k1, color: '--drop', note: 'no reference gene and filtered out: a known false positive or a shifted copy of a reference gene' },
        { id: 'r4', label: 'Case 4 calls replaced', kind: 'minus', value: f.c4 - f.k4, color: '--drop', note: 'the reference\'s coordinates won the conflict' },
        { id: 'r23', label: 'Case 2 and 3 calls dropped later', kind: 'minus', value: f.c2 - f.k2 + f.c3 - f.k3, color: '--drop', note: 'lost in the last deduplication or the stop codon check' },
        { id: 'R', label: 'Genes placed from the reference over a call', kind: 'plus', value: f.R, color: '--ref', note: 'the reference protein\'s coordinates, where Bakta called something that disagreed' },
        { id: 'A', label: 'Genes added where Bakta called nothing', kind: 'plus', value: f.A, color: '--ref', note: 'reference genes Bakta missed' },
        f.other ? { id: 'other', label: 'Other genes', kind: 'plus', value: f.other, color: '--ref' } : null,
        { id: 'amap', label: 'AMAP\'s genes', kind: 'total', value: final, color: '--amap',
          parts: [['k', kept, 'Bakta\'s calls kept'], ['ref', f.R + f.A + f.other, 'placed from the reference']] },
      ].filter(Boolean);
      current = { steps, total, final };
      clear(figBox).append(chartBox);
      chart.redraw();
      const sw = (color, opacity, label) => h('span', { class: 'st', style: { '--c': 'var(' + color + ')' } }, h('span', { class: 'sw', style: { opacity: String(opacity) } }), label);
      figBox.append(h('div', { class: 'legend', style: { margin: '4px 0 0' } },
        h('span', { class: 'muted' }, 'Bakta\'s calls:'), sw('--bakta', 1, 'case 2'), sw('--bakta', 0.82, 'case 3'), sw('--bakta', 0.64, 'case 4'), sw('--bakta', 0.46, 'case 1'),
        h('span', { class: 'muted', style: { 'margin-left': '8px' } }, 'AMAP\'s genes:'), sw('--amap', 1, 'Bakta\'s calls kept'), sw('--ref', 1, 'from the reference')));
      const parts = [];
      parts.push(n === 1 ? (selected >= 0 ? OK[selected].id : visible().find((x) => x.flow).id) + ': ' : plural(n, 'genome') + ' added up: ');
      parts.push(fmt(total) + ' Bakta calls became ' + fmt(final) + ' genes. ' + fmt(f.k2 + f.k3) + ' matched a reference gene (cases 2 and 3) and stand as called.');
      if (f.twice) parts.push(' ' + plural(f.twice, 'call') + ' sat in two case tables at once and ' + (f.twice === 1 ? 'is' : 'are') + ' counted once, under the better case.');
      figBox.append(h('p', { class: 'map-hint' }, parts.join('')));
      if (selected < 0) listBox.append(h('p', { class: 'map-hint' }, 'Pick one genome above to list the genes behind each bar.'));
      fillTable();
    }

    function listFor(id, what) {
      if (selected < 0) return;
      const si = selected;
      const g = genesOf(si);
      const gone = goneOf(si);
      clear(listBox);
      let rows = [], cols = [];
      const finalCols = [
        { key: 'id', label: 'Locus tag' },
        { key: 'rvs', label: 'Reference gene', value: (r) => r.rvs.join(', '), wrap: true },
        { key: 'contig', label: 'Contig' },
        { key: 'start', label: 'Start', num: true }, { key: 'end', label: 'End', num: true }, { key: 'strand', label: 'Strand' },
        { key: 'origin', label: 'How AMAP got it', value: (r) => ORIGIN_LABEL[r.origin] || r.origin },
      ];
      const pickOrigin = (set) => {
        const out = [];
        for (let row = 0; row < g.n; row++) if (set.includes(g.origin[row])) out.push(geneRow(si, row));
        return out;
      };
      const goneCols = [
        { key: 'tag', label: 'Bakta\'s locus tag' },
        { key: 'case', label: 'Case' },
        { key: 'call', label: 'Bakta\'s call', value: (r) => r.contig + ':' + fmt(r.start) + '-' + fmt(r.end) + ' ' + r.strand },
        { key: 'by', label: 'AMAP\'s gene in its place', value: (r) => (r.row >= 0 ? g.ids[r.row] + ' (' + rvsOf(g, r.row).join(', ') + ')' : 'none') },
        { key: 'amap', label: 'AMAP\'s coordinates', value: (r) => (r.row >= 0 ? fmt(g.start[r.row]) + '-' + fmt(g.end[r.row]) : '') },
        { key: 'ds', label: 'Start moved', num: true, value: (r) => (r.row >= 0 ? (g.plus[r.row] ? g.start[r.row] - r.start : r.end - g.end[r.row]) : null),
          render: (r) => (r.row >= 0 ? signed(g.plus[r.row] ? g.start[r.row] - r.start : r.end - g.end[r.row]) : ''), title: 'How far upstream (+) or downstream (-) of Bakta\'s start codon AMAP starts the gene' },
      ];
      if (id.startsWith('bakta:')) {
        const c = id.split(':')[1];
        rows = pickOrigin([c]).concat((gone || []).filter((x) => x.case === c).map((x) => ({ id: x.tag, rvs: [], contig: x.contig, start: x.start, end: x.end, strand: x.strand, origin: 'gone' })));
        cols = finalCols.slice(0, 6).concat([{ key: 'kept', label: 'In the final annotation', value: (r) => (r.origin === 'gone' ? 'no' : 'yes') }]);
      } else if (id === 'r1' || id === 'r4' || id === 'r23') {
        if (!gone) {
          listBox.append(callout('absent', 'Not embedded.', 'This run has more genomes than the report lists calls for; ' + OK[si].name + '\'s case tables are in CDS_CASES/' + OK[si].name + '/.'));
          return;
        }
        const want = id === 'r1' ? ['1'] : id === 'r4' ? ['4'] : ['2', '3'];
        rows = gone.filter((x) => want.includes(x.case));
        cols = goneCols;
      } else if (id === 'R' || id === 'A' || id === 'other') {
        rows = pickOrigin([id === 'other' ? '?' : id]);
        cols = finalCols;
      } else if (id === 'amap:k') {
        rows = pickOrigin(['1', '2', '3', '4']);
        cols = finalCols;
      } else if (id === 'amap:ref') {
        rows = pickOrigin(['R', 'A', '?']);
        cols = finalCols;
      } else return;
      const el = card(what + ' in ' + OK[si].id, plural(rows.length, 'gene') + '.', null,
        table({ columns: cols, rows, unit: 'gene', csv: OK[si].name + '_' + id.replace(':', '_') + '.csv', pageSize: 15 }));
      el.style.marginTop = '14px';
      listBox.append(el);
    }

    function fillTable() {
      clear(tableBox);
      const rows = visible().filter((x) => x.flow);
      const el = card('Every genome', 'The same steps, one genome per row.', null, table({
        columns: [
          { key: 'genome', label: 'Genome', value: (r) => r.id, render: (r) => h('b', null, r.id) },
          { key: 'lineage', label: 'Lineage', value: (r) => r.lineage, render: (r) => lineageChip(r.lineage) },
          { key: 'calls', label: 'Bakta calls', num: true, value: (r) => r.flow.c1 + r.flow.c2 + r.flow.c3 + r.flow.c4 },
          { key: 'c2', label: 'Case 2', num: true, value: (r) => r.flow.c2 },
          { key: 'c3', label: 'Case 3', num: true, value: (r) => r.flow.c3 },
          { key: 'c4', label: 'Case 4', num: true, value: (r) => r.flow.c4 },
          { key: 'k4', label: 'Case 4 kept', num: true, value: (r) => r.flow.k4 },
          { key: 'c1', label: 'Case 1', num: true, value: (r) => r.flow.c1 },
          { key: 'k1', label: 'Case 1 kept', num: true, value: (r) => r.flow.k1 },
          { key: 'R', label: 'Over a call', num: true, value: (r) => r.flow.R, title: 'Genes placed from the reference over a Bakta call' },
          { key: 'A', label: 'Added', num: true, value: (r) => r.flow.A, title: 'Genes added where Bakta called nothing' },
          { key: 'genes', label: 'AMAP genes', num: true },
        ],
        rows, unit: 'genome', csv: 'amap_correction.csv', pageSize: 25, sort: 'lineage',
        onRow: (r) => { selected = r.i; fillPick(); draw(); body.scrollIntoView({ behavior: 'smooth' }); },
      }));
      el.style.marginTop = '16px';
      tableBox.append(el);
    }

    correctionApi.select = (si) => { selected = si; fillPick(); draw(); };
    onFilter(() => { fillPick(); draw(); });
    fillPick();
    draw();
  }

  const signed = (v) => (v > 0 ? '+' : '') + fmt(v) + ' bp';
  const baseName = (path) => String(path || '').split(/[\\/]/).pop();

  // ── 5. Pseudogenes ─────────────────────────────────────────────────────

  function renderPseudoSection() {
    const body = clear(document.getElementById('pseudo-body'));
    const box = h('div');
    body.append(box);

    function draw() {
      clear(box);
      const vis = visible();
      const code = ST.P.code;
      const byLineage = new Map();
      for (const smp of vis) {
        if (!byLineage.has(smp.lineage)) byLineage.set(smp.lineage, []);
        byLineage.get(smp.lineage).push(smp);
      }
      // Each gene's pattern across the genomes shown: specific to some
      // lineages (a pseudogene in all their genomes and in no other), in every
      // genome, in all but one or two, or scattered.
      const rows = [];
      for (let gi = 0; gi < NREF; gi++) {
        let n = 0;
        const reasons = new Map();
        for (const smp of vis) {
          if (statusAt(smp.i, gi) === code) {
            n++;
            const r = D.notes[+detailAt(smp.i, gi)] || '';
            reasons.set(r, (reasons.get(r) || 0) + 1);
          }
        }
        if (!n) continue;
        const whole = [], partly = [];
        for (const [l, list] of byLineage) {
          const k = list.filter((smp) => statusAt(smp.i, gi) === code).length;
          if (k === list.length) whole.push(l);
          else if (k) partly.push(l);
        }
        let kind, label;
        if (n === vis.length) {
          kind = 1;
          label = 'all';
        } else if (!partly.length && whole.some((l) => byLineage.get(l).length >= 2)) {
          kind = 0;
          label = whole.join(', ') + ' only';
        } else if (vis.length - n <= 2 && vis.length >= 5) {
          kind = 2;
          label = 'all but ' + vis.filter((smp) => statusAt(smp.i, gi) !== code).map((x) => x.id).join(', ');
        } else {
          kind = 3;
          label = n + '/' + vis.length;
        }
        rows.push({
          gi, n, kind, label, lineages: kind === 0 ? whole : [],
          reason: [...reasons.entries()].sort((a, b) => b[1] - a[1])[0][0],
          reasons: [...reasons.entries()].sort((a, b) => b[1] - a[1]),
          refp: REF.refpseudo ? REF.refpseudo[gi] : -1,
        });
      }
      rows.sort((a, b) => a.kind - b.kind || (a.kind === 0 ? a.lineages.length - b.lineages.length : 0) || b.n - a.n || a.gi - b.gi);
      if (!rows.length) {
        box.append(callout('', 'No pseudogenes.', 'No genome shown has a gene AMAP flags as a pseudogene.'));
        return;
      }
      const nSpecific = rows.filter((r) => r.kind === 0).length;
      const nEvery = rows.filter((r) => r.kind === 1).length;
      const nAllBut = rows.filter((r) => r.kind === 2).length;
      const lines = [plural(rows.length, 'gene') + ' ' + (rows.length === 1 ? 'is' : 'are') + ' a pseudogene in at least one of the ' + plural(vis.length, 'genome') + ' shown: ' +
        fmt(nSpecific) + ' in every genome of some lineages and in no other, ' + fmt(nEvery) + ' in every genome' +
        (nAllBut ? ', ' + fmt(nAllBut) + ' in all but one or two' : '') + '.'];
      if (D.meta.modes.mycobrowser) {
        const refP = [];
        for (let gi = 0; gi < NREF; gi++) if (REF.refpseudo[gi] === 1) refP.push(gi);
        const flaggedEverywhere = refP.filter((gi) => rows.some((r) => r.gi === gi && r.kind === 1)).length;
        const notRef = rows.filter((r) => r.kind === 1 && r.refp === 0).length;
        lines.push(' Mycobrowser calls ' + plural(refP.length, 'reference gene') + ' pseudogenes; AMAP flags ' + fmt(flaggedEverywhere) + ' of them in every genome. ' +
          fmt(notRef) + ' of the genes flagged in every genome are not pseudogenes in Mycobrowser: broken in all these genomes but not in the reference, or flagged by a rule the reference does not meet.');
      }
      box.append(h('p', { class: 'muted', style: { 'max-width': '60em' } }, lines.join('')));

      // The matrix: genes against genomes
      const fig = h('div', { class: 'canvas-fig' });
      const canvas = h('canvas', { 'aria-label': 'Pseudogenes by gene and genome' });
      const scroll = h('div', { class: 'canvas-scroll' }, canvas);
      fig.append(scroll);
      const legend = h('div', { class: 'legend', style: { margin: '0 0 10px' } },
        stateChip('P'), stateChip('M'), stateChip('T'), h('span', { class: 'st', style: { '--c': 'var(--st-intact)' } }, h('span', { class: 'sw' }), 'Intact or other'),
        h('span', { class: 'muted' }, '▸ specific to the lineages named on the right'));
      box.append(card('Pseudogenes by genome', 'Rows are genes, columns genomes grouped by lineage. First the genes broken in every genome of some lineages and in no other, then those broken in every genome, in all but one or two, and the rest by how many genomes have them. Hover for the reason, click for the gene.', [svgButton(() => matrixSvg(), 'amap_pseudogenes.svg')], legend, fig));
      let lay = null;
      function layout() {
        const width = Math.max(320, fig.clientWidth || 800);
        const narrow = width < 640;
        const labelW = narrow ? 110 : 190, sideW = narrow ? 40 : 120;
        const n = vis.length;
        const plotW = width - labelW - sideW - 12;
        const cellW = Math.min(26, plotW / Math.max(1, n + (byLineage.size - 1) * 0.5));
        const rowH = 15, top = 30;
        const xs = [];
        let x = labelW, prev = null;
        for (const smp of vis) {
          if (prev && smp.lineage !== prev) x += cellW * 0.5;
          xs.push(x);
          x += cellW;
          prev = smp.lineage;
        }
        return { width, labelW, sideW, cellW, rowH, top, xs, end: x, height: top + rows.length * rowH + 8 };
      }
      function paint() {
        lay = layout();
        const L = lay;
        const ctx = sizeCanvas(canvas, L.width, L.height);
        ctx.font = '11px ' + css('--font');
        const colors = { M: colorOf('M'), P: colorOf('P'), T: colorOf('T') };
        const intact = css('--st-intact');
        // Lineage headers
        let prev = null;
        vis.forEach((smp, k) => {
          if (smp.lineage !== prev) {
            ctx.fillStyle = lineageColor(smp.lineage);
            const count = vis.filter((x) => x.lineage === smp.lineage).length;
            ctx.fillRect(L.xs[k], L.top - 8, count * L.cellW - 1, 3);
            ctx.fillStyle = css('--ink-soft');
            ctx.textAlign = 'left';
            if (count * L.cellW > 16) ctx.fillText(smp.lineage, L.xs[k], L.top - 13);
          }
          prev = smp.lineage;
        });
        rows.forEach((r, k) => {
          const y = L.top + k * L.rowH;
          ctx.fillStyle = r.kind === 0 ? css('--ink') : css('--ink-soft');
          ctx.textAlign = 'right';
          ctx.font = (r.kind === 0 ? '600 ' : '') + '11px ' + css('--font');
          ctx.fillText((r.kind === 0 ? '▸ ' : '') + (REF.name[r.gi] || REF.names[r.gi]), L.labelW - 8, y + 11);
          vis.forEach((smp, c) => {
            const letter = BY_CODE[statusAt(smp.i, r.gi)];
            ctx.fillStyle = letter && colors[letter] ? colors[letter] : intact;
            ctx.fillRect(L.xs[c], y + 1, Math.max(1, L.cellW - (L.cellW > 4 ? 2 : 0)), L.rowH - 2);
          });
          ctx.textAlign = 'left';
          ctx.font = '11px ' + css('--font');
          ctx.fillStyle = css('--muted');
          ctx.fillText(trimText2(ctx, r.label + (r.refp === 1 ? ' · ref' : ''), L.sideW - 4), L.end + 8, y + 11);
        });
      }
      canvas.addEventListener('pointermove', (ev) => {
        if (!lay) return;
        const rc = canvas.getBoundingClientRect();
        const x = ev.clientX - rc.left, y = ev.clientY - rc.top;
        const k = Math.floor((y - lay.top) / lay.rowH);
        if (k < 0 || k >= rows.length) { tip.hide(); return; }
        const r = rows[k];
        const c = lay.xs.findIndex((x0) => x >= x0 && x < x0 + lay.cellW);
        const lines = [h('div', null, h('b', null, geneLabel(r.gi))), REF.product[r.gi] ? h('div', null, REF.product[r.gi]) : null];
        if (c >= 0) {
          const smp = vis[c];
          const letter = BY_CODE[statusAt(smp.i, r.gi)];
          lines.push(h('div', { style: { 'margin-top': '4px' } }, h('b', null, smp.id), ' (' + smp.lineage + '): ', stateChip(letter), stateDetail(letter, detailAt(smp.i, r.gi))));
        } else {
          lines.push(h('div', null, 'Pseudogene in ' + r.n + ' of ' + vis.length + ' (' + r.label + ')'));
          for (const [reason, n] of r.reasons.slice(0, 3)) lines.push(h('div', { class: 'muted' }, n + ' × ' + (reason || 'no reason given')));
        }
        tip.show(ev.clientX, ev.clientY, lines);
      });
      canvas.addEventListener('pointerleave', () => tip.hide());
      canvas.addEventListener('click', (ev) => {
        const rc = canvas.getBoundingClientRect();
        const k = Math.floor((ev.clientY - rc.top - lay.top) / lay.rowH);
        if (k >= 0 && k < rows.length) lookupApi.show(rows[k].gi);
      });
      function matrixSvg() {
        const L = lay;
        const svg = s('svg', { viewBox: '0 0 ' + L.width + ' ' + L.height, width: L.width, height: L.height });
        const colors = { M: colorOf('M'), P: colorOf('P'), T: colorOf('T') };
        rows.forEach((r, k) => {
          const y = L.top + k * L.rowH;
          svg.append(s('text', { x: L.labelW - 8, y: y + 11, 'text-anchor': 'end', 'font-size': '11', fill: css('--ink-soft') }, REF.name[r.gi] || REF.names[r.gi]));
          vis.forEach((smp, c) => {
            const letter = BY_CODE[statusAt(smp.i, r.gi)];
            svg.append(s('rect', { x: L.xs[c].toFixed(2), y: y + 1, width: Math.max(1, L.cellW - 2).toFixed(2), height: L.rowH - 2, fill: letter && colors[letter] ? colors[letter] : css('--st-intact') }, s('title', null, smp.id + ': ' + (letter ? ST[letter].label : 'intact'))));
          });
        });
        vis.forEach((smp, c) => svg.append(s('text', { x: (L.xs[c] + L.cellW / 2).toFixed(1), y: L.top - 6, 'text-anchor': 'middle', 'font-size': '9', fill: lineageColor(smp.lineage) }, '|')));
        return svg;
      }
      new ResizeObserver(() => paint()).observe(fig);
      paint();

      const tbl = table({
        columns: [
          { key: 'locus', label: 'Locus', value: (r) => REF.names[r.gi] },
          { key: 'name', label: 'Gene', value: (r) => REF.name[r.gi] || '' },
          { key: 'product', label: 'Product', value: (r) => REF.product[r.gi] || '', wrap: true },
          { key: 'n', label: 'Genomes', num: true, render: (r) => fmt(r.n) + ' of ' + fmt(vis.length) },
          { key: 'label', label: 'Pattern', value: (r) => r.label, wrap: true },
          { key: 'refp', label: 'Reference pseudogene', value: (r) => (r.refp === 1 ? 'yes' : r.refp === 0 ? 'no' : 'not known') },
          { key: 'reason', label: 'Most common reason', wrap: true },
        ],
        rows, unit: 'gene', csv: 'amap_pseudogenes.csv', pageSize: 20,
        onRow: (r) => lookupApi.show(r.gi),
      });
      tbl.style.marginTop = '16px';
      box.append(tbl);
    }
    onFilter(draw);
    draw();
  }

  // ── 6. Beyond the reference ────────────────────────────────────────────

  function renderBeyondSection() {
    const body = clear(document.getElementById('beyond-body'));
    const box = h('div');
    body.append(box);

    function draw() {
      clear(box);
      const vis = visible();
      const visNames = new Set(vis.map((x) => x.name));
      const lineageOf = new Map(OK.map((x) => [x.name, x.lineage]));

      // Reference genes across the genomes shown
      const n = vis.length;
      const hist = new Array(n + 1).fill(0);
      const accessory = [];
      for (let gi = 0; gi < NREF; gi++) {
        let present = 0;
        const absentIn = [];
        for (const smp of vis) {
          if (statusAt(smp.i, gi) !== ST.M.code) present++;
          else absentIn.push(smp);
        }
        hist[present]++;
        if (present > 0 && present < n) accessory.push({ gi, present, absentIn });
      }
      const spectrum = (W) => {
        const H = 220, padL = 52, padB = 34, padT = 14;
        const maxH = Math.max(1, ...hist);
        const bw = (W - padL - 10) / (n + 1);
        const ylog = (v) => (v <= 0 ? 0 : Math.log10(v + 1) / Math.log10(maxH + 1));
        const svg = s('svg', { viewBox: '0 0 ' + W + ' ' + H, width: W, height: H, role: 'img', 'aria-label': 'How many genomes carry each reference gene' });
        for (const tk of [1, 10, 100, 1000, 10000].filter((v) => v <= maxH * 1.2)) {
          const y = H - padB - ylog(tk) * (H - padB - padT);
          svg.append(s('g', { class: 'grid' }, s('line', { x1: padL, x2: W - 10, y1: y, y2: y })), s('text', { x: padL - 6, y: y + 4, 'text-anchor': 'end', class: 'axis' }, fmt(tk)));
        }
        const every = Math.max(1, Math.ceil((n + 1) / Math.max(4, Math.floor((W - padL) / 28))));
        hist.forEach((v, k) => {
          const x = padL + k * bw;
          const hh = ylog(v) * (H - padB - padT);
          const color = k === n ? '--amap' : k === 0 ? '--st-missing' : '--ref';
          svg.append(s('g', { class: 'hit' },
            s('rect', { class: 'mark', x: x + 1, y: H - padB - hh, width: Math.max(1, bw - 2), height: hh, rx: 2, fill: 'var(' + color + ')' }),
            s('title', null, plural(v, 'reference gene') + ' in ' + k + ' of ' + n + ' genomes')));
          if (k % every === 0 || k === n) svg.append(s('text', { x: x + bw / 2, y: H - padB + 14, 'text-anchor': 'middle', class: 'axis' }, String(k)));
        });
        svg.append(s('text', { x: padL + (W - padL) / 2, y: H - 4, 'text-anchor': 'middle', class: 'axis' }, 'genomes carrying the gene (genes on a log scale)'));
        return svg;
      };
      const fig = h('div', { class: 'fig' });
      responsive(fig, spectrum);
      box.append(card('Reference genes across the genomes', fmt(hist[n]) + ' reference genes are in all ' + plural(n, 'genome') + ' shown, ' + fmt(accessory.length) + ' in some, ' + fmt(hist[0]) + ' in none. A gene counts as carried unless it is missing; a pseudogene is carried.',
        [svgButton(() => fig.querySelector('svg'), 'amap_reference_spectrum.svg')], fig,
        accessory.length ? table({
          columns: [
            { key: 'locus', label: 'Locus', value: (r) => REF.names[r.gi] },
            { key: 'name', label: 'Gene', value: (r) => REF.name[r.gi] || '' },
            { key: 'product', label: 'Product', value: (r) => REF.product[r.gi] || '', wrap: true },
            { key: 'present', label: 'Genomes with it', num: true, render: (r) => fmt(r.present) + ' of ' + fmt(n) },
            { key: 'absent', label: 'Missing from', wrap: true, value: (r) => (r.absentIn.length <= 6 ? r.absentIn.map((x) => x.id).join(', ') : summariseLineages(r.absentIn, vis)) },
          ],
          rows: accessory, unit: 'gene', csv: 'amap_reference_accessory.csv', pageSize: 10, sort: 'present', dir: 1,
          onRow: (r) => lookupApi.show(r.gi),
        }) : null));

      // Genes not in the reference
      const groups = D.noref.map(([product, byName]) => {
        let total = 0;
        const names = [];
        for (const [name, c] of Object.entries(byName)) if (visNames.has(name)) { total += c; names.push(name); }
        return { product, total, genomes: names.length, lineages: [...new Set(names.map((x) => lineageOf.get(x)))].sort((a, b) => LINEAGES.indexOf(a) - LINEAGES.indexOf(b)) };
      }).filter((x) => x.total);
      const perGenome = vis.map((x) => x.noref);
      const el2 = card('Genes not in the reference', 'Genes that name no reference gene: Bakta calls with nothing in the reference to compare (case 1 calls the filter kept). ' +
        (groups.length ? 'Grouped by the product Bakta or eggNOG gave them, so the same element in several genomes is one row; "hypothetical protein" groups unrelated genes.' : ''),
        null,
        h('p', { class: 'muted' }, 'Per genome: median ' + fmt(median(perGenome)) + ', from ' + fmt(Math.min(...perGenome)) + ' to ' + fmt(Math.max(...perGenome)) + '. Each genome\'s own list is in its row of section 1.'),
        groups.length ? table({
          columns: [
            { key: 'product', label: 'Product', wrap: true },
            { key: 'total', label: 'Genes', num: true },
            { key: 'genomes', label: 'Genomes', num: true },
            { key: 'lineages', label: 'Lineages', value: (r) => r.lineages.join(', '), render: (r) => r.lineages.map((l) => h('span', { style: { 'margin-right': '8px' } }, lineageChip(l))) },
          ],
          rows: groups, unit: 'product', csv: 'amap_not_in_reference.csv', pageSize: 10, sort: 'genomes', dir: -1,
        }) : callout('absent', '', 'No genome shown has a gene outside the reference.'));
      el2.style.marginTop = '16px';
      box.append(el2);

      // Copies
      const copies = [];
      for (let gi = 0; gi < NREF; gi++) {
        const list = [];
        for (const smp of vis) if (statusAt(smp.i, gi) === ST.C.code) list.push([smp, +detailAt(smp.i, gi)]);
        if (list.length) copies.push({ gi, list, max: Math.max(...list.map((x) => x[1])) });
      }
      const el3 = card('Genes found more than once', 'Reference genes with more than one copy in a genome: AMAP keeps a further copy when it is at least 90% as long as the first and clear of the copies kept. The copy number is what differs between genomes.', null,
        copies.length ? table({
          columns: [
            { key: 'locus', label: 'Locus', value: (r) => REF.names[r.gi] },
            { key: 'name', label: 'Gene', value: (r) => REF.name[r.gi] || '' },
            { key: 'product', label: 'Product', value: (r) => REF.product[r.gi] || '', wrap: true },
            { key: 'genomes', label: 'Genomes with copies', num: true, value: (r) => r.list.length },
            { key: 'max', label: 'Most copies', num: true },
            { key: 'which', label: 'Copies by genome', wrap: true, value: (r) => r.list.map(([smp, k]) => smp.id + ' ×' + k).join(', ') },
          ],
          rows: copies, unit: 'gene', csv: 'amap_copies.csv', pageSize: 10, sort: 'genomes', dir: -1,
          onRow: (r) => lookupApi.show(r.gi),
        }) : callout('absent', '', 'No reference gene has more than one copy in the genomes shown.'));
      el3.style.marginTop = '16px';
      box.append(el3);

      // Pangenome
      let el4;
      if (D.pangene) {
        const rows = vis.filter((x) => D.pangene[x.name]).map((smp) => {
          const p = D.pangene[smp.name];
          let total = 0, all = 0, none = 0;
          for (const [k, v] of Object.entries(p.hist)) {
            total += v;
            if (+k === p.others) all += v;
            if (+k === 0) none += v;
          }
          return { smp, total, all, none, some: total - all - none, others: p.others };
        });
        el4 = card('Pangenome', 'What PANGENE found: each genome\'s genes aligned to every other genome. A gene in all of them is core to this run; one in none is this genome\'s own.', null,
          table({
            columns: [
              { key: 'genome', label: 'Genome', value: (r) => r.smp.id },
              { key: 'lineage', label: 'Lineage', value: (r) => r.smp.lineage, render: (r) => lineageChip(r.smp.lineage) },
              { key: 'total', label: 'Genes', num: true },
              { key: 'all', label: 'In every other genome', num: true },
              { key: 'some', label: 'In some', num: true },
              { key: 'none', label: 'In no other genome', num: true },
            ],
            rows, unit: 'genome', csv: 'amap_pangenome.csv', pageSize: 25,
          }));
      } else {
        el4 = card('Pangenome', null, null, callout('absent', 'Not in this run.', D.meta.modes.skip_functional ?
          'PANGENE runs with the functional annotation, and this run used --skip_functional_annotation.' :
          'The report was written without PANGENE\'s matrix.'));
      }
      el4.style.marginTop = '16px';
      box.append(el4);

      // Duplication communities
      let el5;
      if (D.communities && D.communities.length) {
        const keys = Object.keys(D.communities[0]);
        el5 = card('Homologous communities', 'Duplicated genes across genomes, grouped into communities by GET_ALLDUP (DUPLICATIONS/homologous_communities.tsv).', null,
          table({ columns: keys.map((k) => ({ key: k, label: k.replace(/_/g, ' '), num: D.communities.every((r) => r[k] === '' || !Number.isNaN(+r[k])), value: (r) => (r[k] !== '' && !Number.isNaN(+r[k]) ? +r[k] : r[k]), wrap: /node|gene|member/i.test(k) })), rows: D.communities, unit: 'community', csv: 'amap_homologous_communities.csv', pageSize: 10 }));
      } else {
        el5 = card('Homologous communities', null, null, callout('absent', 'Not in this run.', D.meta.modes.skip_functional ?
          'The duplication analysis runs with the functional annotation, and this run used --skip_functional_annotation. The copies above come from the annotation itself.' :
          'GET_ALLDUP found no duplicated genes shared between genomes, or the report was written without its table.'));
      }
      el5.style.marginTop = '16px';
      box.append(el5);
    }
    onFilter(draw);
    draw();
  }

  function summariseLineages(list, vis) {
    const total = new Map();
    for (const x of vis) total.set(x.lineage, (total.get(x.lineage) || 0) + 1);
    const by = new Map();
    for (const x of list) by.set(x.lineage, (by.get(x.lineage) || 0) + 1);
    return [...by.entries()].map(([l, n]) => l + ' ' + n + '/' + total.get(l)).join(', ');
  }

  // ── 7. Functional profile ──────────────────────────────────────────────

  const COG_NAMES = {
    A: 'RNA processing and modification', B: 'Chromatin structure and dynamics', C: 'Energy production and conversion',
    D: 'Cell cycle control, cell division', E: 'Amino acid transport and metabolism', F: 'Nucleotide transport and metabolism',
    G: 'Carbohydrate transport and metabolism', H: 'Coenzyme transport and metabolism', I: 'Lipid transport and metabolism',
    J: 'Translation, ribosomal structure', K: 'Transcription', L: 'Replication, recombination and repair',
    M: 'Cell wall, membrane, envelope biogenesis', N: 'Cell motility', O: 'Post-translational modification, chaperones',
    P: 'Inorganic ion transport', Q: 'Secondary metabolites biosynthesis', R: 'General function prediction only',
    S: 'Function unknown', T: 'Signal transduction mechanisms', U: 'Intracellular trafficking, secretion',
    V: 'Defense mechanisms', W: 'Extracellular structures', X: 'Mobilome: prophages, transposons', Z: 'Cytoskeleton', '-': 'No COG category',
  };

  function renderFunctionSection() {
    const body = clear(document.getElementById('function-body'));
    const box = h('div');
    body.append(box);

    // Categories by lineage, the categories that differ most first: a dot per
    // lineage on one line per category, at the lineage's mean per genome.
    function dotPlot(title, desc, categories, valueOf, label, fileName) {
      const vis = visible();
      const lins = LINEAGES.filter((l) => vis.some((x) => x.lineage === l));
      const rows = categories.map((cat) => {
        const means = lins.map((l) => {
          const v = vis.filter((x) => x.lineage === l).map((x) => valueOf(x, cat)).filter((x) => x != null);
          return v.length ? v.reduce((a, b) => a + b, 0) / v.length : null;
        });
        const vals = means.filter((x) => x != null);
        const spread = vals.length ? Math.max(...vals) - Math.min(...vals) : 0;
        const all = vis.map((x) => valueOf(x, cat)).filter((x) => x != null);
        return { cat, means, spread, mean: all.length ? all.reduce((a, b) => a + b, 0) / all.length : 0 };
      }).filter((r) => r.mean > 0).sort((a, b) => b.spread - a.spread || b.mean - a.mean);
      const top = rows.slice(0, 15);
      // Category names to the left of the dots, or above them on a phone
      const plot = (W) => {
        const narrow = W < 600;
        const labelW = narrow ? 8 : Math.min(290, Math.round(W * 0.38)), rowH = narrow ? 40 : 24, padT = 26;
        const H = padT + top.length * rowH + 10;
        const max = Math.max(1, ...top.map((r) => Math.max(...r.means.filter((x) => x != null))));
        const x = (v) => labelW + v / max * (W - labelW - 20);
        const svg = s('svg', { viewBox: '0 0 ' + W + ' ' + H, width: W, height: H, role: 'img', 'aria-label': title });
        const step = niceStep(max / (narrow ? 3 : 5));
        for (let t2 = 0; t2 <= max; t2 += step) {
          svg.append(s('g', { class: 'grid' }, s('line', { x1: x(t2), x2: x(t2), y1: padT - 6, y2: H - 6 })), s('text', { x: x(t2), y: padT - 10, 'text-anchor': 'middle', class: 'axis' }, fmt(Math.round(t2 * 10) / 10)));
        }
        top.forEach((r, k) => {
          const y = padT + k * rowH + (narrow ? 28 : rowH / 2);
          if (narrow) svg.append(s('text', { x: 0, y: y - 12 }, label(r.cat)));
          else svg.append(s('text', { x: labelW - 10, y: y + 4, 'text-anchor': 'end' }, label(r.cat)));
          const vals = r.means.filter((v) => v != null);
          if (vals.length > 1) svg.append(s('line', { x1: x(Math.min(...vals)), x2: x(Math.max(...vals)), y1: y, y2: y, stroke: 'var(--line-strong)', 'stroke-width': 2 }));
          r.means.forEach((v, li) => {
            if (v == null) return;
            svg.append(s('g', { class: 'hit' }, s('circle', { class: 'mark', cx: x(v), cy: y, r: 5, fill: lineageColor(lins[li]), stroke: 'var(--bg-elev)', 'stroke-width': 2 }),
              s('title', null, lins[li] + ', ' + label(r.cat) + ': ' + (Math.round(v * 10) / 10) + ' per genome')));
          });
        });
        return svg;
      };
      const figBox = h('div', { class: 'fig' });
      responsive(figBox, plot);
      const legend = h('div', { class: 'legend', style: { margin: '0 0 8px' } }, lins.map((l) => lineageChip(l)));
      const tbl = table({
        columns: [{ key: 'cat', label: 'Category', value: (r) => label(r.cat), wrap: true }]
          .concat(lins.map((l, li) => ({ key: 'l' + li, label: l, num: true, value: (r) => (r.means[li] == null ? null : Math.round(r.means[li] * 10) / 10) })))
          .concat([{ key: 'spread', label: 'Spread', num: true, value: (r) => Math.round(r.spread * 10) / 10 }]),
        rows, unit: 'category', csv: fileName + '.csv', pageSize: 10, sort: 'spread', dir: -1,
      });
      tbl.style.marginTop = '12px';
      return card(title, desc, [svgButton(() => figBox.querySelector('svg'), fileName + '.svg')], legend, figBox, tbl);
    }

    function draw() {
      clear(box);
      const m = D.meta.modes;
      const fx = D.functional || {};
      let any = false;
      if (m.mycobrowser && fx.myco_cats) {
        any = true;
        // Lost genes (missing or pseudogene) per category, per genome
        const catOf = REF.cat;
        const cache = new Map();
        const lost = (smp, cat) => {
          let key = smp.i;
          if (!cache.has(key)) {
            const c = new Map();
            for (let gi = 0; gi < NREF; gi++) {
              const st = statusAt(smp.i, gi);
              if (st === ST.M.code || st === ST.P.code) c.set(catOf[gi], (c.get(catOf[gi]) || 0) + 1);
            }
            cache.set(key, c);
          }
          return cache.get(key).get(cat) || 0;
        };
        box.append(dotPlot('Lost genes by Mycobrowser category', 'Reference genes missing or flagged as pseudogenes, per genome, by the reference gene\'s functional category in Mycobrowser. The categories where the lineages differ most come first; a lineage\'s dot is its mean per genome.',
          fx.myco_cats, lost, (c) => c, 'amap_lost_by_category'));
      }
      if (m.functional && fx.cog) {
        any = true;
        const cats = [...new Set(Object.values(fx.cog).flatMap((c) => Object.keys(c)))].sort();
        const el = dotPlot('COG categories', 'Genes per COG category from eggNOG, per genome. A gene with several letters counts in each. The categories where the lineages differ most come first.',
          cats, (smp, cat) => (fx.cog[smp.name] ? fx.cog[smp.name][cat] || 0 : null), (c) => c + ' · ' + (COG_NAMES[c] || c), 'amap_cog');
        el.style.marginTop = '16px';
        box.append(el);
        if (fx.kegg) {
          const paths = [...new Set(Object.values(fx.kegg).flatMap((c) => Object.keys(c)))];
          const el2 = dotPlot('KEGG pathways', 'Genes per KEGG pathway from eggNOG, per genome, the pathways that differ most between lineages first. Pathways are given by their KEGG map number.',
            paths, (smp, p) => (fx.kegg[smp.name] ? fx.kegg[smp.name][p] || 0 : null), (p) => p, 'amap_kegg');
          el2.style.marginTop = '16px';
          box.append(el2);
        }
      }
      if (!m.functional) {
        const c = callout('absent', 'No functional annotation.', m.skip_functional ?
          'This run used --skip_functional_annotation, so there are no COG categories or KEGG pathways from eggNOG.' :
          'The report was written without master tables, so there are no COG categories or KEGG pathways.');
        if (any) c.style.marginTop = '16px';
        box.append(c);
      }
      if (!m.mycobrowser) {
        const c = callout('absent', 'No Mycobrowser categories.', m.skip_mycobrowser ? 'Mycobrowser describes MTBC genes only, and this run skipped it (--skip_mycobrowser).' : 'The report was written without the Mycobrowser table of the reference.');
        c.style.marginTop = '10px';
        box.append(c);
      }
      if (LINEAGES.length < 2 && any) {
        const c = callout('', 'One lineage.', 'With every genome in one lineage the categories are compared with nothing; the table still gives each category\'s mean per genome.');
        c.style.marginTop = '10px';
        box.append(c);
      }
    }
    onFilter(draw);
    draw();
  }

  // ── 8. Methods ─────────────────────────────────────────────────────────

  function renderMethodsSection() {
    const body = clear(document.getElementById('methods-body'));
    const m = D.meta;
    const p = m.params || {};
    const run = m.run || {};
    const ref = m.reference;
    const kv = (pairs) => h('dl', { class: 'kv' }, pairs.filter((x) => x && x[1] != null && x[1] !== '').map(([k, v]) => [h('dt', null, k), h('dd', null, v)]));
    body.append(h('div', { class: 'grid-2' },
      card('This run', null, null, kv([
        ['AMAP', m.version ? 'v' + String(m.version).replace(/^v/, '') : 'not given'],
        ['Report written', m.generated],
        ['Run started', run.run || ''],
        ['Nextflow', run.nextflow || ''],
        ['Profile', run.profile || ''],
        ['Command', run.command ? h('code', null, run.command) : ''],
        ['Genomes', fmt(OK.length) + (FAILED.length ? ' annotated, ' + fmt(FAILED.length) + ' without annotation' : '')],
      ])),
      card('The reference', null, null, kv([
        ['Reference', ref.custom === true ? 'custom: ' + (ref.annotation || '') : ref.custom === false ? 'H37Rv, bundled in the container' : 'not recorded: the report was written without the run\'s params.yml'],
        ref.custom === true ? ['Genome', ref.genome] : null,
        ['Genes', fmt(ref.genes) + ' proteins in the reference proteome'],
        ['Length', ref.length ? kb(ref.length) + ' (the last gene\'s end)' : ''],
        ['Organism', [p.genus, p.species].filter(Boolean).join(', ')],
        ['Mycobrowser', m.modes.mycobrowser ? 'used for names, products, categories and drug resistance' : 'not used'],
      ]))));

    const defs = [
      [stateChip(null), 'A gene of the genome names the reference gene and nothing sets it apart: same length, one copy, not a pseudogene, not cut by a contig end.'],
      [stateChip('M'), 'No gene of the genome names the reference gene. A region of difference shows as a block of these.'],
      [stateChip('P'), 'The best gene naming it is flagged as a pseudogene: an internal stop codon, a lost stop codon or a frame that runs on past the reference protein. The note says which.'],
      [stateChip('T'), 'The best gene naming it runs into a contig end and is cut there (draft assemblies). Not a pseudogene.'],
      [stateChip('C'), 'More than one gene names it alone: AMAP kept further copies, each at least 90% as long as the first.'],
      [stateChip('L'), 'One intact copy, of a different length from the reference gene: an indel, or another start codon.'],
      [stateChip('S'), 'Only a gene that names several reference genes at once names it, such as an IS6110 copy that matches all seventeen. Counted as found.'],
    ];
    body.append(h('div', { class: 'card', style: { 'margin-top': '16px' } }, h('h3', { style: { 'margin-bottom': '12px' } }, 'What each state means'),
      h('div', { class: 'defs' }, defs.map(([a, b]) => h('div', null, h('div', null, a), h('div', null, b)))),
      h('p', { class: 'muted', style: { 'margin-top': '14px' } }, 'Each reference gene gets one state per genome, from the genes that name it in that genome\'s final annotation (', h('code', null, 'AMAP/<name>/<name>.AMAP.tsv'), '): a pseudogene first, then a partial gene, then copies, then a length change. "Found" counts every state but missing.')));

    const params = Object.entries(p);
    if (params.length) {
      const el = card('Parameters', 'Every parameter of the run, as DUMP_PARAMS wrote them (pipeline_info/params.yml).', null,
        table({ columns: [{ key: 'k', label: 'Parameter' }, { key: 'v', label: 'Value', wrap: true }], rows: params.map(([k, v]) => ({ k, v })), unit: 'parameter', csv: 'amap_params.csv', pageSize: 15 }));
      el.style.marginTop = '16px';
      body.append(el);
    } else {
      const c = callout('absent', 'No parameters.', 'The report was written without pipeline_info/params.yml.');
      c.style.marginTop = '16px';
      body.append(c);
    }
    if (m.versions && m.versions.length) {
      const el = card('Software', 'The version of every tool, as SOFTWARE_VERSIONS found them in the container (pipeline_info/software_versions.tsv).', null,
        table({ columns: [{ key: 't', label: 'Tool' }, { key: 'v', label: 'Version' }], rows: m.versions.map(([t, v]) => ({ t, v })), unit: 'tool', csv: 'amap_software_versions.csv', pageSize: 20, search: false }));
      el.style.marginTop = '16px';
      body.append(el);
    } else {
      const c = callout('absent', 'No software versions.', 'The report was written without pipeline_info/software_versions.tsv.');
      c.style.marginTop = '10px';
      body.append(c);
    }
    const c2 = callout('', 'How to cite.', 'Ruiz-Rodriguez P, Coscolla M. AMAP: Assisted MTBC Annotation Pipeline. Zenodo, doi:10.5281/zenodo.8024446. Cite Bakta, miniprot and, if used, eggNOG-mapper and pangene as well.');
    c2.style.marginTop = '16px';
    body.append(c2);
  }

  // ── Boot ───────────────────────────────────────────────────────────────

  function setupChrome() {
    tip.el = document.getElementById('tip');
    document.getElementById('theme').addEventListener('click', () => {
      const next = theme() === 'light' ? 'dark' : 'light';
      if (next === 'light') document.documentElement.dataset.theme = 'light';
      else delete document.documentElement.dataset.theme;
      try { localStorage.setItem('amap.theme', next); } catch (e) { /* private window: the choice lasts the visit */ }
      renderFilterBar();
      for (const fn of listeners) fn();
    });
    document.getElementById('sample-filter').addEventListener('input', (e) => {
      state.text = e.target.value.trim();
      filtersChanged();
    });
    document.getElementById('filters-reset').addEventListener('click', () => {
      state.lineages = new Set(LINEAGES);
      state.text = '';
      document.getElementById('sample-filter').value = '';
      filtersChanged();
    });
    document.getElementById('drawer-close').addEventListener('click', closeDrawer);
    document.getElementById('drawer').addEventListener('click', (e) => { if (e.target.id === 'drawer') closeDrawer(); });
    document.addEventListener('keydown', (e) => { if (e.key === 'Escape' && document.getElementById('drawer').classList.contains('open')) closeDrawer(); });
    // The section in view, marked in the top bar
    const links = [...document.querySelectorAll('.nav a')];
    const io = new IntersectionObserver((entries) => {
      for (const en of entries) {
        if (!en.isIntersecting) continue;
        for (const a of links) a.classList.toggle('on', a.getAttribute('href') === '#' + en.target.id);
      }
    }, { rootMargin: '-45% 0px -50% 0px' });
    for (const sec of document.querySelectorAll('main section[id]')) io.observe(sec);
  }

  function fail(err) {
    const main = document.getElementById('main');
    clear(main).append(h('div', { class: 'wrap', style: { padding: '60px 0' } },
      h('div', { class: 'callout warn' }, h('p', null, h('b', null, 'The report could not be drawn. '), String(err && err.message ? err.message : err)))));
  }

  async function boot() {
    setupChrome();
    let data;
    try {
      data = await loadData();
    } catch (err) {
      fail(err);
      return;
    }
    buildModel(data);
    state.lineages = new Set(LINEAGES);
    renderFilterBar();
    renderHero();
    if (!OK.length) {
      // Nothing was annotated: say so once, and keep the methods
      clear(document.getElementById('hero-tiles'));
      document.getElementById('hero-lead').textContent = 'No genome of this run produced an annotation, so there is nothing to compare.' +
        (FAILED.length ? ' The genomes with a case summary but no final annotation are named below.' : '');
      for (const id of ['samples-body', 'map-body', 'lookup-body', 'correction-body', 'pseudo-body', 'beyond-body', 'function-body']) {
        clear(document.getElementById(id)).append(callout('absent', 'No genome to show.', ''));
      }
      renderMethodsSection();
      return;
    }
    const sections = [renderSamples, renderMapSection, renderLookupSection, renderCorrectionSection, renderPseudoSection, renderBeyondSection, renderFunctionSection, renderMethodsSection];
    for (const fn of sections) {
      try {
        fn();
      } catch (err) {
        console.error(err);
      }
    }
    onFilter(renderSamples);
    if (location.hash && document.querySelector(location.hash)) document.querySelector(location.hash).scrollIntoView();
    window.AMAP_REPORT = { data: D, genesOf, goneOf, statusAt, detailAt, ok: OK };
  }

  if (typeof module === 'object' && module.exports) {
    // Loaded by node, for the tests: the decoders, and no page to draw
    module.exports = { unrle, decodeGenes, decodeGone };
  } else if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', boot);
  else boot();
})();
