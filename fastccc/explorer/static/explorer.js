/* FastCCC result explorer — renders the embedded payload with D3 (no server). */
(function () {
  'use strict';
  const P = window.FASTCCC_PAYLOAD;
  const d3 = window.d3;
  const KIND = P.kind;                       // standard | split | reference
  const CT = P.celltypes;                    // [{name, n_cells, color}]
  const LR = P.lr;                           // [{id, ligand, receptor, classification, directionality}]
  const R = P.rows;                          // columnar arrays
  const N = R.sender.length;
  const TREND = { 0: 'Not in reference', 1: 'Up', 2: 'Down', 3: 'Both Sig', 4: 'Both NS' };
  const TREND_COLOR = { 1: '#d1495b', 2: '#3b6fb6', 3: '#8c7b6b', 4: '#c9c9c9', 0: '#e3e3e3' };
  const COND_COLORS = ['#8c7b6b', '#d1495b', '#3b6fb6', '#e9a44b', '#5aa9a0'];
  const PATTERN_COLORS = ['#d1495b', '#3b6fb6', '#e9a44b', '#5aa9a0', '#9b8bc4', '#a8c36b', '#c45b8a', '#7b9ac9'];
  const SINGLE = '#4b6fa5';
  const LIG_INK = '#2f5d9e', REC_INK = '#c2385a';
  const FONT = 'system-ui,-apple-system,"Segoe UI",Helvetica,Arial,sans-serif';
  const fmtP = p => (p == null ? '–' : p < 1e-4 ? p.toExponential(1) : p.toFixed(4));
  const fmt2 = v => (v == null ? '–' : (+v).toFixed(2));
  const pct = v => (v == null ? '–' : Math.round(v * 100) + '%');
  const nlogp = p => (p == null ? 0 : Math.min(16, -Math.log10(Math.max(p, 1e-16))));
  const esc = s => String(s).replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
  const classOf = lr => lr.classification || 'Unclassified';
  const dbLabel = m => {                   // "CPDBv5.0.0 (human, 2,916 interactions)" from the database metadata
    const d = m.database_info || {}; const name = m.database || d.name || '';
    const ver = d.version == null ? null : String(d.version);
    const parts = [d.species, ver && !name.includes(ver) ? (/^v?\d+(\.\d+)*$/.test(ver) ? 'v' + ver.replace(/^v/i, '') : ver) : null,
                   d.n_interactions ? d.n_interactions.toLocaleString() + ' interactions' : null].filter(Boolean);
    return name + (parts.length ? ` (${parts.join(', ')})` : '');
  };

  // ── derived per-LR statistics (global, independent of filters) ─────────────
  const lrSigPairs = new Array(LR.length).fill(0);      // specificity: # cell-type pairs where significant
  for (let i = 0; i < N; i++) {
    const k = R.lr[i];
    if (KIND === 'reference') { if (R.trend[i] === 1 || R.trend[i] === 2 || R.trend[i] === 3) lrSigPairs[k]++; }
    else if (R.cond == null || R.cond[i] === 0) lrSigPairs[k]++;
  }
  const classes = Array.from(new Set(LR.map(classOf))).sort();
  const DB_HAS_CLASSES = LR.some(l => l.classification);
  const DB_SPECIES = (P.meta.database_info && P.meta.database_info.species) || 'human';
  // gene-set membership of LR pairs: a pair belongs to every set that contains any of its ligand or receptor genes
  const lrSetsCache = {};
  function lrSets(libIdx) {
    if (lrSetsCache[libIdx]) return lrSetsCache[libIdx];
    const lib = P.genesets.libraries[libIdx]; const gpos = new Map(P.genesets.genes.map((g, i) => [g, i]));
    const geneLR = new Map(); LR.forEach((l, k) => l.ligand.split(',').concat(l.receptor.split(',')).forEach(g => { const gi = gpos.get(g.trim()); if (gi != null) { if (!geneLR.has(gi)) geneLR.set(gi, new Set()); geneLR.get(gi).add(k); } }));
    const out = LR.map(() => []);
    lib.sets.forEach((s, si) => { const hit = new Set(); s.genes.forEach(gi => { const ks = geneLR.get(gi); if (ks) ks.forEach(k => hit.add(k)); }); hit.forEach(k => out[k].push(si)); });
    return (lrSetsCache[libIdx] = out);
  }
  const lrInSet = (k, libIdx, si) => lrSets(libIdx)[k].includes(si);
  const modules = P.modules;
  const nModules = modules ? modules.names.length : 0;
  const conditions = KIND === 'split' ? P.meta.conditions.map(c => c.name) : null;
  // scoring modules: "<summary>_<complex aggregation>_<LR combination>", e.g. Quantile_0.9_Minimum_Arithmetic
  const MOD = modules ? modules.names.map(n => { const parts = n.split('_'); const comb = parts.pop(); const agg = parts.pop(); return { summary: parts.join('_'), agg, comb }; }) : [];
  // factor levels come from the run itself: FastCCC accepts Mean, Median, Q2, Q3 and any Quantile_<q> summary
  const sumQ = s => (s === 'Median' || s === 'Q2') ? 0.5 : s === 'Q3' ? 0.75 : s.startsWith('Quantile_') ? parseFloat(s.slice(9)) : null;
  const present = key => Array.from(new Set(MOD.map(m => m[key])));
  const ordered = (known, xs) => known.filter(k => xs.includes(k)).concat(xs.filter(x => !known.includes(x)));
  const SUMMARIES = present('summary').sort((a, b) => (sumQ(a) ?? -1) - (sumQ(b) ?? -1));
  const AGGS = ordered(['Minimum', 'Average'], present('agg')), COMBS = ordered(['Arithmetic', 'Geometric'], present('comb'));
  const sumLabel = s => s === 'Mean' ? 'Mean' : s === 'Median' ? 'Median' : sumQ(s) != null ? 'Q' + Math.round(sumQ(s) * 100) : s;
  const sumLong = s => s === 'Mean' ? 'mean' : s === 'Median' ? 'median' : sumQ(s) != null ? `${Math.round(sumQ(s) * 100)}th percentile` : s;
  const sumHelp = s => { if (s === 'Mean') return 'average over all cells; a few very high cells can carry it'; const q = sumQ(s); if (q == null) return ''; if (q === 0.5) return 'a typical cell; most cells must express (broad signal)'; return `${Math.round(q * 100)}th percentile; the top ${Math.round((1 - q) * 100)}% of the cells suffice${q >= 0.9 ? ' (sub-population signal)' : ''}`; };
  const AGG_SHORT = { Minimum: 'Min', Average: 'Avg' }, AGG_HELP = { Minimum: 'every subunit of a complex must pass', Average: 'subunits are averaged' };
  const COMB_SHORT = { Arithmetic: 'Arith', Geometric: 'Geo' }, COMB_HELP = { Arithmetic: 'one strongly expressed side can carry the pair', Geometric: 'both ligand and receptor must be expressed' };
  // the signal-type reading needs the full design: a mean, a median-type (q = 0.5) and an upper-quantile summary, both aggregations, both combinations
  const FULL_DESIGN = SUMMARIES.includes('Mean') && SUMMARIES.some(s => sumQ(s) === 0.5) && SUMMARIES.some(s => (sumQ(s) ?? 0) > 0.5) && AGGS.length === 2 && COMBS.length === 2;
  const SIGNAL_TYPES = ['broad', 'sub-population', 'mean-driven', 'weak'];
  const SIGNAL_COLOR = { 'broad': '#4b6fa5', 'sub-population': '#d1495b', 'mean-driven': '#e9a44b', 'weak': '#bdbab2' };
  const SIGNAL_HELP = {
    'broad': 'the median-based modules are significant: a typical sender cell expresses the ligand and a typical receiver cell the receptor',
    'sub-population': 'only the upper-quantile modules are significant: the signal is carried by a fraction of the cells',
    'mean-driven': 'only the mean-based modules are significant: a few very high-expressing cells pull the mean',
    'weak': 'few modules agree; the Cauchy combination is significant but no single summary dominates',
  };
  function signalTags(i) {
    if (!modules) return null;
    const ps = modules.p[i]; const sig = k => ps[k] != null && ps[k] < 0.05;
    const count = pred => MOD.reduce((a, m, k) => a + ((pred(m) && sig(k)) ? 1 : 0), 0);
    const nMed = count(m => sumQ(m.summary) === 0.5), nQ = count(m => (sumQ(m.summary) ?? 0) > 0.5), nMean = count(m => m.summary === 'Mean');
    const nMin = count(m => m.agg === 'Minimum'), nAvg = count(m => m.agg === 'Average'), nAri = count(m => m.comb === 'Arithmetic'), nGeo = count(m => m.comb === 'Geometric');
    const lr = LR[R.lr[i]]; const isComplex = lr.ligand.includes(',') || lr.receptor.includes(',');
    const tags = [nMed >= 2 ? 'broad' : (nQ >= 2 ? 'sub-population' : (nMean >= 1 ? 'mean-driven' : 'weak'))];
    if (isComplex && nAvg >= 3 && nMin <= 1) tags.push('complex-limited');
    if (nAri >= 3 && nGeo <= 1) tags.push('one-sided');
    return tags;
  }
  const SIG = modules && FULL_DESIGN ? Array.from({ length: N }, (_, i) => signalTags(i)) : null;
  function quantile(values, q) { const a = values.filter(v => v != null && !Number.isNaN(v)).sort((x, y) => x - y); if (!a.length) return 0; const pos = (a.length - 1) * q, lo = Math.floor(pos), hi = Math.ceil(pos); return a[lo] + (a[hi] - a[lo]) * (pos - lo); }
  const DIST = (() => {
    const cs = R.cs, pairs = R.lr.map(k => lrSigPairs[k]);
    const pctRow = [];
    for (let i = 0; i < N; i++) {
      if (KIND === 'reference') pctRow.push(Math.max(Math.min(R.lig_pct[i] ?? 0, R.rec_pct[i] ?? 0), Math.min(R.lig_pct_ref[i] ?? 0, R.rec_pct_ref[i] ?? 0)));
      else pctRow.push(null);
    }
    // d3.max rather than Math.max(...array): spreading the per-interaction arrays overflows the call stack above ~65k interactions
    return { csMax: quantile(cs, 0.99) || 1, pairsMax: Math.max(1, d3.max(pairs) || 0), pairs, pctRow,
             nlogpMax: KIND === 'reference' ? 0 : Math.ceil(Math.max(1.3, d3.max(R.p, p => nlogp(p)) || 0)) };
  })();

  // ── state ──────────────────────────────────────────────────────────────────
  const state = {
    senders: new Set(CT.map((_, i) => i)), receivers: new Set(CT.map((_, i) => i)),
    trends: new Set(KIND === 'reference' ? (R.trend.some(t => t === 0) ? [0, 1, 2, 3] : [1, 2, 3]) : []), cond: 0,
    gene: '', classes: new Set(classes), minCS: 0, maxPairs: Infinity, minModules: 0, minPct: KIND === 'reference' ? 0 : 0.1, maxP: 0.05, modSummary: new Set(SUMMARIES), modAgg: new Set(AGGS), modComb: new Set(COMBS), modSel: new Set(MOD.map((_, k) => k)), modMask: 0,
    pair: null, sort: { key: 'cs', dir: -1 }, page: 0, selected: null, pinned: new Set(),
    chordBy: KIND === 'reference' ? 'trend' : 'sender', dotTop: 30, dotCols: 'auto', onlyPinned: false,
    rolesBy: 'count', patDir: 'outgoing', patK: 'auto', enrLib: 0, pathBy: 'cpdb', pathRows: KIND === 'reference' ? 'percell' : 'overall', pathRank: 'count', pathTop: 30, setFilter: null, scatterBy: 'fold', cellsGene: null, diffA: 1, diffB: 2, diffPair: null,
    signals: new Set(SIGNAL_TYPES.concat(['complex-limited', 'one-sided'])), minZ: 0, minBS: 0, minRep: 0,
  };
  const HAS_Z = !!(R.z && R.z.some(v => v != null)), HAS_BS = !!(R.bs && R.bs.some(v => v != null)), HAS_REP = !!(R.rep && R.rep.some(v => v != null));
  const pageSize = 50;
  // scoring-module selection: a bit per module; modCount(i) = significant modules among the selected ones
  const maskOf = sel => { let m = 0; sel.forEach(k => { m |= (1 << k); }); return m; };
  const deriveModSel = () => new Set(MOD.map((m, k) => (state.modSummary.has(m.summary) && state.modAgg.has(m.agg) && state.modComb.has(m.comb)) ? k : -1).filter(k => k >= 0));
  state.modMask = maskOf(state.modSel);
  const sigMask = modules ? modules.p.map(ps => ps.reduce((m, p, k) => (p != null && p < 0.05) ? (m | (1 << k)) : m, 0)) : null;
  const popcount = x => { let c = 0; while (x) { x &= x - 1; c++; } return c; };
  const modCount = i => (sigMask ? popcount(sigMask[i] & state.modMask) : null);
  const nSel = () => state.modSel.size;

  // ── expression lookup ──────────────────────────────────────────────────────
  const geneIndex = new Map(); const ctIndex = new Map();
  if (P.genes) { P.genes.genes.forEach((g, i) => geneIndex.set(g, i)); P.genes.celltypes.forEach((c, i) => ctIndex.set(c, i)); }
  function genePct(gene, ctName) { const gi = geneIndex.get(gene), ci = ctIndex.get(ctName); return (gi == null || ci == null) ? null : P.genes.pct[ci][gi]; }
  function complexPct(genes, ctName) { const vals = genes.split(',').map(g => genePct(g.trim(), ctName)).filter(v => v != null); return vals.length ? Math.min(...vals) : null; }
  function exprOf(i) {
    if (!P.genes) return null;
    const lr = LR[R.lr[i]];
    const lig = complexPct(lr.ligand, CT[R.sender[i]].name), rec = complexPct(lr.receptor, CT[R.receiver[i]].name);
    return (lig == null && rec == null) ? null : { lig: lig ?? 0, rec: rec ?? 0 };
  }

  // ── filtering ──────────────────────────────────────────────────────────────
  function passesCore(i, ignoreCond, ignoreTrend) {
    if (!state.senders.has(R.sender[i]) || !state.receivers.has(R.receiver[i])) return false;
    if (KIND === 'reference' && !ignoreTrend && !state.trends.has(R.trend[i])) return false;
    if (KIND === 'split' && !ignoreCond && R.cond[i] !== state.cond) return false;
    if (KIND !== 'reference' && R.p[i] != null && R.p[i] > state.maxP) return false;
    const lr = LR[R.lr[i]];
    if (!state.classes.has(classOf(lr))) return false;
    if (state.setFilter && !lrInSet(R.lr[i], state.setFilter.lib, state.setFilter.set)) return false;
    if (state.gene) { const g = state.gene.toLowerCase(); if (!(lr.ligand.toLowerCase().includes(g) || lr.receptor.toLowerCase().includes(g))) return false; }
    if (R.cs[i] != null && R.cs[i] < state.minCS) return false;
    if (lrSigPairs[R.lr[i]] > state.maxPairs) return false;
    if (modules && state.minModules > 0 && nSel() > 0 && modCount(i) < state.minModules) return false;  // module support (0 = off): at least minModules of the selected modules significant on their own
    if (SIG && !SIG[i].every(t => state.signals.has(t))) return false;
    if (HAS_Z && state.minZ > 0 && (R.z[i] == null || R.z[i] < state.minZ)) return false;
    if (HAS_BS && state.minBS > 0 && (R.bs[i] == null || R.bs[i] < state.minBS)) return false;
    if (HAS_REP && state.minRep > 0 && (R.rep[i] == null || R.rep[i] < state.minRep)) return false;
    if (state.minPct > 0) {
      if (KIND === 'reference') { const q = Math.min(R.lig_pct[i] ?? 0, R.rec_pct[i] ?? 0), rf = Math.min(R.lig_pct_ref[i] ?? 0, R.rec_pct_ref[i] ?? 0); if (Math.max(q, rf) < state.minPct) return false; }
      else if (P.genes) { const e = exprOf(i); if (e && Math.min(e.lig, e.rec) < state.minPct) return false; }
    }
    return true;
  }
  function passes(i) {
    if (state.onlyPinned) return state.pinned.has(i);
    if (!passesCore(i, false)) return false;
    if (state.pair && (R.sender[i] !== state.pair[0] || R.receiver[i] !== state.pair[1])) return false;
    return true;
  }
  let filtered = [], shown = [];
  // `shown`: rows passing the filters (Both NS pairs included when ticked) -> table, evidence, query-vs-reference scatter.
  // `filtered`: the same minus Both NS -> every panel that counts or ranks interactions (a Both NS pair is not a communication).
  function refilter() { filtered = []; shown = []; for (let i = 0; i < N; i++) if (passes(i)) { shown.push(i); if (!(KIND === 'reference' && R.trend[i] === 4)) filtered.push(i); } }
  function rowPct(i) { if (KIND === 'reference') return DIST.pctRow[i]; const e = exprOf(i); return e ? Math.min(e.lig, e.rec) : null; }
  const PCT_Q = (() => { const v = []; for (let i = 0; i < N; i += Math.max(1, Math.floor(N / 20000))) { const x = rowPct(i); if (x != null) v.push(x); } return q => quantile(v, q); })();
  const PAIRS_Q = q => quantile(DIST.pairs, q);
  function setSlider(id, vid, value, fmt) { const s = document.getElementById(id); s.value = value; document.getElementById(vid).textContent = fmt(s.value); }
  function applyHighlights() {
    const levels = [[0.05, 0.7, 1], [0.1, 0.6, 0.9], [0.25, 0.5, 0.75], [0.35, 0.4, 0.6], [0.5, 0.3, 0.5], [0.75, 0.2, 0.5], [1, 0, 0]];
    for (const [pq, eq, mf] of levels) {
      state.maxPairs = Math.max(1, Math.ceil(PAIRS_Q(pq)));
      state.minPct = Math.max(DEFAULT.minPct, Math.floor(PCT_Q(eq) * 20) / 20);
      state.minModules = modules ? Math.ceil(nSel() * mf) : 0;
      state.minBS = HAS_BS ? (mf >= 0.75 ? 0.8 : mf >= 0.5 ? 0.6 : 0) : 0;
      if (KIND === 'reference') assign(state.trends, [1, 2]);
      refilter();
      if (filtered.length >= 20) break;
    }
  }
  function highlightNote() {
    return `for this dataset: LR pair significant in ≤ ${state.maxPairs} cell-type pair${state.maxPairs === 1 ? '' : 's'}, ≥ ${Math.round(state.minPct * 100)}% expressing cells` +
      (modules ? `, ≥ ${state.minModules}/${nSel()} modules` : '') + (HAS_BS && state.minBS ? `, bootstrap support ≥ ${Math.round(state.minBS * 100)}%` : '') + (KIND === 'reference' ? ', Up or Down only' : '');
  }
  // ── presets: Default = the run's own criteria; Highlights = Default + adaptive thresholds; Custom = the last manual settings ──
  const FILTER_KEYS = ['gene', 'minCS', 'maxPairs', 'minModules', 'minPct', 'maxP', 'minZ', 'minBS', 'minRep', 'cond', 'pair', 'setFilter'];
  const SET_KEYS = ['senders', 'receivers', 'trends', 'classes', 'signals', 'modSummary', 'modAgg', 'modComb'];
  function snapshot() { const s = {}; FILTER_KEYS.forEach(k => { s[k] = state[k]; }); SET_KEYS.forEach(k => { s[k] = new Set(state[k]); }); return s; }
  const assign = (target, values) => { target.clear(); values.forEach(v => target.add(v)); return target; };  // mutate in place: the checkbox handlers hold references to these Sets
  function applySnapshot(s) { FILTER_KEYS.forEach(k => { state[k] = s[k]; }); SET_KEYS.forEach(k => { assign(state[k], s[k]); }); state.modSel = deriveModSel(); state.modMask = maskOf(state.modSel); }
  const DEFAULT = snapshot();
  let preset = 'default', custom = null, applying = false, syncControls = () => {};
  function userChange() { if (!applying) { preset = 'custom'; custom = snapshot(); } update(); renderPresets(); }
  function setPreset(name) {
    applying = true;
    if (name === 'reset') { custom = null; name = 'default'; }
    applySnapshot(name === 'custom' ? (custom || DEFAULT) : DEFAULT);
    if (name === 'highlights') applyHighlights();
    syncControls(); preset = name; applying = false; update(); renderPresets();
  }
  function renderPresets() {
    ['default', 'highlights', 'custom'].forEach(n => document.getElementById('f-' + n).classList.toggle('primary', preset === n));
    document.getElementById('f-custom').classList.toggle('dim', !custom);
    const n = `${filtered.length.toLocaleString()} interaction${filtered.length === 1 ? '' : 's'}`;
    let note;
    if (preset === 'highlights') note = `${highlightNote()} → ${n}`;
    else if (preset === 'custom') note = `your own filters → ${n}`;
    else if (KIND === 'reference') note = `Up, Down and Both Sig interactions of all cell types → ${n}`;
    else note = `all significant interactions (p < 0.05; ligand and receptor in ≥ ${Math.round(DEFAULT.minPct * 100)}% of the cells) → ${n}`;
    document.getElementById('f-preset-note').textContent = note;
  }

  // ── tooltip & svg helpers ──────────────────────────────────────────────────
  const tip = d3.select('body').append('div').attr('class', 'tip').style('display', 'none');
  function showTip(ev, html) { tip.style('display', 'block').html(html).style('left', (ev.pageX + 14) + 'px').style('top', (ev.pageY + 12) + 'px'); }
  function hideTip() { tip.style('display', 'none'); }
  function rowTip(i) {
    const lr = LR[R.lr[i]];
    let h = `<b>${esc(CT[R.sender[i]].name)}</b> → <b>${esc(CT[R.receiver[i]].name)}</b><br><span class="lig">${esc(lr.ligand)}</span> – <span class="rec">${esc(lr.receptor)}</span><br>`;
    if (KIND === 'reference') h += `trend: <b>${TREND[R.trend[i]]}</b> · CS ${fmt2(R.cs[i])}` + (R.ref_cs[i] != null ? ` (reference ${fmt2(R.ref_cs[i])})` : '') + `<br>ligand ${pct(R.lig_pct[i])} vs ref ${pct(R.lig_pct_ref[i])} · receptor ${pct(R.rec_pct[i])} vs ref ${pct(R.rec_pct_ref[i])}`;
    else h += `p = ${fmtP(R.p[i])} · CS ${fmt2(R.cs[i])}` + (modules ? ` · ${modules.n_sig[i]}/${nModules} modules` : '');
    h += `<br>${esc(classOf(lr))} · significant in ${lrSigPairs[R.lr[i]]} cell-type pair${lrSigPairs[R.lr[i]] === 1 ? '' : 's'}`;
    return h;
  }
  const SVGNS = 'http://www.w3.org/2000/svg';
  const SVG_CSS = `text{font-family:${FONT};font-size:11px;fill:#333}.axis-lab{font-size:10.5px;fill:#444}.plot-title{font-size:12px;font-weight:600;fill:#333}.chord-lab{font-size:10.5px;fill:#333}`;
  function newSvg(host, W, H, id) {
    const svg = host.append('svg').attr('viewBox', `0 0 ${W} ${H}`).attr('width', W).attr('height', H).attr('xmlns', 'http://www.w3.org/2000/svg');
    if (id) svg.attr('id', id);
    svg.append('style').text(SVG_CSS);
    return svg;
  }
  function axisStyle(g) { g.selectAll('text').attr('class', 'axis-lab'); g.selectAll('line,path').attr('stroke', '#c9c7c0'); }
  // Fit every chart to its drawn content: labels that run past the nominal width/height (rotated axis labels, legends,
  // long titles) enlarge the SVG instead of being clipped. Runs after each render via a mutation observer.
  function fitSvg(s) {
    try {
      if (!s.isConnected || s.getClientRects().length === 0) return;
      const b = s.getBBox(); if (!isFinite(b.width) || (b.width === 0 && b.height === 0)) return;
      const pad = 6, vb = (s.getAttribute('viewBox') || '').trim().split(/[\s,]+/).map(Number), hasVb = vb.length === 4 && vb.every(Number.isFinite);
      const vx = hasVb ? vb[0] : 0, vy = hasVb ? vb[1] : 0, W = hasVb ? vb[0] + vb[2] : (+s.getAttribute('width') || 0), H = hasVb ? vb[1] + vb[3] : (+s.getAttribute('height') || 0);
      const x0 = Math.min(vx, Math.floor(b.x - pad)), y0 = Math.min(vy, Math.floor(b.y - pad)), x1 = Math.max(W, Math.ceil(b.x + b.width + pad)), y1 = Math.max(H, Math.ceil(b.y + b.height + pad));
      if (x0 !== vx || y0 !== vy || x1 !== W || y1 !== H) {
        s.setAttribute('viewBox', `${x0} ${y0} ${x1 - x0} ${y1 - y0}`);
        if (/%$/.test(s.getAttribute('width') || '')) s.removeAttribute('height');  // responsive charts (width 100%): height follows the viewBox
        else { s.setAttribute('width', x1 - x0); s.setAttribute('height', y1 - y0); }
      }
    } catch (e) { /* detached or hidden */ }
  }
  let fitPending = false;
  function fitSvgs() { fitPending = false; document.querySelectorAll('main svg').forEach(fitSvg); }
  new MutationObserver(() => { if (!fitPending) { fitPending = true; requestAnimationFrame(fitSvgs); } }).observe(document.querySelector('main'), { childList: true, subtree: true });
  // ── exports ────────────────────────────────────────────────────────────────
  // The link is attached to the page while it is clicked and the object URL is kept for a minute: some browsers
  // ignore a detached link, and cancel a download whose URL is revoked before the file has been read.
  function downloadBlob(blob, name) {
    const a = document.createElement('a'); const url = URL.createObjectURL(blob); a.href = url; a.download = name; a.style.display = 'none';
    document.body.appendChild(a); a.click(); a.remove(); setTimeout(() => URL.revokeObjectURL(url), 60000);
  }
  function svgBox(node) {
    const vb = (node.getAttribute('viewBox') || '').trim().split(/[\s,]+/).map(Number); if (vb.length === 4 && vb.every(Number.isFinite)) return vb;
    const r = node.getBoundingClientRect(); return [0, 0, +node.getAttribute('width') || r.width || 600, +node.getAttribute('height') || r.height || 400];
  }
  // A chart as it goes into a file: a copy with an absolute size (charts that are width 100% on the page would have
  // none on their own) and without what only matters in the page (element id, layout styles, pointer cursors).
  function exportNode(node) {
    const c = node.cloneNode(true), vb = svgBox(node);
    c.setAttribute('xmlns', SVGNS); c.setAttribute('viewBox', vb.join(' ')); c.setAttribute('width', vb[2]); c.setAttribute('height', vb[3]);
    c.removeAttribute('id'); c.removeAttribute('style');
    c.querySelectorAll('[style]').forEach(e => { e.style.removeProperty('cursor'); if (!e.getAttribute('style')) e.removeAttribute('style'); });
    // Label halos use paint-order on the page. Editors that ignore that property (Illustrator, older librsvg) would
    // paint the white outline over the letters, so in the file the outline is its own element under the label.
    c.querySelectorAll('[paint-order]').forEach(e => {
      const stroke = e.getAttribute('stroke'); e.removeAttribute('paint-order');
      if (stroke && stroke !== 'none') { const halo = e.cloneNode(true); halo.style.fill = 'none'; halo.setAttribute('stroke-linejoin', 'round'); e.parentNode.insertBefore(halo, e); }
      e.removeAttribute('stroke'); e.removeAttribute('stroke-width');
    });
    return c;
  }
  // Several charts in one SVG document, built as DOM nodes so that the serialised file is always well-formed XML.
  // Each chart becomes a <g> placed by a transform (not every editor understands nested <svg> elements), and the
  // text style they share is written once.
  function composeSvg(nodes, letters) {
    const items = nodes.map(node => { const vb = svgBox(node); return { node, vb, w: vb[2], h: vb[3] }; });
    const gap = 30, head = letters ? 20 : 0; let W, H;
    if (letters) {      // figure builder: two columns of at most 620 px, each panel under its letter
      const colW = 620, cols = items.length > 1 ? 2 : 1, ys = [0, 0];
      items.forEach(it => { const c = ys[0] <= ys[1] ? 0 : (cols === 2 ? 1 : 0); it.sc = Math.min(1, colW / it.w); it.x = c * (colW + gap); it.y = ys[c]; ys[c] += it.h * it.sc + gap + head; });
      W = cols === 1 ? Math.ceil(items[0].w * items[0].sc) : cols * colW + (cols - 1) * gap; H = Math.ceil(Math.max(ys[0], ys[1]));
    } else {            // the charts of one card: side by side at their own size, as on the page
      let x = 0; items.forEach(it => { it.sc = 1; it.x = x; it.y = 0; x += it.w + gap; });
      W = Math.ceil(x - gap); H = Math.ceil(Math.max(...items.map(it => it.h)));
    }
    const mk = (tag, attrs) => { const e = document.createElementNS(SVGNS, tag); Object.entries(attrs || {}).forEach(([k, v]) => e.setAttribute(k, v)); return e; };
    const root = mk('svg', { xmlns: SVGNS, viewBox: `0 0 ${W} ${H}`, width: W, height: H });
    root.appendChild(mk('style')).textContent = SVG_CSS;
    root.appendChild(mk('rect', { width: W, height: H, fill: 'white' }));
    items.forEach((it, k) => {
      if (letters) root.appendChild(mk('text', { x: it.x, y: it.y + 14, style: 'font-size:16px;font-weight:700;fill:#111' })).textContent = String.fromCharCode(65 + k);
      const g = root.appendChild(mk('g', { transform: `translate(${it.x - it.vb[0] * it.sc},${it.y + head - it.vb[1] * it.sc}) scale(${it.sc})` }));
      Array.from(exportNode(it.node).childNodes).forEach(ch => { if (String(ch.nodeName).toLowerCase() !== 'style') g.appendChild(ch); });
    });
    return { root, W, H };
  }
  function svgFile(root) { return new Blob(['<?xml version="1.0" encoding="UTF-8"?>\n' + new XMLSerializer().serializeToString(root)], { type: 'image/svg+xml' }); }
  // `spec` is a container id, or several separated by commas when a card holds more than one chart. Returns false when nothing is drawn.
  function downloadSvg(spec, name) {
    const ids = String(spec || '').split(',').map(t => t.trim()).filter(Boolean);
    const nodes = ids.map(id => document.querySelector(`#${id} svg`)).filter(Boolean); if (!nodes.length) return false;
    downloadBlob(svgFile(nodes.length === 1 ? exportNode(nodes[0]) : composeSvg(nodes, false).root), `fastccc_${name || ids[0]}.svg`); return true;
  }
  document.querySelectorAll('.btn.svg[data-svg]').forEach(b => b.addEventListener('click', () => {
    if (downloadSvg(b.dataset.svg, b.dataset.name)) return;
    b.textContent = 'nothing drawn'; setTimeout(() => { b.textContent = 'SVG'; }, 1500);
  }));

  // ── stat tiles ─────────────────────────────────────────────────────────────
  function renderTiles() {
    const el = document.getElementById('tiles'); const t = [];
    // every tile reads "n / m": n = under the current filters, m = in the whole dataset (the current condition on split pages)
    const base = KIND === 'split' ? d3.range(N).filter(i => R.cond[i] === state.cond) : d3.range(N);
    const frac = (n, m) => `${n.toLocaleString()}<span class="tile-m"> / ${m.toLocaleString()}</span>`;
    if (KIND === 'reference') {
      const c = { 1: 0, 2: 0, 3: 0, 4: 0, 0: 0 }, cm = { 1: 0, 2: 0, 3: 0, 4: 0, 0: 0 }; shown.forEach(i => c[R.trend[i]]++); base.forEach(i => cm[R.trend[i]]++);
      t.push(['Up', frac(c[1], cm[1]), TREND_COLOR[1], 'significant in query, not in reference']);
      t.push(['Down', frac(c[2], cm[2]), TREND_COLOR[2], 'significant in reference, not in query']);
      t.push(['Both Sig', frac(c[3], cm[3]), TREND_COLOR[3], 'significant in both']);
      t.push(['Both NS', frac(c[4], cm[4]), '#9a9a9a', 'significant in neither']);
      if (cm[0]) t.push(['Not in reference', frac(c[0], cm[0]), '#777', 'significant in the query; cell types not in the panel']);
      const cov = P.meta.celltype_coverage;
      if (cov) t.push(['Coverage', `${cov.n_covered}/${cov.n_query_types}`, '#777', `query cell types in the reference · ${Math.round(cov.cells_covered_fraction * 100)}% of cells`]);
    } else {
      const lrN = new Set(filtered.map(i => R.lr[i])).size, lrM = new Set(base.map(i => R.lr[i])).size;
      const pairN = new Set(filtered.map(i => R.sender[i] * 1000 + R.receiver[i])).size, pairM = new Set(base.map(i => R.sender[i] * 1000 + R.receiver[i])).size;
      const ctN = new Set(filtered.flatMap(i => [R.sender[i], R.receiver[i]])).size;
      t.push(['Interactions', frac(filtered.length, base.length), SINGLE, 'significant sender–receiver × LR pairs']);
      t.push(['LR pairs', frac(lrN, lrM), SINGLE, 'distinct ligand–receptor pairs']);
      t.push(['Cell-type pairs', frac(pairN, pairM), SINGLE, 'sender–receiver pairs with ≥ 1 interaction']);
      t.push(['Cell types', frac(ctN, CT.length), SINGLE, P.meta.n_cells ? `with ≥ 1 interaction · ${P.meta.n_cells.toLocaleString()} cells` : 'with ≥ 1 interaction']);
    }
    el.innerHTML = t.map(([l, v, col, sub]) => `<div class="tile"><div class="tile-v" style="color:${col}">${v}</div><div class="tile-l">${l}</div><div class="tile-s">${sub}</div></div>`).join('');
  }

  // ── filter panel ───────────────────────────────────────────────────────────
  function checkList(id, items, set, onchange, colorOf) {
    const el = document.getElementById(id);
    el.innerHTML = `<div class="cl-tools"><a href="#" data-a="all">all</a> · <a href="#" data-a="none">none</a></div>` +
      items.map((name, i) => `<label class="cl-item"><input type="checkbox" data-i="${i}" ${set.has(i) ? 'checked' : ''}>` + (colorOf ? `<span class="sw" style="background:${colorOf(i)}"></span>` : '') + `<span>${esc(name)}</span></label>`).join('');
    el.querySelectorAll('input').forEach(cb => cb.addEventListener('change', () => { const i = +cb.dataset.i; cb.checked ? set.add(i) : set.delete(i); onchange(); }));
    el.querySelectorAll('a').forEach(a => a.addEventListener('click', ev => { ev.preventDefault(); set.clear(); if (a.dataset.a === 'all') items.forEach((_, i) => set.add(i)); el.querySelectorAll('input').forEach(cb => cb.checked = set.has(+cb.dataset.i)); onchange(); }));
  }
  function syncClassList() { document.querySelectorAll('#f-classes input').forEach(cb => cb.checked = state.classes.has(cb.dataset.c)); }
  function buildFilters() {
    const names = CT.map(c => c.name);
    checkList('f-senders', names, state.senders, userChange, i => CT[i].color);
    checkList('f-receivers', names, state.receivers, userChange, i => CT[i].color);
    const cls = document.getElementById('f-classes');
    cls.innerHTML = `<div class="cl-tools"><a href="#" data-a="all">all</a> · <a href="#" data-a="none">none</a></div>` + classes.map(c => `<label class="cl-item"><input type="checkbox" data-c="${esc(c)}" checked><span>${esc(c)}</span></label>`).join('');
    cls.querySelectorAll('input').forEach(cb => cb.addEventListener('change', () => { cb.checked ? state.classes.add(cb.dataset.c) : state.classes.delete(cb.dataset.c); userChange(); }));
    cls.querySelectorAll('a').forEach(a => a.addEventListener('click', ev => { ev.preventDefault(); state.classes.clear(); if (a.dataset.a === 'all') classes.forEach(c => state.classes.add(c)); syncClassList(); userChange(); }));
    const tr = document.getElementById('f-trend');
    if (KIND === 'reference') {
      const HAS_T0 = R.trend.some(t => t === 0);
      const TREND_TIP = { 0: 'significant in the query but the cell types are not in the reference panel: counted as interactions, no Up/Down call', 4: 'significant in neither dataset: listed in the table and the query-vs-reference scatter only, never counted as an interaction' };
      tr.innerHTML = (HAS_T0 ? [1, 2, 3, 0, 4] : [1, 2, 3, 4]).map(t => `<label class="cl-item" title="${esc(TREND_TIP[t] || '')}"><input type="checkbox" data-t="${t}" ${state.trends.has(t) ? 'checked' : ''}><span class="sw" style="background:${TREND_COLOR[t] || '#777'}"></span><span>${TREND[t]}</span></label>`).join('');
      tr.querySelectorAll('input').forEach(cb => cb.addEventListener('change', () => { const t = +cb.dataset.t; cb.checked ? state.trends.add(t) : state.trends.delete(t); userChange(); }));
    } else if (KIND === 'split') {
      document.getElementById('f-trend-title').textContent = 'Condition';
      tr.innerHTML = `<select id="f-cond">` + conditions.map((c, i) => `<option value="${i}">${esc(c)}</option>`).join('') + `</select>`;
      document.getElementById('f-cond').addEventListener('change', ev => { state.cond = +ev.target.value; userChange(); });
    } else tr.closest('.frow').style.display = 'none';
    const gene = document.getElementById('f-gene'); gene.addEventListener('input', () => { state.gene = gene.value.trim(); userChange(); });
    const fs = document.getElementById('f-signal');
    if (!SIG) { fs.style.display = 'none'; document.getElementById('f-signal-title').style.display = 'none'; }
    else {
      const all = SIGNAL_TYPES.concat(['complex-limited', 'one-sided']);
      fs.innerHTML = all.map(tname => `<label class="cl-item" title="${esc(SIGNAL_HELP[tname] || (tname === 'complex-limited' ? 'a complex whose Minimum-aggregated modules fail while Average-aggregated ones pass: one subunit limits it' : 'Arithmetic modules pass while Geometric ones fail: ligand and receptor levels are unbalanced'))}"><input type="checkbox" data-s="${tname}" checked>${SIGNAL_COLOR[tname] ? `<span class="sw" style="background:${SIGNAL_COLOR[tname]}"></span>` : ''}<span>${tname}</span></label>`).join('');
      fs.querySelectorAll('input').forEach(cb => cb.addEventListener('change', () => { cb.checked ? state.signals.add(cb.dataset.s) : state.signals.delete(cb.dataset.s); userChange(); }));
    }
    const sliders = [
      { id: 'f-mincs', set: v => { state.minCS = +v; }, get: () => state.minCS, fmt: v => (+v).toFixed(1) },
      { id: 'f-maxpairs', set: v => { state.maxPairs = +v >= DIST.pairsMax ? Infinity : +v; }, get: () => (state.maxPairs >= DIST.pairsMax ? DIST.pairsMax : state.maxPairs), fmt: v => (+v >= DIST.pairsMax ? 'any' : '≤ ' + v) },
      { id: 'f-minpct', set: v => { state.minPct = +v / 100; }, get: () => Math.round(state.minPct * 100), fmt: v => '≥ ' + v + '%' },
      { id: 'f-minmod', set: v => { state.minModules = +v; }, get: () => Math.min(state.minModules, nSel()), fmt: v => (+v === 0 ? 'off' : `≥ ${v} of ${nSel()}`) },
      { id: 'f-maxp', set: v => { state.maxP = Math.pow(10, -v); }, get: () => +(-Math.log10(state.maxP)).toFixed(1), fmt: v => '< 1e-' + v },
      { id: 'f-minz', set: v => { state.minZ = +v; }, get: () => state.minZ, fmt: v => '≥ ' + v },
      { id: 'f-minbs', set: v => { state.minBS = +v / 100; }, get: () => Math.round(state.minBS * 100), fmt: v => '≥ ' + v + '%' },
      { id: 'f-minrep', set: v => { state.minRep = +v / 100; }, get: () => Math.round(state.minRep * 100), fmt: v => '≥ ' + v + '%' },
    ];
    if (!HAS_Z) document.getElementById('f-z-row').style.display = 'none'; else document.getElementById('f-minz').max = Math.ceil(quantile(R.z, 0.99) || 10);
    if (!HAS_BS) document.getElementById('f-bs-row').style.display = 'none';
    if (!HAS_REP) document.getElementById('f-rep-row').style.display = 'none';
    document.getElementById('f-mincs').max = Math.ceil(DIST.csMax); document.getElementById('f-mincs').step = DIST.csMax > 20 ? 0.5 : DIST.csMax > 4 ? 0.1 : 0.05;
    document.getElementById('f-maxpairs').max = DIST.pairsMax;
    if (KIND === 'reference') document.getElementById('f-maxp').closest('.frow').style.display = 'none'; else document.getElementById('f-maxp').max = Math.max(2, DIST.nlogpMax);
    sliders.forEach(({ id, set, fmt }) => { const s = document.getElementById(id); const v = document.getElementById(id + '-v'); s.addEventListener('input', () => { set(s.value); v.textContent = fmt(s.value); userChange(); }); });
    function syncSliders() { sliders.forEach(({ id, get, fmt }) => { const s = document.getElementById(id); s.value = get(); document.getElementById(id + '-v').textContent = fmt(s.value); }); }
    // scoring modules as three factors: expression summary × complex aggregation × ligand–receptor combination
    function buildModuleFactors() {
      const host = document.getElementById('f-modfactors');
      const FACTORS = [
        ['modSummary', 'Expression summary', SUMMARIES.map(s => [s, sumLabel(s), sumHelp(s)])],
        ['modAgg', 'Complex aggregation', AGGS.map(a => [a, AGG_SHORT[a] || a, AGG_HELP[a] || ''])],
        ['modComb', 'Ligand–receptor combination', COMBS.map(c => [c, COMB_SHORT[c] || c, COMB_HELP[c] || ''])],
      ];
      host.innerHTML = FACTORS.map(([key, title, levels]) => `<div class="fac"><div class="fac-title">${title}</div>` +
        levels.map(([v, short, help]) => `<label class="cl-item fac-item"><input type="checkbox" data-f="${key}" data-v="${esc(v)}" ${state[key].has(v) ? 'checked' : ''}><b>${esc(short)}</b><span>${esc(help)}</span></label>`).join('') + '</div>').join('') +
        `<div class="fac-count muted" id="f-modcount"></div>`;
      host.querySelectorAll('input').forEach(cb => cb.addEventListener('change', () => { cb.checked ? state[cb.dataset.f].add(cb.dataset.v) : state[cb.dataset.f].delete(cb.dataset.v); modChanged(); }));
    }
    function modChanged() { state.modSel = deriveModSel(); state.modMask = maskOf(state.modSel); if (state.minModules === 0 && nSel() > 0) state.minModules = 1; if (state.minModules > nSel()) state.minModules = nSel(); syncControls(); userChange(); }
    if (modules) buildModuleFactors(); else { document.getElementById('f-mod-row').style.display = 'none'; document.getElementById('f-mod-title').style.display = 'none'; }
    syncControls = function () {
      document.querySelectorAll('#f-senders input').forEach(cb => { cb.checked = state.senders.has(+cb.dataset.i); });
      document.querySelectorAll('#f-receivers input').forEach(cb => { cb.checked = state.receivers.has(+cb.dataset.i); });
      document.querySelectorAll('#f-trend input').forEach(cb => { cb.checked = state.trends.has(+cb.dataset.t); });
      const cond = document.getElementById('f-cond'); if (cond) cond.value = state.cond;
      document.getElementById('f-gene').value = state.gene;
      syncClassList();
      document.querySelectorAll('#f-signal input').forEach(cb => { cb.checked = state.signals.has(cb.dataset.s); });
      document.querySelectorAll('#f-modfactors input').forEach(cb => { cb.checked = state[cb.dataset.f].has(cb.dataset.v); });
      const mc = document.getElementById('f-modcount'); if (mc) mc.textContent = !nSel() ? 'No module selected: module support is not applied.' : state.minModules === 0 ? `${nSel()} of ${nModules} modules selected; support filter off, the run's own call (${nModules}-module Cauchy combination) stands. Changing the selection or the slider keeps only interactions that at least that many of the selected modules call significant on their own.` : `${nSel()} of ${nModules} modules selected. Keeps the interactions that at least ${state.minModules} of them call significant on their own; a stricter filter within the set called by the ${nModules}-module Cauchy combination (0 = off). To test another combination, choose the modules when running FastCCC (--summary, --aggregation, --combination).`;
      document.getElementById('f-minmod').max = Math.max(1, nSel());
      syncSliders();
    };
    ['default', 'highlights', 'custom', 'reset'].forEach(n => document.getElementById('f-' + n).addEventListener('click', ev => { ev.preventDefault(); setPreset(n); }));
    document.getElementById('dot-top').addEventListener('change', ev => { state.dotTop = +ev.target.value; renderDotplot(); });
    document.getElementById('dot-cols').addEventListener('change', ev => { state.dotCols = ev.target.value === 'auto' ? 'auto' : +ev.target.value; renderDotplot(); });
    const cb = document.getElementById('chord-by');
    if (KIND !== 'reference') cb.closest('label').style.display = 'none'; else cb.addEventListener('change', ev => { state.chordBy = ev.target.value; renderChord(); });
    document.getElementById('roles-by').addEventListener('change', ev => { state.rolesBy = ev.target.value; renderRoles(); });
    document.getElementById('pat-dir').addEventListener('change', ev => { state.patDir = ev.target.value; renderPatterns(); });
    document.getElementById('pat-k').addEventListener('change', ev => { state.patK = ev.target.value; renderPatterns(); });
    const pb = document.getElementById('path-by');
    if (P.genesets) P.genesets.libraries.forEach((lib, i) => pb.insertAdjacentHTML('beforeend', `<option value="${i}">${esc(lib.name.replace(/_/g, ' '))}</option>`));
    if (!DB_HAS_CLASSES) { pb.querySelector('option[value="cpdb"]').remove(); if (P.genesets) state.pathBy = '0'; }
    pb.value = state.pathBy; pb.addEventListener('change', ev => { state.pathBy = ev.target.value; renderPathways(); });
    document.getElementById('path-rows').value = state.pathRows;
    document.getElementById('path-rows').addEventListener('change', ev => { state.pathRows = ev.target.value; renderPathways(); });
    const sb = document.getElementById('scatter-by'); if (KIND !== 'reference') sb.closest('label').style.display = 'none'; else sb.addEventListener('change', ev => { state.scatterBy = ev.target.value; renderScatter(); });
    document.getElementById('path-rank').addEventListener('change', ev => { state.pathRank = ev.target.value; renderPathways(); });
    document.getElementById('path-top').addEventListener('change', ev => { state.pathTop = +ev.target.value; renderPathways(); });
    if (KIND !== 'reference') document.getElementById('ref-row').style.display = 'none';
    if (KIND !== 'split') document.getElementById('diff-card').style.display = 'none';
    else {
      const a = document.getElementById('diff-a'), b = document.getElementById('diff-b');
      conditions.forEach((c, i) => { a.insertAdjacentHTML('beforeend', `<option value="${i}">${esc(c)}</option>`); b.insertAdjacentHTML('beforeend', `<option value="${i}">${esc(c)}</option>`); });
      state.diffA = conditions.length > 1 ? 1 : 0; state.diffB = conditions.length > 2 ? 2 : 0; a.value = state.diffA; b.value = state.diffB;
      a.addEventListener('change', ev => { state.diffA = +ev.target.value; state.diffPair = null; renderDiff(); });
      b.addEventListener('change', ev => { state.diffB = +ev.target.value; state.diffPair = null; renderDiff(); });
    }
    if (!P.patterns) document.getElementById('patterns').closest('.card').style.display = 'none';
    if (!P.genesets) document.getElementById('enrichment').closest('.card').style.display = 'none';
    else { const sel = document.getElementById('enr-lib'); P.genesets.libraries.forEach((l, i) => sel.insertAdjacentHTML('beforeend', `<option value="${i}">${esc(l.name.replace(/_/g, ' '))}</option>`)); sel.addEventListener('change', ev => { state.enrLib = +ev.target.value; renderEnrichment(); }); }
    if (!P.cells) document.getElementById('cells-card').style.display = 'none';
    else { const sel = document.getElementById('cells-gene'); P.cells.genes.forEach((g, i) => sel.insertAdjacentHTML('beforeend', `<option value="${i}">${esc(g)}</option>`)); sel.addEventListener('change', ev => { state.cellsGene = +ev.target.value; renderCells(); }); }
  }

  // ── sender × receiver matrices ─────────────────────────────────────────────
  function pairMatrix(indices) {
    const n = CT.length; const mk = () => Array.from({ length: n }, () => new Array(n).fill(0));
    const up = mk(), down = mk(), all = mk();
    indices.forEach(i => { all[R.sender[i]][R.receiver[i]]++; if (KIND === 'reference') { if (R.trend[i] === 1) up[R.sender[i]][R.receiver[i]]++; if (R.trend[i] === 2) down[R.sender[i]][R.receiver[i]]++; } });
    return { up, down, all };
  }
  function matrixHeatmap(host, order, valueOf, countOf, opts) {
    const rows = order, cols = opts.cols || order;  // row order = senders, column order = receivers
    const maxLab = d3.max(rows.concat(cols), d => CT[d].name.length) || 8;
    const avail = Math.max(360, (host.node() && host.node().clientWidth) || 640);  // fit the card width; the fit step adds room for overhanging labels
    const labelW = Math.min(260, 30 + maxLab * 10.5 * 0.62);
    const cell = Math.max(6, Math.min(30, Math.floor((avail - labelW - 24) / Math.max(cols.length, 1)), Math.floor(900 / Math.max(rows.length, 1))));
    const labFont = Math.max(7, Math.min(10.5, cell - 2));
    const left = Math.min(260, 30 + maxLab * labFont * 0.62), top = Math.min(220, 40 + maxLab * labFont * 0.6);
    const W = left + cell * cols.length + 20, H = top + cell * rows.length + 44;
    const svg = newSvg(host, W, H);
    const vals = []; rows.forEach(s => cols.forEach(r => vals.push(valueOf(s, r))));
    const ext = d3.max(vals.map(Math.abs)) || 1;
    const color = opts.diverging ? d3.scaleDiverging(t => d3.interpolateRgbBasis(['#3b6fb6', '#f4f2ee', '#d1495b'])(t)).domain([-ext, 0, ext])
      : d3.scaleSequential(d3.interpolateRgbBasis(['#f6f4ef', '#9fb3d1', '#2f4f7f'])).domain([0, ext]);
    const g = svg.append('g').attr('transform', `translate(${left},${top})`);
    rows.forEach((s, yi) => cols.forEach((r, xi) => {
      const v = valueOf(s, r), isSel = opts.selected && opts.selected[0] === s && opts.selected[1] === r;
      g.append('rect').attr('x', xi * cell).attr('y', yi * cell).attr('width', cell - 1).attr('height', cell - 1).attr('rx', 2)
        .attr('fill', countOf(s, r) === 0 ? '#faf9f6' : color(v)).attr('stroke', isSel ? '#111' : 'none').attr('stroke-width', 1.5).style('cursor', 'pointer')
        .on('mousemove', ev => showTip(ev, opts.tip(s, r))).on('mouseleave', hideTip).on('click', () => opts.click(s, r, isSel));
    }));
    g.selectAll('.rl').data(rows).enter().append('text').attr('class', 'axis-lab').style('font-size', labFont + 'px').attr('x', -8).attr('y', (d, i) => i * cell + cell / 2 + labFont * 0.35).attr('text-anchor', 'end').text(d => CT[d].name);
    g.selectAll('.cl').data(cols).enter().append('text').attr('class', 'axis-lab').style('font-size', labFont + 'px').attr('transform', (d, i) => `translate(${i * cell + cell / 2 + labFont * 0.35},-8) rotate(-60)`).attr('text-anchor', 'start').text(d => CT[d].name);
    g.selectAll('.rs').data(rows).enter().append('rect').attr('x', -5).attr('y', (d, i) => i * cell).attr('width', 3).attr('height', cell - 1).attr('fill', d => CT[d].color);
    g.selectAll('.cs').data(cols).enter().append('rect').attr('x', (d, i) => i * cell).attr('y', -5).attr('width', cell - 1).attr('height', 3).attr('fill', d => CT[d].color);
    svg.append('text').attr('x', left - 8).attr('y', top - 70).attr('text-anchor', 'end').attr('class', 'axis-title').text('sender ↓');
    svg.append('text').attr('x', left).attr('y', 14).attr('class', 'axis-title').text('receiver →');
    const lg = svg.append('g').attr('transform', `translate(${left + 90},${H - 16})`); const lw = 120;
    const gid = 'grad' + Math.random().toString(36).slice(2, 8);
    const grad = svg.append('defs').append('linearGradient').attr('id', gid);
    d3.range(0, 1.01, 0.1).forEach(t => grad.append('stop').attr('offset', t).attr('stop-color', opts.diverging ? color(-ext + 2 * ext * t) : color(ext * t)));
    lg.append('rect').attr('x', 0).attr('y', -8).attr('width', lw).attr('height', 6).attr('fill', `url(#${gid})`);
    lg.append('text').attr('x', -6).attr('y', -3).attr('text-anchor', 'end').attr('class', 'axis-lab').text(opts.diverging ? `${ext} more ${opts.neg}` : '0');
    lg.append('text').attr('x', lw + 6).attr('y', -3).attr('class', 'axis-lab').text(opts.diverging ? `${ext} more ${opts.pos}` : `${ext} interactions`);
    if (opts.diverging) lg.append('text').attr('x', lw / 2).attr('y', 10).attr('text-anchor', 'middle').attr('class', 'axis-lab').text(opts.legend);
    return svg;
  }
  function renderHeatmap() {
    const host = d3.select('#heatmap'); host.selectAll('*').remove();
    const M = pairMatrix(filtered);
    // rows = ticked senders, columns = ticked receivers; empty cells stay visible as grey
    const rows = CT.map((_, i) => i).filter(i => state.senders.has(i)), cols = CT.map((_, i) => i).filter(i => state.receivers.has(i));
    if (!rows.length || !cols.length) { host.append('div').attr('class', 'empty').text('Tick at least one sender and one receiver.'); return; }
    const diverging = KIND === 'reference';
    matrixHeatmap(host, rows, (s, r) => diverging ? M.up[s][r] - M.down[s][r] : M.all[s][r], (s, r) => M.all[s][r], {
      cols, diverging, neg: 'Down', pos: 'Up', legend: 'Up − Down per pair', selected: state.pair,
      tip: (s, r) => `<b>${esc(CT[s].name)}</b> → <b>${esc(CT[r].name)}</b><br>` + (diverging ? `Up ${M.up[s][r]} · Down ${M.down[s][r]} · total ${M.all[s][r]}` : `${M.all[s][r]} interactions`) + '<br><i>click to focus this pair</i>',
      click: (s, r, isSel) => { state.pair = isSel ? null : [s, r]; update(); },
    });
  }

  // ── chord diagram ──────────────────────────────────────────────────────────
  function renderChord() {
    const host = d3.select('#chord'); host.selectAll('*').remove();
    const M = pairMatrix(filtered);
    let active = CT.map((_, i) => i).filter(i => M.all[i].some(v => v) || M.all.some(r => r[i]));
    if (active.length < 2) { host.append('div').attr('class', 'empty').text('Not enough cell-type pairs in the current selection.'); return; }
    const byTrend = KIND === 'reference' && state.chordBy === 'trend';
    const MAXNODES = 36; let folded = 0;
    if (active.length > MAXNODES) {
      const degree = i => d3.sum(M.all[i]) + d3.sum(M.all.map(r => r[i]));
      const keep = active.slice().sort((a, b) => degree(b) - degree(a)).slice(0, MAXNODES - 1);
      const rest = active.filter(i => !keep.includes(i)); folded = rest.length;
      const fold = mat => { const out = keep.map(a => keep.map(b => mat[a][b])); keep.forEach((a, ai) => { out[ai].push(d3.sum(rest.map(b => mat[a][b]))); }); out.push(keep.map(b => d3.sum(rest.map(a => mat[a][b])))); out[keep.length].push(d3.sum(rest.map(a => d3.sum(rest.map(b => mat[a][b]))))); return out; };
      M.allF = fold(M.all); M.upF = fold(M.up); M.downF = fold(M.down);
      active = keep.concat([CT.length]);
    }
    const nameOf = i => i === CT.length ? `other (${folded})` : CT[i].name, colorOf = i => i === CT.length ? '#bdbab2' : CT[i].color;
    const allM = M.allF ? M.allF : active.map(s => active.map(r => M.all[s][r])), upM = M.upF ? M.upF : active.map(s => active.map(r => M.up[s][r])), downM = M.downF ? M.downF : active.map(s => active.map(r => M.down[s][r]));
    const maxLab = d3.max(active, i => nameOf(i).length) || 6;
    const size = active.length > 20 ? 680 : 560, outer = size / 2 - Math.min(150, 24 + maxLab * 6.2), inner = outer - 12;
    const svg = newSvg(host, size, size).attr('width', '100%').style('max-width', `${size}px`).style('display', 'block').style('margin', '0 auto');
    const g = svg.append('g').attr('transform', `translate(${size / 2},${size / 2})`);
    const mat = active.map((s, si) => active.map((r, ri) => byTrend ? upM[si][ri] + downM[si][ri] : allM[si][ri]));
    if (!mat.some(row => row.some(v => v))) { host.selectAll('*').remove(); host.append('div').attr('class', 'empty').text('No Up or Down interactions in the current selection.'); return; }
    const chord = d3.chordDirected().padAngle(0.03).sortSubgroups(d3.descending)(mat);
    const arc = d3.arc().innerRadius(inner).outerRadius(outer);
    const ribbon = d3.ribbonArrow().radius(inner - 2).padAngle(0.005);
    g.append('g').selectAll('path').data(chord.groups).enter().append('path').attr('d', arc).attr('fill', d => colorOf(active[d.index]))
      .on('mousemove', (ev, d) => showTip(ev, `<b>${esc(nameOf(active[d.index]))}</b><br>sends ${d3.sum(mat[d.index])} · receives ${d3.sum(mat.map(r => r[d.index]))}`)).on('mouseleave', hideTip);
    g.append('g').selectAll('text').data(chord.groups).enter().append('text').each(d => { d.angle = (d.startAngle + d.endAngle) / 2; })
      .attr('dy', '0.35em').attr('class', 'chord-lab').attr('transform', d => `rotate(${d.angle * 180 / Math.PI - 90}) translate(${outer + 6}) ${d.angle > Math.PI ? 'rotate(180)' : ''}`)
      .attr('text-anchor', d => d.angle > Math.PI ? 'end' : null).style('font-size', (active.length > 30 ? 9 : 10.5) + 'px').text(d => nameOf(active[d.index]));
    const rib = g.append('g').attr('fill-opacity', 0.72);
    chord.forEach(c => {
      const s = active[c.source.index], r = active[c.target.index]; const isVirtual = s === CT.length || r === CT.length;
      if (byTrend) {
        const up = upM[c.source.index][c.target.index], down = downM[c.source.index][c.target.index]; if (!up && !down) return;
        const tot = up + down; const parts = []; if (up) parts.push([1, up / tot]); if (down) parts.push([2, down / tot]);
        let a0 = c.source.startAngle; const span = c.source.endAngle - c.source.startAngle;
        parts.forEach(([t, f]) => {
          const seg = { source: { startAngle: a0, endAngle: a0 + span * f }, target: c.target };
          rib.append('path').attr('d', ribbon(seg)).attr('fill', TREND_COLOR[t]).attr('stroke', 'none')
            .on('mousemove', ev => showTip(ev, `<b>${esc(nameOf(s))}</b> → <b>${esc(nameOf(r))}</b><br>Up ${up} · Down ${down}`)).on('mouseleave', hideTip)
            .on('click', () => { if (!isVirtual) { state.pair = [s, r]; update(); } });
          a0 += span * f;
        });
      } else {
        rib.append('path').attr('d', ribbon(c)).attr('fill', colorOf(s)).attr('stroke', 'none')
          .on('mousemove', ev => showTip(ev, `<b>${esc(nameOf(s))}</b> → <b>${esc(nameOf(r))}</b><br>${allM[c.source.index][c.target.index]} interactions`)).on('mouseleave', hideTip)
          .on('click', () => { if (!isVirtual) { state.pair = [s, r]; update(); } });
      }
    });
  }

  // ── dot plot ───────────────────────────────────────────────────────────────
  function topLRs(idx, n) {
    const byLR = new Map();
    idx.forEach(i => { const k = R.lr[i]; if (!byLR.has(k)) byLR.set(k, { n: 0, cs: 0 }); const o = byLR.get(k); o.n++; o.cs = Math.max(o.cs, R.cs[i] || 0); });
    return Array.from(byLR.entries()).sort((a, b) => b[1].n - a[1].n || b[1].cs - a[1].cs).slice(0, n).map(e => e[0]);
  }
  function renderDotplot() {
    const host = d3.select('#dotplot'); host.selectAll('*').remove();
    const idx = filtered;
    if (!idx.length) { host.append('div').attr('class', 'empty').text('No interactions match the current filters.'); return; }
    const lrs = topLRs(idx, state.dotTop); const lrPos = new Map(lrs.map((k, i) => [k, i]));
    const rows = idx.filter(i => lrPos.has(R.lr[i]));
    const pairCount = new Map(); rows.forEach(i => { const k = R.sender[i] * 1000 + R.receiver[i]; pairCount.set(k, (pairCount.get(k) || 0) + 1); });
    const colCap = state.dotCols === 'auto' ? Math.max(8, Math.floor((((host.node() && host.node().clientWidth) || 1100) - 230) / 18)) : state.dotCols;  // 'auto': as many columns as fit the card
    const keep = new Set(Array.from(pairCount.entries()).sort((a, b) => b[1] - a[1]).slice(0, colCap).map(e => e[0]));
    const shown = rows.filter(i => keep.has(R.sender[i] * 1000 + R.receiver[i]));
    const senders = Array.from(new Set(shown.map(i => R.sender[i]))).sort((a, b) => a - b);
    const cols = []; senders.forEach(s => { const recs = Array.from(new Set(shown.filter(i => R.sender[i] === s).map(i => R.receiver[i]))).sort((a, b) => a - b); recs.forEach(r => cols.push([s, r])); });
    const colPos = new Map(cols.map((c, i) => [c[0] * 1000 + c[1], i])); const truncated = pairCount.size - cols.length;
    const cw = 18, rh = 18, left = 200, top = 56;
    const recExtent = Math.ceil((d3.max(cols, c => CT[c[1]].name.length) || 4) * 6.2 * Math.sin(Math.PI / 3)) + 14;  // rotated receiver labels
    const groupW = s => cols.filter(c => c[0] === s).length * cw;
    const narrow = senders.some(s => CT[s].name.length * 6.6 > groupW(s) + 10);  // sender names wider than their column group: rotate them too
    const sendExtent = narrow ? Math.ceil((d3.max(senders, s => CT[s].name.length) || 4) * 6.6 * Math.sin(Math.PI / 3)) + 16 : 26;
    const barY = rh * lrs.length + recExtent + 4, bottom = recExtent + sendExtent + 24;
    const W = Math.max(left + cw * cols.length + 30, left + 640), H = top + rh * lrs.length + bottom;
    const svg = newSvg(host, W, H, 'dotplot-svg');
    const g = svg.append('g').attr('transform', `translate(${left},${top})`);
    lrs.forEach((k, yi) => { const t = g.append('text').attr('x', -8).attr('y', yi * rh + rh / 2 + 4).attr('text-anchor', 'end').style('font-style', 'italic'); t.append('tspan').attr('fill', LIG_INK).text(LR[k].ligand); t.append('tspan').attr('fill', '#888').text(' – '); t.append('tspan').attr('fill', REC_INK).text(LR[k].receptor); });
    lrs.forEach((_, yi) => g.append('line').attr('x1', 0).attr('x2', cw * cols.length).attr('y1', yi * rh + rh / 2).attr('y2', yi * rh + rh / 2).attr('stroke', '#efede8'));
    let x0 = 0;
    senders.forEach(s => {
      const nrec = cols.filter(c => c[0] === s).length; const x1 = x0 + nrec * cw;
      if (x0 > 0) g.append('line').attr('x1', x0 - 0.5).attr('x2', x0 - 0.5).attr('y1', -8).attr('y2', rh * lrs.length + 4).attr('stroke', '#c9c7c0').attr('stroke-dasharray', '3,3');
      g.append('line').attr('x1', x0 + 2).attr('x2', x1 - 2).attr('y1', barY).attr('y2', barY).attr('stroke', CT[s].color).attr('stroke-width', 2);
      if (narrow) g.append('text').attr('transform', `translate(${(x0 + x1) / 2 + 4},${barY + 8}) rotate(60)`).attr('fill', CT[s].color).style('font-weight', 600).text(CT[s].name);
      else g.append('text').attr('x', (x0 + x1) / 2).attr('y', barY + 15).attr('text-anchor', 'middle').attr('fill', CT[s].color).style('font-weight', 600).text(CT[s].name);
      x0 = x1;
    });
    cols.forEach(([s, r], xi) => g.append('text').attr('transform', `translate(${xi * cw + cw / 2 + 4},${rh * lrs.length + 8}) rotate(60)`).attr('fill', '#444').text(CT[r].name));
    const sizeVal = i => KIND === 'reference' ? (R.cs[i] || 0) : nlogp(R.p[i]);
    const smax = d3.max(shown, sizeVal) || 1;
    const rad = d3.scaleSqrt().domain([0, smax]).range([1.5, cw / 2 - 2]);
    const fillOf = i => KIND === 'reference' ? TREND_COLOR[R.trend[i]] : KIND === 'split' ? COND_COLORS[R.cond[i] % COND_COLORS.length] : SINGLE;
    g.selectAll('circle').data(shown).enter().append('circle')
      .attr('cx', i => colPos.get(R.sender[i] * 1000 + R.receiver[i]) * cw + cw / 2).attr('cy', i => lrPos.get(R.lr[i]) * rh + rh / 2)
      .attr('r', i => rad(sizeVal(i))).attr('fill', fillOf).attr('fill-opacity', 0.9).attr('stroke', i => state.pinned.has(i) ? '#111' : 'white').attr('stroke-width', 1)
      .style('cursor', 'pointer').on('mousemove', (ev, i) => showTip(ev, rowTip(i))).on('mouseleave', hideTip).on('click', (ev, i) => selectRow(i));
    svg.append('text').attr('x', left).attr('y', 16).attr('class', 'plot-title').text(`Top ${lrs.length} ligand–receptor pairs by number of cell-type pairs · dot size = ${KIND === 'reference' ? 'communication score' : '−log10 p'}` + (truncated > 0 ? ` · ${cols.length} of ${pairCount.size} cell-type pairs shown` : ''));
    const lg = svg.append('g').attr('transform', `translate(${left},34)`); let lx = 0;
    if (KIND === 'reference') [1, 2, 3, 4].filter(t => state.trends.has(t)).forEach(t => { lg.append('circle').attr('cx', lx).attr('cy', 0).attr('r', 5).attr('fill', TREND_COLOR[t]); lg.append('text').attr('x', lx + 9).attr('y', 4).text(TREND[t]); lx += 90; });
    else if (KIND === 'split') { lg.append('circle').attr('cx', 0).attr('cy', 0).attr('r', 5).attr('fill', COND_COLORS[state.cond % COND_COLORS.length]); lg.append('text').attr('x', 9).attr('y', 4).text(conditions[state.cond]); lx = 120; }
    [0.25, 0.5, 1].forEach((f, i) => { lg.append('circle').attr('cx', lx + 40 + i * 44).attr('cy', 0).attr('r', rad(smax * f)).attr('fill', 'none').attr('stroke', '#666'); lg.append('text').attr('x', lx + 40 + i * 44 + 10).attr('y', 4).text((smax * f).toFixed(1)); });
  }

  // ── signalling roles (outgoing vs incoming) ────────────────────────────────
  function renderRoles() {
    const host = d3.select('#roles'); host.selectAll('*').remove();
    if (!filtered.length) { host.append('div').attr('class', 'empty').text('Nothing selected.'); return; }
    const out = new Array(CT.length).fill(0), inc = new Array(CT.length).fill(0);
    filtered.forEach(i => { const w = state.rolesBy === 'strength' ? (R.cs[i] || 0) : 1; out[R.sender[i]] += w; inc[R.receiver[i]] += w; });
    const pts = CT.map((c, i) => ({ i, name: c.name, color: c.color, x: out[i], y: inc[i] })).filter(d => d.x || d.y);
    const W = 520, H = 400, m = { l: 56, r: 20, t: 20, b: 46 };
    const svg = newSvg(host, W, H);
    const x = d3.scaleLinear().domain([0, (d3.max(pts, d => d.x) || 1) * 1.1]).range([m.l, W - m.r]).nice();
    const y = d3.scaleLinear().domain([0, (d3.max(pts, d => d.y) || 1) * 1.1]).range([H - m.b, m.t]).nice();
    svg.append('g').attr('transform', `translate(0,${H - m.b})`).call(d3.axisBottom(x).ticks(6)).call(axisStyle);
    svg.append('g').attr('transform', `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(6)).call(axisStyle);
    svg.append('text').attr('x', (m.l + W - m.r) / 2).attr('y', H - 8).attr('text-anchor', 'middle').attr('class', 'axis-lab').text(`outgoing (${state.rolesBy === 'strength' ? 'summed score' : 'interactions'} as sender)`);
    svg.append('text').attr('transform', `translate(14,${(m.t + H - m.b) / 2}) rotate(-90)`).attr('text-anchor', 'middle').attr('class', 'axis-lab').text('incoming (as receiver)');
    const lim = Math.min(x.domain()[1], y.domain()[1]);
    svg.append('line').attr('x1', x(0)).attr('y1', y(0)).attr('x2', x(lim)).attr('y2', y(lim)).attr('stroke', '#e6e4de').attr('stroke-dasharray', '3,3');
    const size = d3.scaleSqrt().domain([0, d3.max(CT, c => c.n_cells) || 1]).range([4, 14]);
    svg.selectAll('circle').data(pts).enter().append('circle').attr('cx', d => x(d.x)).attr('cy', d => y(d.y)).attr('r', d => size(CT[d.i].n_cells || 0)).attr('fill', d => d.color).attr('fill-opacity', 0.85).attr('stroke', 'white').style('cursor', 'pointer')
      .on('mousemove', (ev, d) => showTip(ev, `<b>${esc(d.name)}</b><br>outgoing ${fmt2(d.x)} · incoming ${fmt2(d.y)}<br>${(CT[d.i].n_cells || 0).toLocaleString()} cells<br><i>click to show this cell type as sender</i>`)).on('mouseleave', hideTip)
      .on('click', (ev, d) => { assign(state.senders, [d.i]); assign(state.receivers, CT.map((_, k) => k)); syncControls(); userChange(); });
    const MAXLAB = 16; const labelled = pts.length > MAXLAB ? new Set(pts.slice().sort((a, b) => (b.x + b.y) - (a.x + a.y)).slice(0, MAXLAB).map(d => d.i)) : null;
    const labs = pts.filter(d => !labelled || labelled.has(d.i)).map(d => ({ d, lx: x(d.x) + size(CT[d.i].n_cells || 0) + 3, ly: y(d.y) + 4 })).sort((a, b) => a.ly - b.ly || a.lx - b.lx);
    if (labelled) svg.append('text').attr('x', W - m.r).attr('y', m.t + 18).attr('text-anchor', 'end').attr('class', 'axis-lab').attr('fill', '#999').text(`labels: ${MAXLAB} busiest cell types · hover for the rest`);
    for (let a = 1; a < labs.length; a++) { for (let b = 0; b < a; b++) { if (Math.abs(labs[a].lx - labs[b].lx) < 90 && Math.abs(labs[a].ly - labs[b].ly) < 11) labs[a].ly = labs[b].ly + 11; } }
    svg.selectAll('.lab').data(labs).enter().append('text').attr('x', l => l.lx).attr('y', l => l.ly).attr('class', 'axis-lab').text(l => l.d.name);
    svg.append('text').attr('x', W - m.r).attr('y', m.t + 4).attr('text-anchor', 'end').attr('class', 'axis-lab').attr('fill', '#999').text('dot size = number of cells');
  }

  // ── sankey: sender → LR → receiver ─────────────────────────────────────────
  function sankeyDraw(host, nodes, links, W, H, labelStyle) {
    // the drawing fills the card; side margins are sized from the longest left ('s') and right ('r') labels
    W = Math.max(420, Math.floor((host.node() && host.node().clientWidth) || W));
    const lw = kind => d3.max(nodes.filter(n => n.kind === kind), n => n.label.length) || 0;
    const leftPad = Math.min(220, 12 + lw('s') * 6.2), rightPad = Math.min(240, 12 + lw('r') * 5.6);
    const svg = newSvg(host, W, H);
    const sk = d3.sankey().nodeId(d => d.index).nodeWidth(10).nodePadding(6).extent([[leftPad, 14], [W - rightPad, H - 12]]).nodeAlign(d3.sankeyJustify);
    const graph = sk({ nodes: nodes.map((n, index) => ({ ...n, index })), links: links.map(l => ({ ...l })) });
    svg.append('g').attr('fill', 'none').selectAll('path').data(graph.links).enter().append('path').attr('d', d3.sankeyLinkHorizontal()).attr('stroke', d => d.color).attr('stroke-opacity', 0.42).attr('stroke-width', d => Math.max(1, d.width))
      .on('mousemove', (ev, d) => showTip(ev, `${esc(d.source.label)} → ${esc(d.target.label)}<br>${d.tip || d.value}`)).on('mouseleave', hideTip);
    svg.append('g').selectAll('rect').data(graph.nodes).enter().append('rect').attr('x', d => d.x0).attr('y', d => d.y0).attr('width', d => d.x1 - d.x0).attr('height', d => Math.max(1, d.y1 - d.y0)).attr('fill', d => d.color).attr('rx', 2);
    svg.append('g').selectAll('text').data(graph.nodes).enter().append('text').attr('x', d => d.kind === 'r' ? d.x1 + 5 : d.kind === 's' ? d.x0 - 5 : (d.x0 + d.x1) / 2).attr('y', d => (d.y0 + d.y1) / 2 + 3.5)
      .attr('text-anchor', d => d.kind === 'r' ? 'start' : d.kind === 's' ? 'end' : 'middle').attr('class', 'axis-lab').style('font-style', d => d.kind === 'l' ? 'italic' : null).style('font-size', d => d.kind === 'l' ? '9.5px' : '10.5px').style('font-weight', d => d.kind === 'p' ? 600 : 400)
      .attr('paint-order', 'stroke').attr('stroke', d => (d.kind === 'l' || d.kind === 'p') ? 'white' : 'none').attr('stroke-width', 3).text(d => d.label);
    return svg;
  }
  function renderSankey() {
    const host = d3.select('#sankey'); host.selectAll('*').remove();
    if (!filtered.length || !d3.sankey) { host.append('div').attr('class', 'empty').text('Nothing selected.'); return; }
    const lrs = topLRs(filtered, 20); const lrSet = new Set(lrs);
    const rows = filtered.filter(i => lrSet.has(R.lr[i]));
    const nodes = [], nodeId = new Map();
    const node = (key, label, color, kind) => { if (!nodeId.has(key)) { nodeId.set(key, nodes.length); nodes.push({ key, label, color, kind }); } return nodeId.get(key); };
    const links = new Map();
    const add = (a, b, color) => { const k = a + '>' + b; if (!links.has(k)) links.set(k, { source: a, target: b, value: 0, color }); links.get(k).value += 1; };
    rows.forEach(i => {
      const s = node('s' + R.sender[i], CT[R.sender[i]].name, CT[R.sender[i]].color, 's');
      const l = node('l' + R.lr[i], `${LR[R.lr[i]].ligand} – ${LR[R.lr[i]].receptor}`, '#c9c7c0', 'l');
      const r = node('r' + R.receiver[i], CT[R.receiver[i]].name, CT[R.receiver[i]].color, 'r');
      add(s, l, CT[R.sender[i]].color); add(l, r, CT[R.receiver[i]].color);
    });
    const ls = Array.from(links.values()).map(l => ({ ...l, tip: `${l.value} interaction${l.value === 1 ? '' : 's'}` }));
    const W = 620, H = Math.max(260, 16 * Math.max(nodes.filter(n => n.kind === 'l').length, nodes.filter(n => n.kind === 's').length, nodes.filter(n => n.kind === 'r').length) + 40);
    const svg = sankeyDraw(host, nodes, ls, W, H);
    svg.append('text').attr('x', 140).attr('y', 10).attr('class', 'axis-lab').attr('fill', '#777').text(`top ${lrs.length} LR pairs of the selection · ribbon width = number of interactions`);
  }

  // ── pathway families × cell types ──────────────────────────────────────────
  function renderPathways() {
    const host = d3.select('#pathways'); host.selectAll('*').remove();
    const focus = document.getElementById('path-focus');
    focus.textContent = state.setFilter ? `filtered to ${state.setFilter.name} (click the row again to clear)` : '';
    if (!filtered.length) { host.append('div').attr('class', 'empty').text('Nothing selected.'); return; }
    const byLib = state.pathBy !== 'cpdb' && P.genesets && P.genesets.libraries[+state.pathBy];
    if (!byLib && !DB_HAS_CLASSES) { host.append('div').attr('class', 'empty').text('This database carries no pathway classification' + (P.genesets ? '; choose a gene-set library above.' : (DB_SPECIES !== 'human' ? ` and the gene-set libraries (Hallmark, KEGG, GO) are human-only, so none apply to a ${DB_SPECIES} database.` : ' and no gene sets were built (run without --no-genesets).'))); return; }
    const libIdx = byLib ? +state.pathBy : null, sets = byLib ? lrSets(libIdx) : null, setName = si => P.genesets.libraries[libIdx].sets[si].name;
    const famsOfLR = k => (byLib ? sets[k].map(setName) : [classOf(LR[k])]);
    // family sizes in the result's LR universe (for enrichment: hypergeometric over distinct LR pairs)
    const famK = new Map(); LR.forEach((_, k) => famsOfLR(k).forEach(c => famK.set(c, (famK.get(c) || 0) + 1)));
    const NN = LR.length;
    // per family: interaction counts per (role, cell type) column, and distinct LR pairs overall and per column
    const cols = []; CT.forEach((_, ct) => cols.push(['s', ct])); CT.forEach((_, ct) => cols.push(['r', ct]));
    const colKey = (role, ct) => role + ct;
    const cnt = new Map(), pairs = new Map(), colPairs = new Map(), colAll = new Map(); const selPairs = new Set();
    filtered.forEach(i => {
      const k = R.lr[i]; selPairs.add(k);
      [['s', R.sender[i]], ['r', R.receiver[i]]].forEach(([role, ct]) => { const ck = colKey(role, ct); if (!colAll.has(ck)) colAll.set(ck, new Set()); colAll.get(ck).add(k); });
      famsOfLR(k).forEach(c => {
        if (!pairs.has(c)) { pairs.set(c, new Set()); cnt.set(c, new Map()); colPairs.set(c, new Map()); }
        pairs.get(c).add(k);
        [['s', R.sender[i]], ['r', R.receiver[i]]].forEach(([role, ct]) => { const ck = colKey(role, ct); cnt.get(c).set(ck, (cnt.get(c).get(ck) || 0) + 1); if (!colPairs.get(c).has(ck)) colPairs.get(c).set(ck, new Set()); colPairs.get(c).get(ck).add(k); });
      });
    });
    if (!pairs.size) { host.append('div').attr('class', 'empty').text('None of the selected ligand–receptor pairs is in a gene set of this library.'); return; }
    const fams = Array.from(pairs.keys());
    const total = c => cnt.has(c) ? d3.sum(Array.from(cnt.get(c).values())) / 2 : 0;  // every interaction was counted once as sender and once as receiver; absent families count 0
    const enrichP = (k, n, K) => (k > 0 ? hyperSf(k, NN, K, n) : 1);
    // overall score per family: interactions, or enrichment q (BH over the families present)
    let score; const wholeUniverse = selPairs.size >= NN;
    if (state.pathRank === 'enrich' && !wholeUniverse) {
      const ps = fams.map(c => [c, enrichP(pairs.get(c).size, selPairs.size, famK.get(c) || pairs.get(c).size)]).sort((a, b) => a[1] - b[1]);
      const m = ps.length; const q = new Map(); let prev = 1; for (let i = m - 1; i >= 0; i--) { prev = Math.min(prev, ps[i][1] * m / (i + 1)); q.set(ps[i][0], prev); }
      score = c => -q.get(c);  // higher = more enriched
      renderPathways.q = q;
    } else { score = total; renderPathways.q = null; }
    const activeCols = cols.filter(([role, ct]) => colAll.has(colKey(role, ct)));
    const ctsActive = Array.from(new Set(activeCols.map(([, ct]) => ct))).sort((a, b) => a - b);
    const percell = state.pathRows === 'percell'; document.getElementById('path-rank').disabled = percell;
    // ---- cell-type specific view: for every (family, cell type, role) test whether the family is over-represented among the
    // cell type's interactions (hypergeometric; k = the cell type's interactions in the family, n = its interactions, K = the
    // family's interactions in that role, N = all selected interactions); shown as log2 enrichment, rows and columns clustered
    const colN = new Map(); filtered.forEach(i => { ['s' + R.sender[i], 'r' + R.receiver[i]].forEach(ck => colN.set(ck, (colN.get(ck) || 0) + 1)); });
    const famRole = (c, role) => d3.sum(ctsActive, ct => cnt.get(c).get(colKey(role, ct)) || 0);
    const Ntot = filtered.length;
    const specMemo = new Map(), famRoleMemo = new Map();
    const famRoleC = (c, role) => { const key = c + '|' + role; if (!famRoleMemo.has(key)) famRoleMemo.set(key, famRole(c, role)); return famRoleMemo.get(key); };
    const spec = (c, role, ct) => { const key = c + '|' + role + ct; if (specMemo.has(key)) return specMemo.get(key); const k = cnt.get(c).get(colKey(role, ct)) || 0, n = colN.get(colKey(role, ct)) || 0, K = famRoleC(c, role); const lift = Math.log2(((k + 0.5) / (n + 1)) / ((K + 0.5) / (Ntot + 1))); const out = { k, n, K, lift: Math.max(-3, Math.min(3, lift)), p: k > 0 ? hyperSf(k, Ntot, K, n) : 1 }; specMemo.set(key, out); return out; };
    const MINK = 3, perColUsed = 3;
    const specificRows = role => {
      const best = new Map();  // family -> best p over the role's columns
      ctsActive.forEach(ct => { if (!colN.has(colKey(role, ct))) return;
        fams.map(c => [c, spec(c, role, ct)]).filter(([, s]) => s.k >= MINK && s.lift > 0).sort((a, b) => a[1].p - b[1].p || b[1].lift - a[1].lift).slice(0, perColUsed)
          .forEach(([c, s]) => { if (!best.has(c) || best.get(c) > s.p) best.set(c, s.p); }); });
      return Array.from(best.entries()).sort((a, b) => a[1] - b[1]).map(e => e[0]).slice(0, state.pathTop);
    };
    // average-linkage hierarchical clustering (Euclidean) -> leaf order
    const hclustOrder = vectors => {
      const n = vectors.length; if (n <= 2) return d3.range(n);
      const dist = (a, b) => Math.sqrt(d3.sum(a.map((v, i) => (v - b[i]) ** 2)));
      let clusters = vectors.map((v, i) => ({ items: [i], vec: v.slice(), size: 1 }));
      const D = new Map(); const key = (a, b) => a < b ? a + '|' + b : b + '|' + a; let id = 0; clusters.forEach(c => { c.id = id++; });
      for (let a = 0; a < n; a++) for (let b = a + 1; b < n; b++) D.set(key(a, b), dist(vectors[a], vectors[b]));
      while (clusters.length > 1) {
        let bi = 0, bj = 1, bd = Infinity;
        for (let a = 0; a < clusters.length; a++) for (let b = a + 1; b < clusters.length; b++) { const d = D.get(key(clusters[a].id, clusters[b].id)); if (d < bd) { bd = d; bi = a; bj = b; } }
        const A = clusters[bi], B = clusters[bj];
        const merged = { items: A.items.concat(B.items), size: A.size + B.size, id: id++ };
        clusters.forEach((c, ci) => { if (ci !== bi && ci !== bj) D.set(key(merged.id, c.id), (A.size * D.get(key(A.id, c.id)) + B.size * D.get(key(B.id, c.id))) / merged.size); });
        clusters = clusters.filter((_, ci) => ci !== bi && ci !== bj); clusters.push(merged);
      }
      return clusters[0].items;
    };
    // top overall: rows ranked by the chosen metric, columns in the dataset's cell-type order (a readable ranked list;
    // the clustered, enrichment-based reading is the job of the cell-type-specific view)
    const overallRows = fams.slice().sort((a, b) => score(b) - score(a) || total(b) - total(a)).slice(0, state.pathTop);
    const blocks = [['as sender', 's'], ['as receiver', 'r']].map(([title, role]) => {
      const colsB = ctsActive.filter(ct => colN.has(colKey(role, ct)) || !percell);
      if (!percell) return { title, role, rows: overallRows, cols: ctsActive };
      const rows = specificRows(role); if (rows.length < 1) return { title, role, rows, cols: colsB };
      const M = rows.map(c => colsB.map(ct => spec(c, role, ct).lift));
      const ro = hclustOrder(M), co = hclustOrder(colsB.map((_, j) => M.map(r => r[j])));
      return { title, role, rows: ro.map(i => rows[i]), cols: co.map(j => colsB[j]) };
    });
    const showQ = !!renderPathways.q && !percell;
    const labelW = rows => Math.min(300, 20 + (d3.max(rows, c => Math.min(c.length, 40)) || 10) * 6.2) + (showQ ? 46 : 0);
    const left1 = labelW(blocks[0].rows), left2 = percell ? labelW(blocks[1].rows) : 0;  // in the specific view each block has its own rows and labels
    const top = 24 + Math.min(150, (d3.max(ctsActive, i => CT[i].name.length) || 6) * 6), gap = 40, rh = 16;
    const avail = Math.max(600, (host.node() && host.node().clientWidth) || 1100);
    const nColsMax = Math.max(blocks[0].cols.length, blocks[1].cols.length, 1);
    let cell = Math.max(6, Math.min(22, Math.floor((avail - left1 - left2 - gap - 30) / (2 * nColsMax))));
    const stacked = percell && cell < 10;  // many cell types: receiver block under the sender block so cells stay readable
    if (stacked) cell = Math.max(6, Math.min(22, Math.floor((avail - Math.max(left1, left2) - 30) / nColsMax)));
    const nRows = Math.max(blocks[0].rows.length, blocks[1].rows.length);
    const W = stacked ? Math.max(left1, left2) + cell * nColsMax + 30 : left1 + left2 + 2 * cell * nColsMax + gap + 30;
    const H = (stacked ? top + rh * blocks[0].rows.length + 30 + top + rh * blocks[1].rows.length : top + rh * nRows) + 60;
    const svg = newSvg(host, W, H);
    const seq = d3.scaleSequential(d3.interpolateRgbBasis(['#f6f4ef', '#9fb3d1', '#2f4f7f'])).domain([0, 1]);
    const div = d3.scaleDiverging(t => d3.interpolateRgbBasis(['#3b6fb6', '#f4f2ee', '#d1495b'])(t)).domain([-3, 0, 3]);
    const vmax = d3.max(blocks.flatMap(b => b.rows), c => d3.max(Array.from(cnt.get(c).values()))) || 1;
    const isSel = c => byLib ? (state.setFilter && state.setFilter.lib === libIdx && state.setFilter.name === c) : (state.classes.size === 1 && state.classes.has(c));
    const clickRow = c => {
      if (byLib) { state.setFilter = isSel(c) ? null : { lib: libIdx, set: P.genesets.libraries[libIdx].sets.findIndex(s => s.name === c), name: c }; }
      else { if (isSel(c)) assign(state.classes, classes); else assign(state.classes, [c]); syncClassList(); }
      userChange();
    };
    const drawLabels = (rows, x, oy = 0) => rows.forEach((c, yi) => {
      const t = svg.append('text').attr('x', x).attr('y', oy + top + yi * rh + 11).attr('text-anchor', 'end').attr('class', 'axis-lab').style('cursor', 'pointer').style('font-weight', isSel(c) ? 700 : 400).on('click', () => clickRow(c));
      t.append('tspan').text((c.length > 40 ? c.slice(0, 38) + '…' : c) + ` (${total(c)})`);
      if (showQ) t.append('tspan').attr('fill', renderPathways.q.get(c) < 0.05 ? '#b8323f' : '#999').text(`  q ${fmtP(renderPathways.q.get(c))}`);
    });
    blocks.forEach(({ title, role, rows, cols: colsB }, bi) => {
      const ox = stacked ? Math.max(left1, left2) : (bi === 0 ? left1 : left1 + cell * nColsMax + gap + left2);
      const oy = stacked && bi === 1 ? top + rh * blocks[0].rows.length + 30 : 0;
      if (bi === 0 || percell) drawLabels(rows, ox - 8, oy);
      const g = svg.append('g').attr('transform', `translate(${ox},${oy + top})`);
      g.append('text').attr('x', 0).attr('y', -top + 14).attr('class', 'plot-title').text(title);
      if (percell && !rows.length) { g.append('text').attr('x', 0).attr('y', 12).attr('class', 'axis-lab').attr('fill', '#999').text(`no family is specific to a cell type ${title} (needs ≥ ${MINK} interactions and enrichment above average)`); return; }
      colsB.forEach((ct, xi) => { g.append('text').attr('transform', `translate(${xi * cell + cell / 2 + 3},-6) rotate(-60)`).attr('class', 'axis-lab').style('font-size', Math.min(10, cell - 1) + 'px').text(CT[ct].name); g.append('rect').attr('x', xi * cell).attr('y', -4).attr('width', cell - 1).attr('height', 3).attr('fill', CT[ct].color); });
      rows.forEach((c, yi) => colsB.forEach((ct, xi) => {
        const v = cnt.get(c).get(colKey(role, ct)) || 0; const s = percell ? spec(c, role, ct) : null;
        g.append('rect').attr('x', xi * cell).attr('y', yi * rh).attr('width', cell - 1).attr('height', rh - 1).attr('rx', 2).attr('fill', v ? (percell ? div(s.lift) : seq(v / vmax)) : '#faf9f6')
          .on('mousemove', ev => showTip(ev, `<b>${esc(c)}</b><br>${esc(CT[ct].name)} ${title}: ${v} interaction${v === 1 ? '' : 's'}` + (s ? `<br>${s.lift >= 0 ? '+' : ''}${s.lift.toFixed(2)} log2 vs the average cell type · ${s.k} of its ${s.n} interactions, family ${s.K} of ${Ntot} · p = ${fmtP(s.p)}` : ''))).on('mouseleave', hideTip);
      }));
    });
    // legend
    const lg = svg.append('g').attr('transform', `translate(${stacked ? Math.max(left1, left2) : left1},${H - 34})`); const gid = 'pw' + Math.random().toString(36).slice(2, 8);
    const grad = svg.append('defs').append('linearGradient').attr('id', gid); d3.range(0, 1.01, 0.1).forEach(t => grad.append('stop').attr('offset', t).attr('stop-color', percell ? div(-3 + 6 * t) : seq(t)));
    lg.append('rect').attr('x', 0).attr('y', -8).attr('width', 120).attr('height', 6).attr('fill', `url(#${gid})`);
    lg.append('text').attr('x', -6).attr('y', -2).attr('text-anchor', 'end').attr('class', 'axis-lab').text(percell ? '−3' : '0');
    lg.append('text').attr('x', 126).attr('y', -2).attr('class', 'axis-lab').text(percell ? '+3 log2 enrichment of the family in the cell type vs the average cell type' : `${vmax} interactions`);
    const line1 = byLib ? `${P.genesets.libraries[libIdx].name.replace(/_/g, ' ')} · a pair counts in every set containing its ligand or receptor` : `CellPhoneDB classification · ${total('Unclassified') || 0} interactions without a classification`;
    const line2 = percell ? `rows: the ${perColUsed} most specific families of every cell type (hypergeometric enrichment within its interactions, ≥ ${MINK} interactions), separately as sender and as receiver · rows and columns ordered by hierarchical clustering`
      : `rows: top ${overallRows.length} of ${fams.length} · colour: interactions` + (state.pathRank === 'enrich' ? (wholeUniverse ? ' · enrichment needs a narrower selection (all LR pairs of the result are selected); ranked by interactions instead' : ` · ranked by enrichment of the selection's LR pairs against all ${NN} pairs of the result (hypergeometric, BH)`) : ' · ranked by number of interactions');
    svg.append('text').attr('x', 0).attr('y', H - 18).attr('class', 'axis-lab').attr('fill', '#777').text(line1);
    svg.append('text').attr('x', 0).attr('y', H - 6).attr('class', 'axis-lab').attr('fill', '#777').text(line2);
  }

  // ── reference: query vs reference scatter ──────────────────────────────────
  function renderScatter() {
    const host = d3.select('#scatter'); host.selectAll('*').remove();
    if (!P.meta.has_ref_cs) { host.append('div').attr('class', 'empty').text('Reference scores are not in this result (re-run the comparison with this FastCCC version).'); return; }
    const by = state.scatterBy; const k = (P.calibration && P.calibration.k) || 3, tol = 1.96 / k;
    // axes: combined score (query vs reference), or one side: the query's ligand/receptor summary vs the reference's;
    // the reference side value is the centre of its interval (lo, hi) = ref x (1 -/+ 1.96/k)
    const EPS = 1e-3, FLOOR = -4, lf = (v, n) => (v == null || n == null || n <= 0 ? null : Math.max(FLOOR, Math.log2((v + EPS) / (n + EPS))));  // floored: a score of 0 (not expressed) sits on the floor
    const xOf = by === 'lig' ? i => (R.lig_lo[i] != null && R.lig_hi[i] != null ? (R.lig_lo[i] + R.lig_hi[i]) / 2 : null) : by === 'rec' ? i => (R.rec_lo[i] != null && R.rec_hi[i] != null ? (R.rec_lo[i] + R.rec_hi[i]) / 2 : null) : by === 'fold' ? i => lf(R.ref_cs[i], R.null_cs[i]) : i => R.ref_cs[i];
    const yOf = by === 'lig' ? i => R.lig_cs[i] : by === 'rec' ? i => R.rec_cs[i] : by === 'fold' ? i => lf(R.cs[i], R.null_cs[i]) : i => R.cs[i];
    const pts = shown.filter(i => xOf(i) != null && yOf(i) != null);  // includes Both NS pairs when ticked: the background of non-significant pairs
    if (!pts.length) { host.append('div').attr('class', 'empty').text('Nothing selected.'); return; }
    const W = 520, H = 440, m = { l: 50, r: 16, t: 18, b: 64 };
    const svg = newSvg(host, W, H);
    const fold = by === 'fold';
    const lim = Math.max(d3.max(pts, yOf), d3.max(pts, xOf)) * 1.05 || 1, lo = fold ? Math.min(0, Math.min(d3.min(pts, yOf), d3.min(pts, xOf))) * 1.05 : 0;
    const x = d3.scaleLinear().domain([lo, lim]).range([m.l, W - m.r]), y = d3.scaleLinear().domain([lo, lim]).range([H - m.b, m.t]);
    svg.append('g').attr('transform', `translate(0,${H - m.b})`).call(d3.axisBottom(x).ticks(6)).call(axisStyle);
    svg.append('g').attr('transform', `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(6)).call(axisStyle);
    const what = by === 'lig' ? 'ligand summary in the sender' : by === 'rec' ? 'receptor summary in the receiver' : fold ? 'log2 score / null expectation' : 'communication score';
    svg.append('text').attr('x', (m.l + W - m.r) / 2).attr('y', H - 28).attr('text-anchor', 'middle').attr('class', 'axis-lab').text(`${what} in the reference`);
    svg.append('text').attr('transform', `translate(14,${(m.t + H - m.b) / 2}) rotate(-90)`).attr('text-anchor', 'middle').attr('class', 'axis-lab').text(`${what} in ${P.meta.name}`);
    if (by === 'cs' || fold) svg.append('line').attr('x1', x(lo)).attr('y1', y(lo)).attr('x2', x(lim)).attr('y2', y(lim)).attr('stroke', '#d6d4ce').attr('stroke-dasharray', '3,3');
    if (fold) {  // origin = the null expectation: specific in the query above y = 0, specific in the reference right of x = 0
      svg.append('line').attr('x1', x(lo)).attr('x2', x(lim)).attr('y1', y(0)).attr('y2', y(0)).attr('stroke', '#bdbab2'); svg.append('line').attr('y1', y(lo)).attr('y2', y(lim)).attr('x1', x(0)).attr('x2', x(0)).attr('stroke', '#bdbab2');
      svg.append('text').attr('x', x(lim) - 4).attr('y', y(0) - 4).attr('text-anchor', 'end').attr('class', 'axis-lab').attr('fill', '#999').text('= null expectation');
      svg.append('text').attr('x', x(0) + 4).attr('y', y(lim) + 10).attr('class', 'axis-lab').attr('fill', '#999').text('= null');
      svg.append('text').attr('x', x(lo) + 4).attr('y', y(FLOOR) - 4).attr('class', 'axis-lab').attr('fill', '#999').text(`${FLOOR} = floor (not expressed)`);
    }
    else {
      // the reference interval of the side: query values above the upper line exceed the reference, below the lower line fall short
      const xe = lim / (1 + tol);
      svg.append('path').attr('d', `M${x(0)},${y(0)} L${x(xe)},${y(xe * (1 + tol))} L${x(lim)},${y(lim)} L${x(lim)},${y(lim * (1 - tol))} Z`).attr('fill', '#8c7b6b').attr('fill-opacity', 0.08);
      svg.append('line').attr('x1', x(0)).attr('y1', y(0)).attr('x2', x(xe)).attr('y2', y(xe * (1 + tol))).attr('stroke', '#b8323f').attr('stroke-dasharray', '3,3');
      svg.append('line').attr('x1', x(0)).attr('y1', y(0)).attr('x2', x(lim)).attr('y2', y(lim * (1 - tol))).attr('stroke', '#3b6fb6').attr('stroke-dasharray', '3,3');
      svg.append('text').attr('x', x(xe) - 4).attr('y', y(xe * (1 + tol)) + 12).attr('text-anchor', 'end').attr('class', 'axis-lab').attr('fill', '#b8323f').text(`above the reference interval (×${(1 + tol).toFixed(2)})`);
      svg.append('text').attr('x', x(lim) - 4).attr('y', y(lim * (1 - tol)) + 12).attr('text-anchor', 'end').attr('class', 'axis-lab').attr('fill', '#3b6fb6').text(`below (×${(1 - tol).toFixed(2)})`);
    }
    const order = [4, 3, 0, 1, 2]; const sorted = pts.slice().sort((a, b) => order.indexOf(R.trend[a]) - order.indexOf(R.trend[b]));
    svg.append('g').selectAll('circle').data(sorted).enter().append('circle').attr('cx', i => x(xOf(i))).attr('cy', i => y(yOf(i))).attr('r', pts.length > 5000 ? 1.6 : pts.length > 1000 ? 2.4 : 3.5)
      .attr('fill', i => TREND_COLOR[R.trend[i]]).attr('fill-opacity', i => (R.trend[i] === 1 || R.trend[i] === 2) ? 0.85 : 0.45)
      .attr('stroke', i => state.selected === i ? '#111' : (fold && R.thr_hi[i] != null && R.cs[i] > R.thr_hi[i] ? '#333' : 'none')).attr('stroke-width', i => state.selected === i ? 1.5 : 0.8)
      .style('cursor', 'pointer').on('mousemove', (ev, i) => showTip(ev, rowTip(i))).on('mouseleave', hideTip).on('click', (ev, i) => selectRow(i));
    const legT = [1, 2, 3, 4, 0].filter(t => sorted.some(i => R.trend[i] === t)); let lx = W - m.r - 6 - legT.length * 78; const ly = H - m.b - 14;
    legT.forEach(t => { svg.append('circle').attr('cx', lx).attr('cy', ly).attr('r', 4).attr('fill', TREND_COLOR[t]); svg.append('text').attr('x', lx + 7).attr('y', ly + 4).attr('class', 'axis-lab').text(TREND[t]); lx += 78; });
    if (by === 'cs' || fold) svg.append('text').attr('x', x(lim) - 4).attr('y', y(lim) + 12).attr('text-anchor', 'end').attr('class', 'axis-lab').attr('fill', '#999').text('query = reference');
    svg.append('text').attr('x', 0).attr('y', H - 6).attr('class', 'axis-lab').attr('fill', '#777').text(fold
      ? 'origin = no cell-type specificity (null expectation); diagonal = same as the reference; outlined = query score above its significance band'
      : by === 'cs' ? 'the diagonal is not the decision line: Up/Down are called on the ligand and receptor sides (switch the axes) and the two datasets differ in depth'
      : `Up needs a side above its reference interval (k = ${k.toFixed(2)} calibration), Down a side below; a point inside the band here was called on the other side`);
  }

  // ── reference: what drives the change ──────────────────────────────────────
  function driverOf(i) {
    const up = R.trend[i] === 1, down = R.trend[i] === 2; if (!up && !down) return null;
    const lig = up ? (R.lig_cs[i] != null && R.lig_hi[i] != null && R.lig_cs[i] > R.lig_hi[i]) : (R.lig_cs[i] != null && R.lig_lo[i] != null && R.lig_cs[i] < R.lig_lo[i]);
    const rec = up ? (R.rec_cs[i] != null && R.rec_hi[i] != null && R.rec_cs[i] > R.rec_hi[i]) : (R.rec_cs[i] != null && R.rec_lo[i] != null && R.rec_cs[i] < R.rec_lo[i]);
    return lig && rec ? 'both' : lig ? 'ligand' : rec ? 'receptor' : 'combined';
  }
  function renderDrivers() {
    const host = d3.select('#drivers'); host.selectAll('*').remove();
    const rows = filtered.filter(i => R.trend[i] === 1 || R.trend[i] === 2);
    if (!rows.length) { host.append('div').attr('class', 'empty').text('No Up or Down interactions in the selection.'); return; }
    const cats = ['ligand', 'receptor', 'both', 'combined']; const catColor = { ligand: LIG_INK, receptor: REC_INK, both: '#7b61c0', combined: '#bdbab2' };
    const tally = {};
    rows.forEach(i => { const d = driverOf(i); const t = R.trend[i]; const addTo = ct => { const k = t + '|' + ct; tally[k] = tally[k] || { ligand: 0, receptor: 0, both: 0, combined: 0 }; tally[k][d]++; }; if (d === 'ligand') addTo(R.sender[i]); else if (d === 'receptor') addTo(R.receiver[i]); else { addTo(R.sender[i]); if (R.receiver[i] !== R.sender[i]) addTo(R.receiver[i]); } });
    const active = CT.map((_, i) => i).filter(i => tally['1|' + i] || tally['2|' + i]);
    const rh = 13, left = Math.min(200, 20 + (d3.max(active, i => CT[i].name.length) || 6) * 6.2), W = 520, H = 40 + rh * active.length + 30;
    const svg = newSvg(host, W, H);
    const half = (W - left - 20) / 2, gap = 24;
    const maxv = d3.max(active, i => Math.max(...[1, 2].map(t => d3.sum(cats.map(c => (tally[t + '|' + i] || {})[c] || 0))))) || 1;
    const x = d3.scaleLinear().domain([0, maxv]).range([0, half - gap]);
    [[1, 'Up'], [2, 'Down']].forEach(([t, lab], k) => {
      const g = svg.append('g').attr('transform', `translate(${left + k * half},34)`);
      g.append('text').attr('x', 0).attr('y', -16).attr('class', 'plot-title').attr('fill', TREND_COLOR[t]).text(lab);
      active.forEach((ct, yi) => { let x0 = 0; cats.forEach(c => { const v = (tally[t + '|' + ct] || {})[c] || 0; if (!v) return; g.append('rect').attr('x', x(x0)).attr('y', yi * rh).attr('width', Math.max(0.5, x(x0 + v) - x(x0) - 1)).attr('height', rh - 3).attr('rx', 1.5).attr('fill', catColor[c]).on('mousemove', ev => showTip(ev, `<b>${esc(CT[ct].name)}</b> · ${lab}<br>${v} ${c === 'combined' ? 'combined only' : c + '-side'}`)).on('mouseleave', hideTip); x0 += v; }); });
      g.append('g').attr('transform', `translate(0,${rh * active.length + 2})`).call(d3.axisBottom(x).ticks(4)).call(axisStyle);
    });
    active.forEach((ct, yi) => svg.append('text').attr('x', left - 8).attr('y', 34 + yi * rh + 9).attr('text-anchor', 'end').attr('class', 'axis-lab').text(CT[ct].name));
    let lx = left; cats.forEach(c => { svg.append('rect').attr('x', lx).attr('y', H - 12).attr('width', 10).attr('height', 10).attr('rx', 2).attr('fill', catColor[c]); svg.append('text').attr('x', lx + 13).attr('y', H - 3).attr('class', 'axis-lab').text(c === 'combined' ? 'combined only' : c + ' side'); lx += c === 'combined' ? 110 : 90; });
    host.append('div').attr('class', 'd-note').text('Ligand side: the sender\'s ligand summary left its reference interval (counted for the sender). Receptor side: the receiver\'s receptor summary left its interval (counted for the receiver). Combined only: neither component alone, but their combination crossed the band.');
  }

  // ── split: difference between two conditions ───────────────────────────────
  function condRows(c) { const out = []; for (let i = 0; i < N; i++) if (R.cond[i] === c && passesCore(i, true)) out.push(i); return out; }
  function renderDiff() {
    const host = d3.select('#diffmap'); host.selectAll('*').remove(); const list = d3.select('#difflist'); list.selectAll('*').remove();
    if (KIND !== 'split') return;
    const A = condRows(state.diffA), B = condRows(state.diffB);
    const MA = pairMatrix(A).all, MB = pairMatrix(B).all;
    const rowsD = CT.map((_, i) => i).filter(i => state.senders.has(i)), colsD = CT.map((_, i) => i).filter(i => state.receivers.has(i));
    if (!rowsD.length || !colsD.length) { host.append('div').attr('class', 'empty').text('Tick at least one sender and one receiver.'); return; }
    const na = conditions[state.diffA], nb = conditions[state.diffB];
    matrixHeatmap(host, rowsD, (s, r) => MA[s][r] - MB[s][r], (s, r) => MA[s][r] + MB[s][r], {
      cols: colsD,
      diverging: true, neg: `in ${nb}`, pos: `in ${na}`, legend: `${na} − ${nb}`, selected: state.diffPair,
      tip: (s, r) => `<b>${esc(CT[s].name)}</b> → <b>${esc(CT[r].name)}</b><br>${esc(na)}: ${MA[s][r]} · ${esc(nb)}: ${MB[s][r]}<br><i>click for the LR pairs</i>`,
      click: (s, r, isSel) => { state.diffPair = isSel ? null : [s, r]; renderDiff(); },
    });
    if (!state.diffPair) { list.append('div').attr('class', 'empty').text('Click a cell to list the ligand–receptor pairs present only in one condition.'); return; }
    const [s, r] = state.diffPair;
    const la = new Map(A.filter(i => R.sender[i] === s && R.receiver[i] === r).map(i => [R.lr[i], i])), lb = new Map(B.filter(i => R.sender[i] === s && R.receiver[i] === r).map(i => [R.lr[i], i]));
    const only = (m1, m2) => Array.from(m1.entries()).filter(([k]) => !m2.has(k)).sort((a, b) => (R.cs[b[1]] || 0) - (R.cs[a[1]] || 0));
    const both = Array.from(la.keys()).filter(k => lb.has(k));
    const item = ([k, i]) => `<li><span class="lig">${esc(LR[k].ligand)}</span> – <span class="rec">${esc(LR[k].receptor)}</span> <span class="muted">CS ${fmt2(R.cs[i])}</span></li>`;
    list.html(`<div class="d-title">${esc(CT[s].name)} → ${esc(CT[r].name)}</div>` +
      `<h4>Only in ${esc(na)} (${la.size - both.length})</h4><ul class="plain">${only(la, lb).slice(0, 40).map(item).join('') || '<li class="muted">none</li>'}</ul>` +
      `<h4>Only in ${esc(nb)} (${lb.size - both.length})</h4><ul class="plain">${only(lb, la).slice(0, 40).map(item).join('') || '<li class="muted">none</li>'}</ul>` +
      `<div class="muted">${both.length} shared</div>`);
  }

  // ── communication patterns (NMF) ───────────────────────────────────────────
  // seeded PRNG (mulberry32) so that the factorisation is reproducible for a given selection
  function rng(seed) { let a = seed >>> 0; return () => { a = (a + 0x6D2B79F5) >>> 0; let t = a; t = Math.imul(t ^ (t >>> 15), t | 1); t ^= t + Math.imul(t ^ (t >>> 7), t | 61); return ((t ^ (t >>> 14)) >>> 0) / 4294967296; }; }
  // NMF by multiplicative updates (same scheme as the Python reference implementation); X is n x m (array of Float64Array rows)
  function nmf(X, k, nIter = 300, seed = 0) {
    const n = X.length, m = X[0].length, rnd = rng(seed), eps = 1e-9;
    let W = Array.from({ length: n }, () => Float64Array.from({ length: k }, () => rnd() + 1e-3));
    let Hm = Array.from({ length: k }, () => Float64Array.from({ length: m }, () => rnd() + 1e-3));
    let prevErr = Infinity;
    for (let it = 0; it < nIter; it++) {
      // H <- H * (W^T X) / (W^T W H)
      const WtW = Array.from({ length: k }, (_, a) => Float64Array.from({ length: k }, (_, b) => { let s = 0; for (let i = 0; i < n; i++) s += W[i][a] * W[i][b]; return s; }));
      for (let a = 0; a < k; a++) for (let j = 0; j < m; j++) { let num = 0; for (let i = 0; i < n; i++) num += W[i][a] * X[i][j]; let den = eps; for (let b = 0; b < k; b++) den += WtW[a][b] * Hm[b][j]; Hm[a][j] *= num / den; }
      // W <- W * (X H^T) / (W H H^T)
      const HHt = Array.from({ length: k }, (_, a) => Float64Array.from({ length: k }, (_, b) => { let s = 0; for (let j = 0; j < m; j++) s += Hm[a][j] * Hm[b][j]; return s; }));
      for (let i = 0; i < n; i++) for (let a = 0; a < k; a++) { let num = 0; for (let j = 0; j < m; j++) num += X[i][j] * Hm[a][j]; let den = eps; for (let b = 0; b < k; b++) den += W[i][b] * HHt[b][a]; W[i][a] *= num / den; }
      if (it % 20 === 19) { let err = 0; for (let i = 0; i < n; i++) for (let j = 0; j < m; j++) { let v = 0; for (let a = 0; a < k; a++) v += W[i][a] * Hm[a][j]; err += (X[i][j] - v) ** 2; } if (prevErr - err < 1e-4 * prevErr) break; prevErr = err; }
    }
    return [W, Hm];
  }
  let lastPatterns = null;  // {k, n, note} for the methods text
  function computePatterns(dir) {
    // rows of the current selection; reference comparisons use the interactions significant in the query (Up, Both Sig),
    // or the Down (reference-only) ones when the selection holds nothing else
    let use = filtered, note = '';
    if (KIND === 'reference') { const q = filtered.filter(i => R.trend[i] === 1 || R.trend[i] === 3); if (q.length) use = q; else { use = filtered.filter(i => R.trend[i] === 2); note = 'Down (reference-only) interactions'; } }
    if (use.length < 30) return { error: `Not enough interactions in the selection for pattern discovery (${use.length}, need at least 30).` };
    const role = dir === 'outgoing' ? R.sender : R.receiver;
    const cts = Array.from(new Set(use.map(i => role[i]))).sort((a, b) => a - b), lrs = Array.from(new Set(use.map(i => R.lr[i]))).sort((a, b) => a - b);
    if (cts.length < 3 || lrs.length < 3) return { error: 'Pattern discovery needs at least 3 cell types and 3 ligand–receptor pairs in the selection.' };
    const ci = new Map(cts.map((c, i) => [c, i])), li = new Map(lrs.map((l, j) => [l, j]));
    const X = cts.map(() => new Float64Array(lrs.length)); let mx = 0;
    use.forEach(i => { const w = (KIND === 'reference' && note ? R.ref_cs[i] : R.cs[i]) ?? 1; X[ci.get(role[i])][li.get(R.lr[i])] += w; });
    X.forEach(r => r.forEach(v => { if (v > mx) mx = v; })); if (mx > 0) X.forEach(r => { for (let j = 0; j < r.length; j++) r[j] /= mx; });
    const kAuto = Math.min(6, Math.max(2, Math.round(Math.sqrt(cts.length)))); const kReq = state.patK === 'auto' ? ((P.patterns && P.patterns.k) || kAuto) : +state.patK;
    const k = Math.max(1, Math.min(kReq, cts.length - 1, lrs.length - 1));
    const [W, Hm] = nmf(X, k);
    const Wn = W.map(r => { const s = d3.sum(r) + 1e-12; return Array.from(r, v => v / s); }), Hn = Hm.map(r => { const s = d3.sum(r) + 1e-12; return Array.from(r, v => v / s); });
    lastPatterns = { k, n: use.length, note };
    return { celltypes: cts, lr: lrs, W: Wn, H: Hn, k, n: use.length, note };
  }
  let patternsTimer = null, patternsRenderedAt = -Infinity;
  function schedulePatterns() {
    if (!P.patterns) return;
    clearTimeout(patternsTimer);
    if (performance.now() - patternsRenderedAt > 800) { renderPatterns(); return; }  // first render and isolated changes: immediately
    const host = d3.select('#patterns'); host.selectAll('*').remove(); host.append('div').attr('class', 'empty').text('Computing patterns for the current selection…');
    patternsTimer = setTimeout(renderPatterns, 350);  // rapid changes (slider drags): wait for the last one
  }
  function renderPatterns() {
    const host = d3.select('#patterns'); host.selectAll('*').remove(); patternsRenderedAt = performance.now();
    if (!P.patterns || !d3.sankey) return;
    const pat = computePatterns(state.patDir);
    if (pat.error) { host.append('div').attr('class', 'empty').text(pat.error); return; }
    const k = pat.k;
    const nodes = [], links = [];
    pat.celltypes.forEach(ct => nodes.push({ key: 'c' + ct, label: CT[ct].name, color: CT[ct].color, kind: 's' }));
    for (let p = 0; p < k; p++) nodes.push({ key: 'p' + p, label: `pattern ${p + 1}`, color: PATTERN_COLORS[p % PATTERN_COLORS.length], kind: 'p' });
    const topPerPattern = [];
    for (let p = 0; p < k; p++) { const order = pat.lr.map((lr, j) => [lr, pat.H[p][j]]).sort((a, b) => b[1] - a[1]).slice(0, 6); topPerPattern.push(order); order.forEach(([lr]) => { if (!nodes.some(n => n.key === 'l' + lr)) nodes.push({ key: 'l' + lr, label: `${LR[lr].ligand} – ${LR[lr].receptor}`, color: '#999', kind: 'r' }); }); }
    const id = key => nodes.findIndex(n => n.key === key);
    pat.celltypes.forEach((ct, ci) => { for (let p = 0; p < k; p++) { const w = pat.W[ci][p]; if (w > 0.08) links.push({ source: id('c' + ct), target: id('p' + p), value: w, color: PATTERN_COLORS[p % PATTERN_COLORS.length], tip: `loading ${w.toFixed(2)}` }); } });
    for (let p = 0; p < k; p++) topPerPattern[p].forEach(([lr, h]) => links.push({ source: id('p' + p), target: id('l' + lr), value: Math.max(0.02, h * 3), color: PATTERN_COLORS[p % PATTERN_COLORS.length], tip: `share ${(h * 100).toFixed(1)}% of the pattern` }));
    const W = 620, H = Math.max(280, 15 * Math.max(pat.celltypes.length, nodes.filter(n => n.kind === 'r').length) + 40);
    const svg = sankeyDraw(host, nodes, links, W, H);
    svg.selectAll('text').filter(d => d && d.kind === 'r').style('font-style', 'italic').style('font-size', '9.5px');
    svg.append('text').attr('x', 0).attr('y', 10).attr('class', 'axis-lab').attr('fill', '#777').text(`${state.patDir} · k = ${k} · ${pat.n.toLocaleString()} selected ${pat.note || 'interactions'} · cell types by loading → pattern → defining LR pairs`);
    fitSvgs();
  }

  // ── enrichment of the selection's genes (hypergeometric, in page) ──────────
  function lgamma(z) { const g = 7, c = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7]; if (z < 0.5) return Math.log(Math.PI / Math.sin(Math.PI * z)) - lgamma(1 - z); z -= 1; let x = c[0]; for (let i = 1; i < g + 2; i++) x += c[i] / (z + i); const t = z + g + 0.5; return 0.5 * Math.log(2 * Math.PI) + (z + 0.5) * Math.log(t) - t + Math.log(x); }
  const lchoose = (n, k) => lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1);
  function hyperSf(k, N, K, n) { let p = 0; const hi = Math.min(K, n); for (let x = k; x <= hi; x++) p += Math.exp(lchoose(K, x) + lchoose(N - K, n - x) - lchoose(N, n)); return Math.min(1, p); }
  function renderEnrichment() {
    const host = d3.select('#enrichment'); host.selectAll('*').remove();
    if (!P.genesets) return;
    const lib = P.genesets.libraries[state.enrLib]; const uni = P.genesets.genes; const upos = new Map(uni.map((g, i) => [g, i]));
    const sel = new Set(); filtered.forEach(i => { const lr = LR[R.lr[i]]; lr.ligand.split(',').concat(lr.receptor.split(',')).forEach(g => { const j = upos.get(g.trim()); if (j != null) sel.add(j); }); });
    if (sel.size < 3) { host.append('div').attr('class', 'empty').text('Select at least a few interactions.'); return; }
    const Nn = uni.length, n = sel.size;
    const res = lib.sets.map(s => { const k = s.genes.filter(j => sel.has(j)).length; return { name: s.name, k, K: s.genes.length, p: k ? hyperSf(k, Nn, s.genes.length, n) : 1, genes: s.genes.filter(j => sel.has(j)).map(j => uni[j]) }; }).filter(r => r.k >= 2);
    res.sort((a, b) => a.p - b.p); const m = res.length;
    res.forEach((r, i) => { r.q = Math.min(1, r.p * m / (i + 1)); }); for (let i = m - 2; i >= 0; i--) res[i].q = Math.min(res[i].q, res[i + 1].q);
    const top = res.slice(0, 15);
    window.__topEnrichment = top.length && top[0].q < 0.05 ? { lib: lib.name.replace(/_/g, ' '), name: top[0].name, k: top[0].k, q: top[0].q } : null;
    if (!top.length) { host.append('div').attr('class', 'empty').text('No gene set overlaps the selected genes.'); return; }
    const W = 560, rh = 18, left = 300, H = 30 + rh * top.length + 24;
    const svg = newSvg(host, W, H);
    const x = d3.scaleLinear().domain([0, Math.max(2, d3.max(top, r => -Math.log10(Math.max(r.q, 1e-300))))]).range([0, W - left - 60]);
    const g = svg.append('g').attr('transform', `translate(${left},24)`);
    top.forEach((r, yi) => {
      const v = -Math.log10(Math.max(r.q, 1e-300)); const sig = r.q < 0.05;
      g.append('text').attr('x', -8).attr('y', yi * rh + 12).attr('text-anchor', 'end').attr('class', 'axis-lab').style('font-weight', sig ? 600 : 400).text(r.name.length > 48 ? r.name.slice(0, 46) + '…' : r.name);
      g.append('rect').attr('x', 0).attr('y', yi * rh + 3).attr('width', Math.max(1, x(v))).attr('height', rh - 6).attr('rx', 2).attr('fill', sig ? SINGLE : '#cfcdc6')
        .on('mousemove', ev => showTip(ev, `<b>${esc(r.name)}</b><br>${r.k} of ${r.K} set genes among ${n} selected genes (universe ${Nn})<br>p = ${fmtP(r.p)} · q = ${fmtP(r.q)}<br><i>${esc(r.genes.slice(0, 25).join(', '))}${r.genes.length > 25 ? ' …' : ''}</i>`)).on('mouseleave', hideTip);
      g.append('text').attr('x', x(v) + 4).attr('y', yi * rh + 12).attr('class', 'axis-lab').text(`${r.k}/${r.K}`);
    });
    g.append('line').attr('x1', x(-Math.log10(0.05))).attr('x2', x(-Math.log10(0.05))).attr('y1', -4).attr('y2', rh * top.length).attr('stroke', '#c2385a').attr('stroke-dasharray', '3,3');
    g.append('g').attr('transform', `translate(0,${rh * top.length + 2})`).call(d3.axisBottom(x).ticks(5)).call(axisStyle);
    svg.append('text').attr('x', 0).attr('y', 12).attr('class', 'axis-lab').attr('fill', '#777').text(`−log10 adjusted p · ${n} of ${Nn} database genes in the selection · dashed: q = 0.05`);
  }

  // ── cells: embedding overlay and violins ───────────────────────────────────
  let cellExpr = null;
  function cellMatrix() { if (!cellExpr && P.cells) { const bin = atob(P.cells.expr_b64); cellExpr = new Uint8Array(bin.length); for (let i = 0; i < bin.length; i++) cellExpr[i] = bin.charCodeAt(i); } return cellExpr; }
  function geneValues(gi) { const M = cellMatrix(); const nG = P.cells.genes.length, nC = P.cells.n_cells, mx = P.cells.gene_max[gi]; const out = new Float32Array(nC); for (let c = 0; c < nC; c++) out[c] = M[c * nG + gi] / 255 * mx; return out; }
  function renderCells() {
    if (!P.cells) return;
    const uh = d3.select('#umap'); uh.selectAll('*').remove(); const vh = d3.select('#violin'); vh.selectAll('*').remove();
    const gi = state.cellsGene == null ? 0 : state.cellsGene; const gene = P.cells.genes[gi]; document.getElementById('cells-gene').value = gi;
    const vals = geneValues(gi); const cts = P.cells.celltypes; const colorOf = name => (CT.find(c => c.name === name) || {}).color || '#999';
    const rh = 22, left = 150, W = 560;
    const groups = cts.map((name, ci) => { const v = []; for (let c = 0; c < P.cells.n_cells; c++) if (P.cells.cell_ct[c] === ci) v.push(vals[c]); return { name, v }; }).filter(g => g.v.length >= 5).sort((a, b) => d3.mean(b.v) - d3.mean(a.v));
    const H = 30 + rh * groups.length + 30; const svg = newSvg(vh, W, H);
    const vmax = d3.max(vals) || 1; const x = d3.scaleLinear().domain([0, vmax]).range([left, W - 20]);
    const xs = d3.range(0, vmax + 1e-9, vmax / 40);
    svg.append('text').attr('x', left).attr('y', 14).attr('class', 'plot-title').style('font-style', 'italic').text(`${gene}`).append('tspan').style('font-style', 'normal').style('font-weight', 400).text(`  log-normalised expression per cell type (${P.cells.n_cells.toLocaleString()} cells)`);
    groups.forEach((gr, yi) => {
      const y0 = 30 + yi * rh + rh / 2; const bw = Math.max(vmax / 25, 1.06 * (d3.deviation(gr.v) || vmax / 10) * Math.pow(gr.v.length, -0.2));
      const dens = xs.map(xx => [xx, d3.mean(gr.v, v => Math.exp(-0.5 * ((xx - v) / bw) ** 2)) / bw]);
      const dmax = d3.max(dens, d => d[1]) || 1; const half = rh / 2 - 2;
      const area = d3.area().x(d => x(d[0])).y0(d => y0 - d[1] / dmax * half).y1(d => y0 + d[1] / dmax * half).curve(d3.curveBasis);
      svg.append('path').attr('d', area(dens)).attr('fill', colorOf(gr.name)).attr('fill-opacity', 0.8);
      const med = d3.median(gr.v); svg.append('line').attr('x1', x(med)).attr('x2', x(med)).attr('y1', y0 - half).attr('y2', y0 + half).attr('stroke', '#222').attr('stroke-width', 1.2);
      svg.append('text').attr('x', left - 8).attr('y', y0 + 4).attr('text-anchor', 'end').attr('class', 'axis-lab').text(gr.name);
      svg.append('text').attr('x', W - 18).attr('y', y0 + 4).attr('text-anchor', 'end').attr('class', 'axis-lab').attr('fill', '#888').text(`${Math.round(100 * gr.v.filter(v => v > 0).length / gr.v.length)}%`);
    });
    svg.append('g').attr('transform', `translate(0,${H - 22})`).call(d3.axisBottom(x).ticks(5)).call(axisStyle);
    const emb = P.cells.embedding;
    if (!emb) { uh.append('div').attr('class', 'empty').text('No 2-D embedding (obsm X_umap / X_tsne / X_pca) in the input file.'); return; }
    const Wu = 520, Hu = 440, m = 16; const su = newSvg(uh, Wu, Hu);
    const ex = d3.scaleLinear().domain(d3.extent(emb.x)).range([m, Wu - m]), ey = d3.scaleLinear().domain(d3.extent(emb.y)).range([Hu - m, m + 10]);
    const col = d3.scaleSequential(d3.interpolateRgbBasis(['#e9e7e1', '#9fb3d1', '#2f4f7f', '#1b2a44'])).domain([0, vmax]);
    const order = d3.range(P.cells.n_cells).sort((a, b) => vals[a] - vals[b]);
    su.append('g').selectAll('circle').data(order).enter().append('circle').attr('cx', c => ex(emb.x[c])).attr('cy', c => ey(emb.y[c])).attr('r', P.cells.n_cells > 3000 ? 1.6 : 2.2).attr('fill', c => col(vals[c])).attr('fill-opacity', 0.9);
    cts.forEach((name, ci) => { const xs2 = [], ys2 = []; for (let c = 0; c < P.cells.n_cells; c++) if (P.cells.cell_ct[c] === ci) { xs2.push(emb.x[c]); ys2.push(emb.y[c]); } if (xs2.length < 5) return; su.append('text').attr('x', ex(d3.median(xs2))).attr('y', ey(d3.median(ys2))).attr('text-anchor', 'middle').attr('class', 'axis-lab').attr('paint-order', 'stroke').attr('stroke', 'white').attr('stroke-width', 3).style('font-weight', 600).attr('fill', colorOf(name)).text(name); });
    su.append('text').attr('x', m).attr('y', 12).attr('class', 'plot-title').style('font-style', 'italic').text(gene).append('tspan').style('font-style', 'normal').style('font-weight', 400).text(`  on ${emb.key.replace('X_', '').toUpperCase()}`);
    const lg = su.append('g').attr('transform', `translate(${Wu - 140},${Hu - 14})`); const gid = 'ug' + Math.random().toString(36).slice(2, 7); const grad = su.append('defs').append('linearGradient').attr('id', gid);
    d3.range(0, 1.01, 0.2).forEach(t => grad.append('stop').attr('offset', t).attr('stop-color', col(vmax * t)));
    lg.append('rect').attr('x', 0).attr('y', -6).attr('width', 100).attr('height', 6).attr('fill', `url(#${gid})`); lg.append('text').attr('x', 0).attr('y', 8).attr('class', 'axis-lab').text('0'); lg.append('text').attr('x', 100).attr('y', 8).attr('text-anchor', 'end').attr('class', 'axis-lab').text(vmax.toFixed(1));
  }

  // ── table ──────────────────────────────────────────────────────────────────
  const columns = (() => {
    const base = [
      { key: 'pin', label: '', w: 28 },
      { key: 'sender', label: 'Sender', val: i => CT[R.sender[i]].name },
      { key: 'receiver', label: 'Receiver', val: i => CT[R.receiver[i]].name },
      { key: 'ligand', label: 'Ligand', val: i => LR[R.lr[i]].ligand, cls: 'lig' },
      { key: 'receptor', label: 'Receptor', val: i => LR[R.lr[i]].receptor, cls: 'rec' },
    ];
    if (KIND === 'reference') base.push({ key: 'trend', label: 'Trend', val: i => TREND[R.trend[i]], num: i => R.trend[i] });
    base.push({ key: 'cs', label: 'CS', val: i => fmt2(R.cs[i]), num: i => R.cs[i] ?? -1 });
    if (KIND === 'reference') {
      if (P.meta.has_ref_cs) base.push({ key: 'ref_cs', label: 'CS (ref)', val: i => fmt2(R.ref_cs[i]), num: i => R.ref_cs[i] ?? -1 });
      base.push({ key: 'driver', label: 'Driver', val: i => driverOf(i) || '–' });
      base.push({ key: 'lig_pct', label: 'Lig % (q / ref)', val: i => `${pct(R.lig_pct[i])} / ${pct(R.lig_pct_ref[i])}`, num: i => (R.lig_pct[i] ?? 0) - (R.lig_pct_ref[i] ?? 0) });
      base.push({ key: 'rec_pct', label: 'Rec % (q / ref)', val: i => `${pct(R.rec_pct[i])} / ${pct(R.rec_pct_ref[i])}`, num: i => (R.rec_pct[i] ?? 0) - (R.rec_pct_ref[i] ?? 0) });
    } else {
      base.push({ key: 'p', label: 'p', val: i => fmtP(R.p[i]), num: i => R.p[i] ?? 1 });
      if (modules) base.push({ key: 'mod', label: 'Modules', val: i => modCount(i), num: i => modCount(i), title: 'how many of the scoring modules selected in the sidebar call the interaction significant on their own' });
      if (SIG) base.push({ key: 'signal', label: 'Signal', val: i => SIG[i].join(', ') });
      if (HAS_Z) base.push({ key: 'z', label: 'Effect', val: i => R.z[i] == null ? '–' : R.z[i].toFixed(1), num: i => R.z[i] ?? -1e9, title: 'standardised effect size of the Mean/Average/Arithmetic module: (score − null mean) / null SD' });
      if (HAS_BS) base.push({ key: 'bs', label: 'Support', val: i => pct(R.bs[i]), num: i => R.bs[i] ?? -1, title: 'fraction of bootstrap resamples in which the interaction stays significant' });
      if (HAS_REP) base.push({ key: 'rep', label: 'Replication', val: i => pct(R.rep[i]), num: i => R.rep[i] ?? -1, title: 'fraction of samples in which the interaction is significant on its own' });
      if (P.genes) base.push({ key: 'expr', label: 'Lig % · Rec %', val: i => { const e = exprOf(i); return e ? `${pct(e.lig)} · ${pct(e.rec)}` : '–'; }, num: i => { const e = exprOf(i); return e ? Math.min(e.lig, e.rec) : -1; } });
    }
    base.push({ key: 'spec', label: 'Pairs', val: i => lrSigPairs[R.lr[i]], num: i => lrSigPairs[R.lr[i]], title: 'number of cell-type pairs in which this LR pair is significant (specificity)' });
    base.push({ key: 'class', label: 'Classification', val: i => classOf(LR[R.lr[i]]) });
    return base;
  })();
  function sortedFiltered() { const c = columns.find(c => c.key === state.sort.key); const arr = shown.slice(); if (!c) return arr; if (c.num) arr.sort((a, b) => (c.num(a) - c.num(b)) * state.sort.dir); else if (c.val) arr.sort((a, b) => String(c.val(a)).localeCompare(String(c.val(b))) * state.sort.dir); return arr; }
  function renderTable() {
    const host = document.getElementById('table');
    const arr = sortedFiltered(); const pages = Math.max(1, Math.ceil(arr.length / pageSize));
    state.page = Math.min(state.page, pages - 1);
    const slice = arr.slice(state.page * pageSize, (state.page + 1) * pageSize);
    let h = '<table><thead><tr>' + columns.map(c => `<th data-k="${c.key}" title="${c.title || ''}" class="${state.sort.key === c.key ? 'sorted' : ''}">${c.label}${state.sort.key === c.key ? (state.sort.dir < 0 ? ' ▾' : ' ▴') : ''}</th>`).join('') + '</tr></thead><tbody>';
    slice.forEach(i => {
      h += `<tr data-i="${i}" class="${state.selected === i ? 'sel' : ''}">` + columns.map(c => {
        if (c.key === 'pin') return `<td><input type="checkbox" data-pin="${i}" ${state.pinned.has(i) ? 'checked' : ''} title="add to shortlist"></td>`;
        if (c.key === 'trend') return `<td><span class="sw" style="background:${TREND_COLOR[R.trend[i]]}"></span>${esc(c.val(i))}</td>`;
        return `<td class="${c.cls || ''}">${esc(c.val(i))}</td>`;
      }).join('') + '</tr>';
    });
    h += '</tbody></table>';
    h += `<div class="pager"><span>${arr.length.toLocaleString()} rows · page ${state.page + 1}/${pages}</span> <button data-p="-1">‹</button> <button data-p="1">›</button> <button id="exp-csv">export CSV (filtered)</button></div>`;
    host.innerHTML = h;
    host.querySelectorAll('th').forEach(th => th.addEventListener('click', () => { const k = th.dataset.k; if (k === 'pin') return; state.sort = { key: k, dir: state.sort.key === k ? -state.sort.dir : -1 }; renderTable(); }));
    host.querySelectorAll('tbody tr').forEach(tr => tr.addEventListener('click', ev => { if (ev.target.tagName === 'INPUT') return; selectRow(+tr.dataset.i); }));
    host.querySelectorAll('input[data-pin]').forEach(cb => cb.addEventListener('change', () => { const i = +cb.dataset.pin; cb.checked ? state.pinned.add(i) : state.pinned.delete(i); renderShortlist(); renderDotplot(); }));
    host.querySelectorAll('.pager button[data-p]').forEach(b => b.addEventListener('click', () => { state.page = Math.max(0, Math.min(pages - 1, state.page + (+b.dataset.p))); renderTable(); }));
    document.getElementById('exp-csv').addEventListener('click', () => downloadCSV(arr, 'fastccc_filtered.csv'));
  }
  function rowRecord(i) {
    const lr = LR[R.lr[i]]; const rec = { sender: CT[R.sender[i]].name, receiver: CT[R.receiver[i]].name, ligand: lr.ligand, receptor: lr.receptor, LRI_ID: lr.id };
    if (KIND === 'reference') Object.assign(rec, { trend_vs_ref: TREND[R.trend[i]], comm_score: R.cs[i], ref_comm_score: R.ref_cs[i], driver: driverOf(i), ligand_pct: R.lig_pct[i], ligand_pct_ref: R.lig_pct_ref[i], receptor_pct: R.rec_pct[i], receptor_pct_ref: R.rec_pct_ref[i] });
    else Object.assign(rec, { p_value: R.p[i], comm_score: R.cs[i], effect_size: HAS_Z ? R.z[i] : null, bootstrap_support: HAS_BS ? R.bs[i] : null, replication: HAS_REP ? R.rep[i] : null, modules_significant: modules ? modCount(i) : null, modules_selected: modules ? nSel() : null, signal_type: SIG ? SIG[i].join('; ') : null, condition: conditions ? conditions[R.cond[i]] : null });
    rec.n_celltype_pairs = lrSigPairs[R.lr[i]]; rec.classification = lr.classification; rec.directionality = lr.directionality;
    return rec;
  }
  function downloadCSV(indices, name) {
    const recs = indices.map(rowRecord); if (!recs.length) return; const keys = Object.keys(recs[0]);
    const csv = [keys.join(',')].concat(recs.map(r => keys.map(k => { const v = r[k]; return v == null ? '' : /[",\n]/.test(String(v)) ? '"' + String(v).replace(/"/g, '""') + '"' : v; }).join(','))).join('\n');
    downloadBlob(new Blob([csv], { type: 'text/csv' }), name);
  }

  // ── detail / evidence panel ────────────────────────────────────────────────
  function selectRow(i) {
    state.selected = i; renderDetail(); renderTable(); if (SIG) renderModuleGrid();
    if (P.cells) { const lr = LR[R.lr[i]]; const g = lr.ligand.split(',')[0].trim(); const gi = P.cells.genes.indexOf(g); if (gi >= 0) { state.cellsGene = gi; renderCells(); } }
    if (KIND === 'reference') renderScatter();
  }
  function barChart(sel, title, highlightCt, valueOf, color) {
    const data = CT.map(c => ({ name: c.name, v: valueOf(c.name) })).filter(d => d.v != null).sort((a, b) => b.v - a.v);
    if (!data.length) { sel.append('div').attr('class', 'empty').text(`${title}: expression data not available`); return; }
    const W = 330, rh = 13, left = 125, H = 26 + rh * data.length;
    const svg = newSvg(sel, W, H);
    svg.append('text').attr('x', 0).attr('y', 12).attr('class', 'plot-title').text(title);
    const x = d3.scaleLinear().domain([0, d3.max(data, d => d.v) || 1]).range([0, W - left - 40]);
    const g = svg.append('g').attr('transform', `translate(${left},22)`);
    data.forEach((d, yi) => { const hl = d.name === highlightCt; g.append('text').attr('x', -6).attr('y', yi * rh + 9).attr('text-anchor', 'end').attr('class', 'axis-lab').style('font-weight', hl ? 700 : 400).text(d.name); g.append('rect').attr('x', 0).attr('y', yi * rh + 1).attr('width', Math.max(0.5, x(d.v))).attr('height', rh - 3).attr('rx', 2).attr('fill', color).attr('fill-opacity', hl ? 1 : 0.4); g.append('text').attr('x', x(d.v) + 4).attr('y', yi * rh + 9).attr('class', 'axis-lab').text(pct(d.v)); });
  }
  function curveFor(ctIdx, mid) { if (!P.nulls) return null; return P.nulls.curves[`${ctIdx}|${mid}`] || null; }
  function convolveMean(a, b) {
    const dx = Math.min(a.dx, b.dx); const grid = c => { const out = new Map(); c.y.forEach((p, i) => { const k = Math.round((c.x0 + i * c.dx) / dx); out.set(k, (out.get(k) || 0) + p); }); return out; };
    const A = grid(a), B = grid(b); const out = new Map();
    A.forEach((pa, ka) => B.forEach((pb, kb) => { const k = ka + kb; out.set(k, (out.get(k) || 0) + pa * pb); }));
    const ks = Array.from(out.keys()).sort((p, q) => p - q); const k0 = ks[0], k1 = ks[ks.length - 1]; const y = []; for (let k = k0; k <= k1; k++) y.push(out.get(k) || 0);
    return { x0: k0 * dx / 2, dx: dx / 2, y };
  }
  function renderDetail() {
    const host = d3.select('#detail'); host.selectAll('*').remove();
    const i = state.selected;
    if (i == null) { host.append('div').attr('class', 'empty').text('Click a row or a dot to see its evidence.'); return; }
    const lr = LR[R.lr[i]]; const s = CT[R.sender[i]].name, r = CT[R.receiver[i]].name;
    host.append('div').attr('class', 'd-title').html(`<span class="lig">${esc(lr.ligand)}</span> – <span class="rec">${esc(lr.receptor)}</span><br><span class="d-sub">${esc(s)} → ${esc(r)} · ${esc(classOf(lr))}${lr.directionality ? ' · ' + esc(lr.directionality) : ''}</span>`);
    const facts = host.append('div').attr('class', 'd-facts'); const f = [];
    if (KIND === 'reference') { f.push(['Trend', TREND[R.trend[i]]]); f.push(['CS query', fmt2(R.cs[i])]); if (R.ref_cs[i] != null) f.push(['CS reference', fmt2(R.ref_cs[i])]); f.push(['Null CS', fmt2(R.null_cs[i])]); f.push(['Significance band', `${fmt2(R.thr_lo[i])} – ${fmt2(R.thr_hi[i])}`]); const d = driverOf(i); if (d) f.push(['Driven by', d === 'combined' ? 'the combination only' : d + ' side']); }
    else { f.push(['p (Cauchy)', fmtP(R.p[i])]); f.push(['CS', fmt2(R.cs[i])]); if (modules) f.push(['Modules significant', `${modules.n_sig[i]} / ${nModules}`]); }
    f.push(['Cell-type pairs with this LR', lrSigPairs[R.lr[i]]]);
    facts.html(f.map(([k, v]) => `<div><span class="k">${k}</span><span class="v">${v}</span></div>`).join(''));
    if (KIND === 'reference') {
      const mids = P.nulls && P.nulls.lr_multidata[lr.id]; const cL = mids ? curveFor(R.sender[i], mids[0]) : null, cR = mids ? curveFor(R.receiver[i], mids[1]) : null;
      const W = 340, left = 70, H = (cL && cR) ? 236 : 124; const svg = newSvg(host, W, H);
      svg.append('text').attr('x', 0).attr('y', 12).attr('class', 'plot-title').text(cL && cR ? 'Reference nulls and the query values' : 'Query statistic vs reference-derived band');
      const vals = [R.cs[i], R.ref_cs[i], R.thr_hi[i], R.lig_hi[i], R.rec_hi[i], R.lig_cs[i], R.rec_cs[i]].filter(v => v != null);
      if (cL) vals.push(cL.x0 + cL.dx * cL.y.length); if (cR) vals.push(cR.x0 + cR.dx * cR.y.length);
      const x = d3.scaleLinear().domain([0, (d3.max(vals) || 1) * 1.08]).range([left, W - 12]);
      if (cL && cR) {
        const cC = convolveMean(cL, cR);
        const rows = [['ligand null', cL, '#f4b9b4', R.lig_cs[i], '#c2385a', R.lig_lo[i], R.lig_hi[i], null], ['receptor null', cR, '#a9ddb8', R.rec_cs[i], '#2e8b57', R.rec_lo[i], R.rec_hi[i], null], ['CS null', cC, '#9cb8dc', R.cs[i], '#2f4f7f', R.thr_lo[i], R.thr_hi[i], R.ref_cs[i]]];
        rows.forEach(([lab, c, fill, v, dot, lo, hi, vref], k) => {
          const y0 = 64 + k * 56, hgt = 38; const ymax = d3.max(c.y) || 1;
          const pts = c.y.map((p, j) => [c.x0 + j * c.dx, p / ymax]);
          const area = d3.area().x(d => x(d[0])).y0(y0).y1(d => y0 - d[1] * hgt).curve(d3.curveBasis);
          svg.append('path').attr('d', area(pts)).attr('fill', fill).attr('fill-opacity', 0.85);
          svg.append('line').attr('x1', x(0)).attr('x2', x.range()[1]).attr('y1', y0).attr('y2', y0).attr('stroke', '#e6e4de');
          if (lo != null && hi != null) svg.append('rect').attr('x', x(lo)).attr('y', y0 + 2).attr('width', Math.max(1, x(hi) - x(lo))).attr('height', 5).attr('rx', 2).attr('fill', fill);
          if (v != null) { svg.append('line').attr('x1', x(v)).attr('x2', x(v)).attr('y1', y0).attr('y2', y0 - hgt - 2).attr('stroke', dot).attr('stroke-width', 1.5); svg.append('circle').attr('cx', x(v)).attr('cy', y0 - hgt - 4).attr('r', 4).attr('fill', dot); }
          if (vref != null) { svg.append('line').attr('x1', x(vref)).attr('x2', x(vref)).attr('y1', y0).attr('y2', y0 - hgt - 2).attr('stroke', '#666').attr('stroke-dasharray', '2,2'); svg.append('text').attr('x', x(vref)).attr('y', y0 - hgt - 8).attr('text-anchor', 'middle').attr('class', 'axis-lab').text('ref'); }
          svg.append('text').attr('x', left - 6).attr('y', y0 - 2).attr('text-anchor', 'end').attr('class', 'axis-lab').text(lab);
        });
        svg.append('g').attr('transform', `translate(0,${64 + 2 * 56 + 6})`).call(d3.axisBottom(x).ticks(5)).call(axisStyle);
        host.append('div').attr('class', 'd-note').text('Curves: the reference-derived null of the ligand summary in the sender, of the receptor summary in the receiver, and of their combined score (CS). The bar under each curve is the reference interval; the lollipop is the query; the dashed line is the reference\'s own score.');
      } else {
        const rowsE = [['CS', R.cs[i], R.ref_cs[i], R.thr_lo[i], R.thr_hi[i], '#9cb8dc', '#2f4f7f'], ['ligand', R.lig_cs[i], null, R.lig_lo[i], R.lig_hi[i], '#f4b9b4', '#c2385a'], ['receptor', R.rec_cs[i], null, R.rec_lo[i], R.rec_hi[i], '#a9ddb8', '#2e8b57']];
        rowsE.forEach(([lab, v, vref, lo, hi, band, dot], k) => {
          const y = 36 + k * 26;
          svg.append('text').attr('x', left - 8).attr('y', y + 4).attr('text-anchor', 'end').attr('class', 'axis-lab').text(lab);
          svg.append('line').attr('x1', x(0)).attr('x2', x.range()[1]).attr('y1', y).attr('y2', y).attr('stroke', '#e6e4de');
          if (lo != null && hi != null) svg.append('rect').attr('x', x(lo)).attr('y', y - 7).attr('width', Math.max(1, x(hi) - x(lo))).attr('height', 14).attr('rx', 3).attr('fill', band).attr('fill-opacity', 0.8);
          if (vref != null) { svg.append('line').attr('x1', x(vref)).attr('x2', x(vref)).attr('y1', y - 10).attr('y2', y + 10).attr('stroke', '#666').attr('stroke-width', 2); svg.append('text').attr('x', x(vref)).attr('y', y - 12).attr('text-anchor', 'middle').attr('class', 'axis-lab').text('ref'); }
          if (v != null) svg.append('circle').attr('cx', x(v)).attr('cy', y).attr('r', 5).attr('fill', dot).attr('stroke', 'white').attr('stroke-width', 1.5);
        });
        svg.append('g').attr('transform', `translate(0,${H - 18})`).call(d3.axisBottom(x).ticks(5)).call(axisStyle);
        host.append('div').attr('class', 'd-note').text('Band = significance range derived from the reference null; the ligand and receptor rows show each component against its reference interval.');
      }
      host.append('div').attr('class', 'd-facts').html(`<div><span class="k">Ligand expressed</span><span class="v">${pct(R.lig_pct[i])} of ${esc(s)} · reference ${pct(R.lig_pct_ref[i])}</span></div><div><span class="k">Receptor expressed</span><span class="v">${pct(R.rec_pct[i])} of ${esc(r)} · reference ${pct(R.rec_pct_ref[i])}</span></div>`);
    } else if (modules) {
      const W = 340, H = 150, left = 60; const svg = newSvg(host, W, H);
      svg.append('text').attr('x', 0).attr('y', 12).attr('class', 'plot-title').text(`Agreement of the ${nModules} scoring modules (−log10 p)`);
      const ps = modules.p[i]; const x = d3.scaleBand().domain(d3.range(nModules)).range([left, W - 10]).padding(0.25);
      const y = d3.scaleLinear().domain([0, Math.max(4, d3.max(ps.map(nlogp)))]).range([H - 30, 24]);
      svg.append('line').attr('x1', left).attr('x2', W - 10).attr('y1', y(nlogp(0.05))).attr('y2', y(nlogp(0.05))).attr('stroke', '#c2385a').attr('stroke-dasharray', '3,3');
      svg.append('text').attr('x', left - 4).attr('y', y(nlogp(0.05)) + 4).attr('text-anchor', 'end').attr('class', 'axis-lab').text('0.05');
      const pWhy = k => { const m = MOD[k]; const q = sumQ(m.summary); return q != null && q >= 0.5 ? `p = 1: the ${sumLong(m.summary)} summary is 0 because fewer than ${Math.round((1 - q) * 100)}% of the cells express the ligand or the receptor` : 'p = 1: not significant in this module'; };
      ps.forEach((p, k) => { const one = p != null && p >= 0.999; svg.append('rect').attr('x', x(k)).attr('y', one ? y(0) - 3 : y(nlogp(p))).attr('width', x.bandwidth()).attr('height', one ? 3 : Math.max(0, y(0) - y(nlogp(p)))).attr('rx', 2).attr('fill', p != null && p < 0.05 ? SINGLE : '#d6d4ce').attr('fill-opacity', one ? 0.6 : 1).attr('stroke', one ? '#b9b6ae' : 'none').attr('stroke-dasharray', one ? '2,2' : null).on('mousemove', ev => showTip(ev, `${esc(modules.names[k])}<br>${one ? pWhy(k) : 'p = ' + fmtP(p)}`)).on('mouseleave', hideTip); });
      svg.append('text').attr('x', left).attr('y', H - 10).attr('class', 'axis-lab').text('hover a bar for the module name');
    }
    if (P.genes) { const wrap = host.append('div').attr('class', 'd-expr'); barChart(wrap.append('div'), `${lr.ligand}: fraction of cells expressing`, s, ct => complexPct(lr.ligand, ct), LIG_INK); barChart(wrap.append('div'), `${lr.receptor}: fraction of cells expressing`, r, ct => complexPct(lr.receptor, ct), REC_INK); }
  }

  // ── shortlist ──────────────────────────────────────────────────────────────
  function renderShortlist() {
    const host = document.getElementById('shortlist'); const arr = Array.from(state.pinned);
    document.getElementById('pin-count').textContent = arr.length;
    if (!arr.length) { host.innerHTML = '<div class="empty">Tick rows in the table to build a shortlist, then export it as a table or as the dot plot.</div>'; state.onlyPinned = false; return; }
    host.innerHTML = '<ul>' + arr.map(i => `<li><span class="sw" style="background:${CT[R.sender[i]].color}"></span>${esc(CT[R.sender[i]].name)} → ${esc(CT[R.receiver[i]].name)}: <span class="lig">${esc(LR[R.lr[i]].ligand)}</span> – <span class="rec">${esc(LR[R.lr[i]].receptor)}</span>` + (KIND === 'reference' ? ` <span class="muted">${TREND[R.trend[i]]}</span>` : ` <span class="muted">p ${fmtP(R.p[i])}</span>`) + ` <a href="#" data-un="${i}" title="remove">×</a></li>`).join('') + '</ul>' +
      `<div class="pager"><button id="exp-short">export CSV (shortlist)</button> <button id="exp-svg">download dot plot SVG</button> <button id="pin-filter">${state.onlyPinned ? 'show all again' : 'show only the shortlist'}</button></div>`;
    host.querySelectorAll('a[data-un]').forEach(a => a.addEventListener('click', ev => { ev.preventDefault(); state.pinned.delete(+a.dataset.un); renderShortlist(); if (state.onlyPinned) update(); else { renderTable(); renderDotplot(); } }));
    document.getElementById('exp-short').addEventListener('click', () => downloadCSV(arr, 'fastccc_shortlist.csv'));
    document.getElementById('exp-svg').addEventListener('click', () => downloadSvg('dotplot'));
    document.getElementById('pin-filter').addEventListener('click', () => { state.onlyPinned = !state.onlyPinned; update(); renderShortlist(); });
  }

  // ── methods section ────────────────────────────────────────────────────────
  function renderMethods() {
    const m = P.meta; const el = document.getElementById('methods'); const rows = [];
    rows.push(['Result folder', m.result_dir]); rows.push(['Database', dbLabel(m)]); if (m.database_info && m.database_info.citation) rows.push(['Database citation', m.database_info.citation]); rows.push(['FastCCC', m.fastccc_version]); rows.push(['Generated', m.created]);
    if (m.n_cells) rows.push(['Cells × genes', `${m.n_cells.toLocaleString()} × ${(m.n_genes || 0).toLocaleString()}`]);
    if (m.celltype_key) rows.push(['Cell-type column', m.celltype_key]);
    if (m.matrix) rows.push(['Input matrix', `${m.matrix.kind} (${m.matrix.confidence})` + (m.preprocessing && m.preprocessing.length ? ` → ${m.preprocessing.join(', ')}` : '')]);
    if (m.parameters) rows.push(['Parameters', Object.entries(m.parameters).filter(([k, v]) => v != null && !(Array.isArray(v) && !v.length)).map(([k, v]) => `${k} = ${Array.isArray(v) ? v.join('/') : v}`).join('; ')]);
    if (m.timing) rows.push(['Timing', Object.entries(m.timing).map(([k, v]) => `${k} ${v}s`).join(', ')]);
    if (conditions) rows.push(['Conditions', m.conditions.map(c => `${c.name}: ${c.n_cells != null ? c.n_cells.toLocaleString() : '?'} cells, ${c.n_significant} significant`).join('; ')]);
    if (P.genes) rows.push(['Expression statistics', `${P.genes.genes.length} ligand/receptor genes × ${P.genes.celltypes.length} cell types from ${P.genes.n_cells_used.toLocaleString()} cells`]);
    if (P.cells) rows.push(['Per-cell data', `${P.cells.n_cells.toLocaleString()} cells (stratified subsample), ${P.cells.genes.length} genes` + (P.cells.embedding ? `, embedding ${P.cells.embedding.key}` : '')]);
    if (P.genesets) rows.push(['Gene sets', P.genesets.libraries.map(l => `${l.name} (${l.sets.length} sets with ≥ 3 database genes)`).join('; ')]);
    if (P.patterns) rows.push(['Communication patterns', lastPatterns ? `NMF on the current selection, k = ${lastPatterns.k}${state.patK === 'auto' ? ' (automatic)' : ''}` : 'NMF on the current selection']);
    let h = '<table class="kv">' + rows.map(([k, v]) => `<tr><th>${esc(k)}</th><td>${esc(v ?? '–')}</td></tr>`).join('') + '</table>';
    if (KIND === 'reference') {
      const rf = m.reference; const cov = m.celltype_coverage;
      h += `<h4>Reference panel</h4><table class="kv"><tr><th>Name</th><td>${esc(rf.display_name || rf.name)} (${esc(rf.source || '')}${rf.tissue ? ', ' + esc(rf.tissue) : ''})</td></tr><tr><th>Size</th><td>${(rf.n_celltypes || 0)} cell types, ${(rf.n_cells || 0).toLocaleString()} cells</td></tr>` + (rf.built_from ? `<tr><th>Built from</th><td>${esc(rf.built_from.obs_column)} == ${esc(rf.built_from.control_value)}</td></tr>` : '') +
        (cov ? `<tr><th>Coverage</th><td>${cov.n_covered}/${cov.n_query_types} query cell types (${Math.round(cov.cells_covered_fraction * 100)}% of cells)${cov.uncovered_query_types.length ? '; not in reference: ' + esc(cov.uncovered_query_types.join(', ')) : ''}</td></tr>` : '') + `<tr><th>Mapping</th><td>${esc(m.celltype_mapping_source || 'none')}</td></tr></table>`;
      const map = m.celltype_mapping || {};
      if (Object.keys(map).length) { const byQ = {}; Object.entries(map).forEach(([ref, q]) => { (byQ[q] = byQ[q] || []).push(ref); }); h += '<h4>Reference cell types pooled into each query cell type</h4><table class="kv">' + Object.entries(byQ).map(([q, refs]) => `<tr><th>${esc(q)}</th><td>${refs.map(rn => `${esc(rn)}${rf.celltypes && rf.celltypes[rn] ? ` <span class="muted">(${rf.celltypes[rn].toLocaleString()})</span>` : ''}`).join('; ')}</td></tr>`).join('') + '</table>'; }
    }
    h += `<h4>How to read</h4><div class="prose">` + (KIND === 'reference'
      ? `<p><b>Up</b>: significant in ${esc(m.name)} but not in the reference. <b>Down</b>: significant in the reference but not in ${esc(m.name)}. The significance band comes from the panel's null distributions calibrated to the query. <b>Pairs</b> counts the cell-type pairs in which an LR pair is significant in either dataset: small numbers mean the interaction is specific. <b>Driver</b> says whether the ligand summary (sender), the receptor summary (receiver) or only their combination left the reference interval.</p>`
      : `<p>Each row is a significant sender–receiver × ligand–receptor pair (Cauchy-combined p < 0.05). <b>Modules</b> counts how many of the scoring modules selected in the sidebar (all ${nModules || 16} by default) call it significant on their own (robustness). <b>Pairs</b> counts the cell-type pairs in which the LR pair is significant (specificity). Use the sliders to require specificity, robustness and expression support; "Highlights" applies a strict, dataset-adaptive preset.</p>`) +
      `<p>Pathway enrichment tests the ligand and receptor genes of the current selection against all genes of the LR database (hypergeometric test, Benjamini–Hochberg), so it asks which signalling families are over-represented among the selected interactions, not among all expressed genes. Communication patterns come from non-negative matrix factorisation of the cell type × LR pair matrix of significant interactions, separately for senders (outgoing) and receivers (incoming).</p></div>`;
    el.innerHTML = h;
  }

  // ── effect size, stability, samples × interactions ─────────────────────────
  function renderEffect() {
    const host = d3.select('#effect'); host.selectAll('*').remove();
    if (!HAS_Z) { host.append('div').attr('class', 'empty').text('Effect sizes need the input file (gene means and variances) and the per-module strength tables.'); return; }
    const pts = filtered.filter(i => R.z[i] != null && R.p[i] != null);
    if (!pts.length) { host.append('div').attr('class', 'empty').text('Nothing selected.'); return; }
    const W = 520, H = 360, m = { l: 50, r: 16, t: 20, b: 44 }; const svg = newSvg(host, W, H);
    const x = d3.scaleLinear().domain([Math.min(0, d3.min(pts, i => R.z[i])), (d3.max(pts, i => R.z[i]) || 1) * 1.05]).range([m.l, W - m.r]).nice();
    const y = d3.scaleLinear().domain([0, 16.5]).range([H - m.b, m.t]);
    svg.append('g').attr('transform', `translate(0,${H - m.b})`).call(d3.axisBottom(x).ticks(6)).call(axisStyle); svg.append('g').attr('transform', `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(6)).call(axisStyle);
    svg.append('text').attr('x', (m.l + W - m.r) / 2).attr('y', H - 8).attr('text-anchor', 'middle').attr('class', 'axis-lab').text('effect size: (score − null mean) / null SD, mean-based module');
    svg.append('text').attr('transform', `translate(14,${(m.t + H - m.b) / 2}) rotate(-90)`).attr('text-anchor', 'middle').attr('class', 'axis-lab').text('−log10 p (Cauchy), capped at 16');
    const atFloor = pts.filter(i => nlogp(R.p[i]) >= 15.9).length;
    svg.append('rect').attr('x', m.l).attr('y', y(16.5)).attr('width', W - m.l - m.r).attr('height', y(15.9) - y(16.5)).attr('fill', '#f3f1ec');
    svg.append('text').attr('x', W - m.r - 4).attr('y', y(16.2) + 4).attr('text-anchor', 'end').attr('class', 'axis-lab').attr('fill', '#999').text(`${atFloor} interactions at the p-value floor: only the effect size ranks them`);
    svg.append('g').selectAll('circle').data(pts).enter().append('circle').attr('cx', i => x(R.z[i])).attr('cy', i => y(nlogp(R.p[i]))).attr('r', pts.length > 3000 ? 1.8 : 2.8)
      .attr('fill', i => KIND === 'split' ? COND_COLORS[R.cond[i] % COND_COLORS.length] : SINGLE).attr('fill-opacity', 0.55).attr('stroke', i => state.selected === i ? '#111' : 'none').style('cursor', 'pointer')
      .on('mousemove', (ev, i) => showTip(ev, rowTip(i) + `<br>effect size ${R.z[i].toFixed(1)}`)).on('mouseleave', hideTip).on('click', (ev, i) => selectRow(i));
  }
  function renderStability() {
    const host = d3.select('#stability'); host.selectAll('*').remove();
    if (!HAS_BS) { host.append('div').attr('class', 'empty').text('Run with --bootstrap N (Python: bootstrap=N) to estimate how often each interaction survives resampling the cells.'); return; }
    const vals = filtered.map(i => R.bs[i]).filter(v => v != null);
    if (!vals.length) { host.append('div').attr('class', 'empty').text('Nothing selected.'); return; }
    const b = P.bootstrap || {}; const n = b.n || 10;
    const bins = d3.range(0, n + 1).map(k => ({ k, v: vals.filter(v => Math.round(v * n) === k).length }));
    const W = 520, H = 312, m = { l: 50, r: 16, t: 48, b: 44 }; const svg = newSvg(host, W, H);
    const x = d3.scaleBand().domain(bins.map(d => d.k)).range([m.l, W - m.r]).padding(0.15); const y = d3.scaleLinear().domain([0, d3.max(bins, d => d.v) || 1]).range([H - m.b, m.t]).nice();
    svg.append('g').attr('transform', `translate(0,${H - m.b})`).call(d3.axisBottom(x).tickFormat(k => `${k}/${n}`)).call(axisStyle); svg.append('g').attr('transform', `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(5)).call(axisStyle);
    svg.selectAll('rect.bar').data(bins).enter().append('rect').attr('class', 'bar').attr('x', d => x(d.k)).attr('y', d => y(d.v)).attr('width', x.bandwidth()).attr('height', d => H - m.b - y(d.v)).attr('rx', 2).attr('fill', d => d.k / n >= 0.8 ? SINGLE : '#cfcdc6')
      .on('mousemove', (ev, d) => showTip(ev, `${d.v} interaction${d.v === 1 ? '' : 's'} significant in ${d.k} of ${n} resamples`)).on('mouseleave', hideTip);
    const stable = vals.filter(v => v >= 0.8).length;
    svg.append('text').attr('x', m.l).attr('y', 16).attr('class', 'plot-title').text(`${stable} of ${vals.length} selected interactions (${Math.round(100 * stable / vals.length)}%) are stable`);
    svg.append('text').attr('x', m.l).attr('y', 32).attr('class', 'axis-lab').text(`stable = significant in ≥ 80% of ${n} resamples of ${Math.round((b.fraction || 0.8) * 100)}% of the cells`);
    svg.append('text').attr('x', (m.l + W - m.r) / 2).attr('y', H - 8).attr('text-anchor', 'middle').attr('class', 'axis-lab').text('resamples in which the interaction is significant');
    host.append('div').attr('class', 'd-note').text('Interactions at low support depend on a few cells or sit at the significance boundary; prefer those with high support when choosing candidates.');
  }
  function renderSamples() {
    const host = d3.select('#samples'); host.selectAll('*').remove(); const ph = d3.select('#samplepca'); ph.selectAll('*').remove();
    const rp = P.replication; if (!rp) return;
    const keys = new Set(); for (let i = 0; i < N; i++) if (passesCore(i, true)) keys.add(`${CT[R.sender[i]].name}|${CT[R.receiver[i]].name}\t${LR[R.lr[i]].id}`);
    const keep = rp.top.map((t, j) => [t, j]).filter(([t]) => keys.has(`${t[0]}\t${t[1]}`)).map(([, j]) => j).slice(0, 150);
    if (!keep.length) { host.append('div').attr('class', 'empty').text('None of the recurrent interactions is in the current selection.'); return; }
    const G = rp.groups, top = keep.map(j => rp.top[j]), M = rp.matrix.map(row => keep.map(j => row[j])); const nG = G.length, nK = top.length;
    // order interactions by prevalence, samples by first principal component (power iteration on the centred binary matrix)
    const B = M.map(row => row.slice(0, nK).map(v => v > 0 ? 1 : 0));
    const colMean = d3.range(nK).map(j => d3.mean(B, r => r[j]));
    const C = B.map(r => r.map((v, j) => v - colMean[j]));
    let w = d3.range(nK).map(() => Math.random() - 0.5); for (let it = 0; it < 60; it++) { const t = C.map(r => d3.sum(r.map((v, j) => v * w[j]))); const nw = d3.range(nK).map(j => d3.sum(C.map((r, i) => r[j] * t[i]))); const nrm = Math.sqrt(d3.sum(nw.map(v => v * v))) || 1; w = nw.map(v => v / nrm); }
    const pc1 = C.map(r => d3.sum(r.map((v, j) => v * w[j])));
    // second component: deflate
    let w2 = d3.range(nK).map(() => Math.random() - 0.5); for (let it = 0; it < 60; it++) { const proj = d3.sum(w2.map((v, j) => v * w[j])); w2 = w2.map((v, j) => v - proj * w[j]); const t = C.map(r => d3.sum(r.map((v, j) => v * w2[j]))); const nw = d3.range(nK).map(j => d3.sum(C.map((r, i) => r[j] * t[i]))); const nrm = Math.sqrt(d3.sum(nw.map(v => v * v))) || 1; w2 = nw.map(v => v / nrm); }
    const pc2 = C.map(r => d3.sum(r.map((v, j) => v * w2[j])));
    const order = d3.range(nG).sort((a, b) => pc1[a] - pc1[b]);
    const cw = Math.max(3, Math.min(10, Math.floor(620 / nK))), rh = 14, left = Math.min(200, 20 + (d3.max(G, g => g.length) || 6) * 6.2), top0 = 30;
    const W = left + cw * nK + 20, H = top0 + rh * nG + 50; const svg = newSvg(host, W, H);
    const vmax = d3.max(M, r => d3.max(r)) || 1; const col = d3.scaleSequential(d3.interpolateRgbBasis(['#f6f4ef', '#9fb3d1', '#2f4f7f'])).domain([0, vmax]);
    order.forEach((gi, yi) => {
      svg.append('text').attr('x', left - 8).attr('y', top0 + yi * rh + 10).attr('text-anchor', 'end').attr('class', 'axis-lab').text(G[gi]);
      for (let j = 0; j < nK; j++) { const v = M[gi][j]; svg.append('rect').attr('x', left + j * cw).attr('y', top0 + yi * rh).attr('width', Math.max(0.5, cw - 1)).attr('height', rh - 2).attr('fill', v > 0 ? col(v) : '#faf9f6')
        .on('mousemove', ev => showTip(ev, `<b>${esc(G[gi])}</b><br>${esc(top[j][0])}<br>${esc(top[j][1])}<br>${v > 0 ? 'CS ' + v.toFixed(2) : 'not significant'}`)).on('mouseleave', hideTip); }
    });
    svg.append('text').attr('x', left).attr('y', 14).attr('class', 'plot-title').text(`${nG} samples × ${nK} most recurrent interactions of the selection (colour = communication score, blank = not significant in that sample)`);
    svg.append('text').attr('x', left).attr('y', H - 8).attr('class', 'axis-lab').attr('fill', '#777').text('samples ordered by the first principal component of their presence/absence profile · hover a cell for the interaction');
    // PCA scatter of samples
    const Wp = 420, Hp = 300, m = { l: 44, r: 14, t: 20, b: 40 }; const sp = newSvg(ph, Wp, Hp);
    const x = d3.scaleLinear().domain(d3.extent(pc1)).range([m.l, Wp - m.r]).nice(), y = d3.scaleLinear().domain(d3.extent(pc2)).range([Hp - m.b, m.t]).nice();
    sp.append('g').attr('transform', `translate(0,${Hp - m.b})`).call(d3.axisBottom(x).ticks(4)).call(axisStyle); sp.append('g').attr('transform', `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(4)).call(axisStyle);
    sp.selectAll('circle').data(d3.range(nG)).enter().append('circle').attr('cx', g => x(pc1[g])).attr('cy', g => y(pc2[g])).attr('r', 5).attr('fill', SINGLE).attr('fill-opacity', 0.8).attr('stroke', 'white').on('mousemove', (ev, g) => showTip(ev, `<b>${esc(G[g])}</b><br>${B[g].filter(v => v).length} of the top ${nK} interactions significant`)).on('mouseleave', hideTip);
    sp.selectAll('text.l').data(d3.range(nG)).enter().append('text').attr('x', g => x(pc1[g]) + 7).attr('y', g => y(pc2[g]) + 4).attr('class', 'axis-lab').text(g => G[g]);
    sp.append('text').attr('x', (m.l + Wp - m.r) / 2).attr('y', Hp - 8).attr('text-anchor', 'middle').attr('class', 'axis-lab').text('PC1 of the communication profile'); sp.append('text').attr('transform', `translate(12,${(m.t + Hp - m.b) / 2}) rotate(-90)`).attr('text-anchor', 'middle').attr('class', 'axis-lab').text('PC2');
  }

  // ── signal types (modules) ─────────────────────────────────────────────────
  function renderSignalTypes() {
    const host = d3.select('#signaltypes'); host.selectAll('*').remove();
    if (!SIG) return;
    if (!filtered.length) { host.append('div').attr('class', 'empty').text('Nothing selected.'); return; }
    const per = new Map(); CT.forEach((_, i) => per.set(i, { broad: 0, 'sub-population': 0, 'mean-driven': 0, weak: 0 }));
    const tot = { broad: 0, 'sub-population': 0, 'mean-driven': 0, weak: 0 }; let nComplex = 0, nOne = 0;
    filtered.forEach(i => { const t = SIG[i][0]; per.get(R.sender[i])[t]++; tot[t]++; if (SIG[i].includes('complex-limited')) nComplex++; if (SIG[i].includes('one-sided')) nOne++; });
    const active = CT.map((_, i) => i).filter(i => d3.sum(Object.values(per.get(i))) > 0);
    const rh = 14, left = Math.min(200, 20 + (d3.max(active, i => CT[i].name.length) || 6) * 6.2), W = 520, H = 60 + rh * active.length + 34;
    const svg = newSvg(host, W, H);
    const maxv = d3.max(active, i => d3.sum(Object.values(per.get(i)))) || 1; const x = d3.scaleLinear().domain([0, maxv]).range([0, W - left - 30]);
    const parts = [[tot.broad, 'broad'], [tot['sub-population'], 'sub-population'], [tot['mean-driven'], 'mean-driven'], [tot.weak, 'weak'], [nComplex, 'complex-limited'], [nOne, 'one-sided']].filter(([n]) => n > 0).map(([n, t]) => `${n} ${t}`);
    svg.append('text').attr('x', 0).attr('y', 14).attr('class', 'plot-title').text(`as sender · ${parts.join(', ') || 'nothing selected'}`);
    const g = svg.append('g').attr('transform', `translate(${left},28)`);
    active.forEach((ct, yi) => { let x0 = 0; SIGNAL_TYPES.forEach(tn => { const v = per.get(ct)[tn]; if (!v) return; g.append('rect').attr('x', x(x0)).attr('y', yi * rh).attr('width', Math.max(0.5, x(x0 + v) - x(x0) - 1)).attr('height', rh - 3).attr('rx', 1.5).attr('fill', SIGNAL_COLOR[tn]).on('mousemove', ev => showTip(ev, `<b>${esc(CT[ct].name)}</b> as sender<br>${v} ${tn} signal${v === 1 ? '' : 's'}<br><i>${esc(SIGNAL_HELP[tn])}</i>`)).on('mouseleave', hideTip); x0 += v; }); svg.append('text').attr('x', left - 8).attr('y', 28 + yi * rh + 10).attr('text-anchor', 'end').attr('class', 'axis-lab').text(CT[ct].name); });
    g.append('g').attr('transform', `translate(0,${rh * active.length + 2})`).call(d3.axisBottom(x).ticks(5)).call(axisStyle);
    let lx = left; SIGNAL_TYPES.forEach(tn => { svg.append('rect').attr('x', lx).attr('y', H - 12).attr('width', 10).attr('height', 10).attr('rx', 2).attr('fill', SIGNAL_COLOR[tn]); svg.append('text').attr('x', lx + 13).attr('y', H - 3).attr('class', 'axis-lab').text(tn); lx += 96; });
  }
  function renderModuleGrid() {
    const host = d3.select('#modulegrid'); host.selectAll('*').remove();
    if (!SIG) return;
    const i = state.selected;
    if (i == null) { host.append('div').attr('class', 'empty').text(`Select an interaction to see its ${nModules} module p-values.`); return; }
    const lr = LR[R.lr[i]]; const ps = modules.p[i];
    const cols = []; AGGS.forEach(a => COMBS.forEach(c => cols.push([a, c])));
    const cw = 84, rh = 26, left = 100, top = 70, W = Math.max(560, left + cw * cols.length + 20), H = top + rh * SUMMARIES.length + 40;
    const svg = newSvg(host, W, H);
    svg.append('text').attr('x', 0).attr('y', 14).attr('class', 'plot-title').text(`${lr.ligand} – ${lr.receptor} · ${CT[R.sender[i]].name} → ${CT[R.receiver[i]].name}`).style('font-style', 'italic');
    svg.append('text').attr('x', 0).attr('y', 32).attr('class', 'axis-lab').text(`signal: ${SIG[i].join(', ')} · ${modules.n_sig[i]}/${nModules} modules significant`);
    const color = v => v == null ? '#f3f1ec' : v < 0.05 ? d3.interpolateRgbBasis(['#cdd9ea', '#4b6fa5', '#1b2a44'])(Math.min(1, nlogp(v) / 12)) : '#f3f1ec';
    cols.forEach(([a, c], xi) => { svg.append('text').attr('x', left + xi * cw + cw / 2).attr('y', top - 20).attr('text-anchor', 'middle').attr('class', 'axis-lab').text(a); svg.append('text').attr('x', left + xi * cw + cw / 2).attr('y', top - 8).attr('text-anchor', 'middle').attr('class', 'axis-lab').attr('fill', '#777').text(c); });
    SUMMARIES.forEach((sm, yi) => {
      svg.append('text').attr('x', left - 8).attr('y', top + yi * rh + rh / 2 + 4).attr('text-anchor', 'end').attr('class', 'axis-lab').text(sumLabel(sm));
      cols.forEach(([a, c], xi) => { const k = MOD.findIndex(m => m.summary === sm && m.agg === a && m.comb === c); const v = k >= 0 ? ps[k] : null;
        svg.append('rect').attr('x', left + xi * cw).attr('y', top + yi * rh).attr('width', cw - 3).attr('height', rh - 3).attr('rx', 3).attr('fill', color(v)).on('mousemove', ev => showTip(ev, `${sm} · ${a} · ${c}<br>` + (v != null && v >= 0.999 && (sumQ(sm) ?? 0) >= 0.5 ? `p = 1: the ${sumLong(sm)} summary is 0 (fewer than ${Math.round((1 - sumQ(sm)) * 100)}% of the cells express the ligand or the receptor)` : `p = ${fmtP(v)}`))).on('mouseleave', hideTip);
        svg.append('text').attr('x', left + xi * cw + (cw - 3) / 2).attr('y', top + yi * rh + rh / 2 + 4).attr('text-anchor', 'middle').attr('class', 'axis-lab').attr('fill', v != null && v < 0.05 ? 'white' : '#999').text(fmtP(v)); });
    });
    svg.append('text').attr('x', 0).attr('y', H - 8).attr('class', 'axis-lab').attr('fill', '#777').text('rows: expression summary · columns: complex aggregation × LR combination');
    host.append('div').attr('class', 'd-note').text(SIG[i].map(tn => `${tn}: ${SIGNAL_HELP[tn] || (tn === 'complex-limited' ? 'one subunit of the complex limits the interaction (Minimum aggregation fails, Average passes)' : 'ligand and receptor levels are unbalanced (Arithmetic passes, Geometric fails)')}`).join(' · '));
  }

  // ── reference: disturbance per cell type and comparability QC ──────────────
  // per cell type (as sender or receiver; once when both): gained = Up / (Up + Both Sig), the share of its query communications
  // that are new versus the reference; lost = Down / (Down + Both Sig), the share of its reference communications missing in
  // the query. Computed on the rows passing every filter except the trend ticks, so the denominators do not move with them.
  function disturbanceStats() {
    const st = CT.map((c, i) => ({ i, name: c.name, color: c.color, up: 0, down: 0, both: 0 }));
    for (let i = 0; i < N; i++) { if (![1, 2, 3].includes(R.trend[i]) || !passesCore(i, false, true)) continue;
      const cts = R.sender[i] === R.receiver[i] ? [R.sender[i]] : [R.sender[i], R.receiver[i]];
      cts.forEach(ct => { if (R.trend[i] === 1) st[ct].up++; else if (R.trend[i] === 2) st[ct].down++; else st[ct].both++; }); }
    st.forEach(d => { d.nq = d.up + d.both; d.nr = d.down + d.both; d.gained = d.nq ? d.up / d.nq : null; d.lost = d.nr ? d.down / d.nr : null; d.score = Math.max(d.gained || 0, d.lost || 0); });
    return st.filter(d => d.nq >= 5 || d.nr >= 5).sort((a, b) => b.score - a.score);
  }
  function renderDisturbance() {
    const host = d3.select('#disturb'); host.selectAll('*').remove();
    const pts = disturbanceStats();
    if (!pts.length) { host.append('div').attr('class', 'empty').text('No cell type with at least 5 matched interactions in the selection.'); return; }
    const rh = 26, bh = 9, left = Math.min(200, 20 + (d3.max(pts, d => d.name.length) || 6) * 6.2), W = 560, H = 30 + rh * pts.length + 54;
    const svg = newSvg(host, W, H); const x = d3.scaleLinear().domain([0, 1]).range([left, W - 150]);
    pts.forEach((d, yi) => {
      const y0 = 26 + yi * rh;
      svg.append('text').attr('x', left - 8).attr('y', y0 + rh / 2 + 4).attr('text-anchor', 'end').attr('class', 'axis-lab').text(d.name);
      [['gained', TREND_COLOR[1], d.gained, d.up, d.nq, 'of its query communications are new vs the reference (Up)', 0], ['lost', TREND_COLOR[2], d.lost, d.down, d.nr, 'of its reference communications are missing in the query (Down)', bh + 2]].forEach(([kind, col, v, k, n, what, dy]) => {
        if (v == null) return;
        svg.append('rect').attr('x', x(0)).attr('y', y0 + 2 + dy).attr('width', Math.max(0, x(v) - x(0))).attr('height', bh).attr('rx', 2).attr('fill', col).attr('fill-opacity', 0.75)
          .on('mousemove', ev => showTip(ev, `<b>${esc(d.name)}</b> · ${kind}<br>${k} of ${n} (${Math.round(v * 100)}%) ${what}`)).on('mouseleave', hideTip);
        svg.append('text').attr('x', x(v) + 4).attr('y', y0 + 2 + dy + bh - 1).attr('class', 'axis-lab').style('font-size', '9.5px').attr('fill', col).text(`${Math.round(v * 100)}% (${k}/${n})`);
      });
    });
    svg.append('g').attr('transform', `translate(0,${H - 44})`).call(d3.axisBottom(x).ticks(5, '.0%')).call(axisStyle);
    const lg = svg.append('g').attr('transform', `translate(${left},10)`);
    [[TREND_COLOR[1], 'gained: Up / (Up + Both Sig)'], [TREND_COLOR[2], 'lost: Down / (Down + Both Sig)']].forEach(([col, t], k) => { lg.append('rect').attr('x', k * 200).attr('y', -7).attr('width', 10).attr('height', 8).attr('rx', 2).attr('fill', col).attr('fill-opacity', 0.75); lg.append('text').attr('x', k * 200 + 14).attr('y', 0).attr('class', 'axis-lab').text(t); });
    svg.append('text').attr('x', 0).attr('y', H - 6).attr('class', 'axis-lab').attr('fill', '#777').text('counted as sender or receiver · independent of the trend ticks');
  }

  function renderQC() {
    const host = d3.select('#qc'); host.selectAll('*').remove();
    const cal = P.calibration;
    // detection comparison over matched rows
    const ratios = [];
    for (let i = 0; i < N; i++) { if (![1, 2, 3, 4].includes(R.trend[i])) continue; [[R.lig_pct[i], R.lig_pct_ref[i]], [R.rec_pct[i], R.rec_pct_ref[i]]].forEach(([q, r]) => { if (q != null && r != null && q > 0 && r > 0) ratios.push(Math.log2(q / r)); }); }
    const med = ratios.length ? quantile(ratios, 0.5) : null;
    const lines = [];
    if (cal) lines.push(`Calibration factor k = ${cal.k.toFixed(2)}${cal.k_raw != null && cal.k_raw < 3 ? ` (raw ${cal.k_raw.toFixed(2)}, floored at 3)` : ''} from ${cal.n_hk_genes} housekeeping genes: the query's gene ranks scatter around the reference's with roughly ±${(100 / cal.k).toFixed(0)}% relative spread; the significance bands are widened accordingly.`);
    if (med != null) lines.push(`Detection: over matched interactions the query detects ligands and receptors in ${med >= 0 ? '' : ''}${(Math.pow(2, med) * 100).toFixed(0)}% as many cells as the reference (median ratio ${Math.pow(2, med).toFixed(2)}).` + (med < -0.5 ? ' The query is shallower than the reference.' : med > 0.5 ? ' The query is deeper than the reference.' : ' Depth is comparable.'));
    if (!cal && med == null) { host.append('div').attr('class', 'empty').text('No calibration diagnostics in this result (re-run the comparison with this FastCCC version).'); return; }
    host.append('div').attr('class', 'prose').html(lines.map(l => `<p>${esc(l)}</p>`).join(''));
    const W = 520, H = 230, m = { l: 50, r: 16, t: 14, b: 40 };
    const svg = newSvg(host, W, H);
    if (cal && cal.hk) {
      const half = (W - m.l - m.r) / 2 - 10;
      const x = d3.scaleLinear().domain([0, d3.max(cal.hk.ref) * 1.05]).range([m.l, m.l + half]), y = d3.scaleLinear().domain([0, d3.max(cal.hk.query) * 1.05]).range([H - m.b, m.t]);
      svg.append('g').attr('transform', `translate(0,${H - m.b})`).call(d3.axisBottom(x).ticks(4)).call(axisStyle); svg.append('g').attr('transform', `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(4)).call(axisStyle);
      svg.append('line').attr('x1', x(0)).attr('y1', y(0)).attr('x2', x(Math.min(x.domain()[1], y.domain()[1]))).attr('y2', y(Math.min(x.domain()[1], y.domain()[1]))).attr('stroke', '#d6d4ce').attr('stroke-dasharray', '3,3');
      svg.append('g').selectAll('circle').data(d3.range(cal.hk.gene.length)).enter().append('circle').attr('cx', j => x(cal.hk.ref[j])).attr('cy', j => y(cal.hk.query[j])).attr('r', 2).attr('fill', '#4b6fa5').attr('fill-opacity', 0.35).on('mousemove', (ev, j) => showTip(ev, `<i>${esc(cal.hk.gene[j])}</i><br>reference rank ${cal.hk.ref[j]} · query rank ${cal.hk.query[j]}`)).on('mouseleave', hideTip);
      svg.append('text').attr('x', m.l + half / 2).attr('y', H - 6).attr('text-anchor', 'middle').attr('class', 'axis-lab').text('housekeeping gene rank · reference');
      svg.append('text').attr('transform', `translate(14,${(m.t + H - m.b) / 2}) rotate(-90)`).attr('text-anchor', 'middle').attr('class', 'axis-lab').text('query');
      if (ratios.length) {
        const gx = d3.scaleLinear().domain([-4, 4]).range([m.l + half + 40, W - m.r]); const bins = d3.bin().domain(gx.domain()).thresholds(32)(ratios);
        const gy = d3.scaleLinear().domain([0, d3.max(bins, b => b.length) || 1]).range([H - m.b, m.t]);
        svg.append('g').selectAll('rect').data(bins).enter().append('rect').attr('x', b => gx(b.x0) + 0.5).attr('y', b => gy(b.length)).attr('width', b => Math.max(0.5, gx(b.x1) - gx(b.x0) - 1)).attr('height', b => H - m.b - gy(b.length)).attr('fill', b => (b.x0 + b.x1) / 2 < 0 ? TREND_COLOR[2] : TREND_COLOR[1]).attr('fill-opacity', 0.75);
        svg.append('line').attr('x1', gx(0)).attr('x2', gx(0)).attr('y1', m.t).attr('y2', H - m.b).attr('stroke', '#666').attr('stroke-dasharray', '3,3');
        svg.append('g').attr('transform', `translate(0,${H - m.b})`).call(d3.axisBottom(gx).ticks(5)).call(axisStyle);
        svg.append('text').attr('x', (gx.range()[0] + gx.range()[1]) / 2).attr('y', H - 6).attr('text-anchor', 'middle').attr('class', 'axis-lab').text('log2 (query / reference % expressing)');
      }
    }
  }

  // ── key findings (rule based) ──────────────────────────────────────────────
  function renderFindings() {
    const el = document.getElementById('findings'); const f = [];
    if (!filtered.length) { el.innerHTML = '<li class="muted">Nothing selected.</li>'; return; }
    const out = new Array(CT.length).fill(0), inc = new Array(CT.length).fill(0); filtered.forEach(i => { out[R.sender[i]]++; inc[R.receiver[i]]++; });
    const topS = d3.maxIndex(out), topR = d3.maxIndex(inc);
    f.push(`<b>${filtered.length.toLocaleString()}</b> interactions in the selection across <b>${new Set(filtered.map(i => R.sender[i] * 1000 + R.receiver[i])).size}</b> cell-type pairs; strongest sender <b>${esc(CT[topS].name)}</b> (${out[topS]}), strongest receiver <b>${esc(CT[topR].name)}</b> (${inc[topR]}).`);
    const lrs = topLRs(filtered, 3); if (lrs.length) f.push(`Most recurrent ligand–receptor pairs: ${lrs.map(k => `<span class="lig">${esc(LR[k].ligand)}</span>–<span class="rec">${esc(LR[k].receptor)}</span>`).join(', ')}.`);
    const spec = filtered.slice().sort((a, b) => lrSigPairs[R.lr[a]] - lrSigPairs[R.lr[b]] || (R.cs[b] || 0) - (R.cs[a] || 0))[0];
    if (spec != null) f.push(`Most specific strong interaction: <b>${esc(CT[R.sender[spec]].name)} → ${esc(CT[R.receiver[spec]].name)}</b> via <span class="lig">${esc(LR[R.lr[spec]].ligand)}</span>–<span class="rec">${esc(LR[R.lr[spec]].receptor)}</span> (significant in ${lrSigPairs[R.lr[spec]]} cell-type pair${lrSigPairs[R.lr[spec]] === 1 ? '' : 's'}, CS ${fmt2(R.cs[spec])}).`);
    const cls = new Map(); filtered.forEach(i => { const c = classOf(LR[R.lr[i]]); cls.set(c, (cls.get(c) || 0) + 1); });
    const ranked = Array.from(cls.entries()).sort((a, b) => b[1] - a[1]); const topC = ranked.find(([c]) => c !== 'Unclassified') || ranked[0]; const nUn = cls.get('Unclassified') || 0;
    if (topC && topC[0] !== 'Unclassified') f.push(`Largest signalling family (CellPhoneDB classification): <b>${esc(topC[0])}</b> (${topC[1]} interactions, ${Math.round(100 * topC[1] / filtered.length)}%)` + (nUn ? `; ${nUn} interactions (${Math.round(100 * nUn / filtered.length)}%) have no classification.` : '.'));
    else if (nUn && DB_HAS_CLASSES) f.push(`None of the ${filtered.length} selected interactions has a CellPhoneDB classification.`);
    if (KIND === 'reference') {
      const c = { 1: 0, 2: 0, 3: 0 }; filtered.forEach(i => { if (c[R.trend[i]] != null) c[R.trend[i]]++; });
      f.push(`Versus the reference: <b style="color:${TREND_COLOR[1]}">${c[1]} Up</b>, <b style="color:${TREND_COLOR[2]}">${c[2]} Down</b>, ${c[3]} unchanged (Both Sig).`);
      const drv = { ligand: 0, receptor: 0, both: 0, combined: 0 }; let nd = 0; filtered.forEach(i => { const d = driverOf(i); if (d) { drv[d]++; nd++; } });
      if (nd) f.push(`Of the changed interactions, ${Math.round(100 * drv.ligand / nd)}% are driven by the ligand side, ${Math.round(100 * drv.receptor / nd)}% by the receptor side, ${Math.round(100 * drv.both / nd)}% by both.`);
      const dist = disturbanceStats()[0];
      if (dist) f.push(`Most changed cell type: <b>${esc(dist.name)}</b>` + (dist.gained != null ? `, ${Math.round(dist.gained * 100)}% of its query communications are new (${dist.up} of ${dist.nq})` : '') + (dist.lost != null ? `${dist.gained != null ? ' and' : ','} ${Math.round(dist.lost * 100)}% of its reference communications are missing (${dist.down} of ${dist.nr})` : '') + '.');
    }
    if (SIG) { const c = { broad: 0, 'sub-population': 0 }; filtered.forEach(i => { if (c[SIG[i][0]] != null) c[SIG[i][0]]++; }); f.push(`Signal types: ${c.broad} broad, ${c['sub-population']} carried by a sub-population of cells.`); }
    if (HAS_BS) { const v = filtered.map(i => R.bs[i]).filter(x => x != null); if (v.length) f.push(`Stability: ${v.filter(x => x >= 0.8).length} of ${v.length} selected interactions stay significant in ≥ 80% of the bootstrap resamples.`); }
    if (HAS_REP) { const v = filtered.map(i => R.rep[i]).filter(x => x != null); if (v.length) f.push(`Replication: ${v.filter(x => x >= 0.5).length} of ${v.length} interactions of the full dataset are significant in at least half of the ${P.replication.n_groups} samples on their own.`); }
    if (HAS_Z) { const best = filtered.filter(i => R.z[i] != null).sort((a, b) => R.z[b] - R.z[a])[0]; if (best != null) f.push(`Largest effect size: <b>${esc(CT[R.sender[best]].name)} → ${esc(CT[R.receiver[best]].name)}</b> via <span class="lig">${esc(LR[R.lr[best]].ligand)}</span>–<span class="rec">${esc(LR[R.lr[best]].receptor)}</span>, ${R.z[best].toFixed(1)} null SDs above expectation.`); }
    if (P.genesets && window.__topEnrichment) f.push(`Top enriched gene set (${esc(window.__topEnrichment.lib)}): <b>${esc(window.__topEnrichment.name)}</b>, ${window.__topEnrichment.k} genes, q = ${fmtP(window.__topEnrichment.q)}.`);
    el.innerHTML = f.map(x => `<li>${x}</li>`).join('');
  }

  // ── methods text ───────────────────────────────────────────────────────────
  function renderMethodsText() {
    const m = P.meta; const el = document.getElementById('methods-text');
    const v = m.fastccc_version || '1.x'; const db = dbLabel(m) || 'CellPhoneDB v5';
    const cells = m.n_cells ? `${m.n_cells.toLocaleString()} cells` : 'the dataset';
    const prep = m.preprocessing && m.preprocessing.length ? `Raw counts were normalised (${m.preprocessing.join(', ')}). ` : '';
    const par = m.parameters || {};
    let txt = `Cell–cell communication was inferred with FastCCC v${v} (Hou, Ma & Zhou, Nat Commun 16, 11428, 2025) using the ${db} ligand–receptor database on ${cells} annotated into ${CT.length} cell types (column "${m.celltype_key || 'cell_type'}"). ${prep}`;
    if (KIND !== 'reference') txt += `FastCCC computes the null distribution of each communication score analytically (fast Fourier transform convolution) and combines ${nModules || 16} scoring modules (${SUMMARIES.map(sumLong).join(', ')} expression summaries × ${AGGS.map(a => a.toLowerCase()).join(' or ')} complex aggregation × ${COMBS.map(c => c.toLowerCase()).join(' or ')} ligand–receptor combination) with the Cauchy combination test; interactions with a combined p < 0.05 and ligand and receptor expressed in at least ${Math.round((par.min_percentile ?? 0.1) * 100)}% of the sender and receiver cells were called significant. `;
    else {
      const rf = m.reference || {}; const cov = m.celltype_coverage || {};
      txt += `Communication scores were compared against the FastCCC ${rf.source === 'healthy_tissue' ? 'healthy ' + (rf.tissue || '').replace('_', ' ') + ' reference panel' : 'reference panel "' + (rf.display_name || rf.name) + '"'} (${(rf.n_cells || 0).toLocaleString()} cells, ${rf.n_celltypes || 0} cell types${rf.built_from ? ', built from the ' + rf.built_from.obs_column + ' = ' + rf.built_from.control_value + ' cells of the same dataset' : ''}). `;
      txt += `Query cell types were matched to ${cov.n_covered ?? '?'} of ${cov.n_query_types ?? '?'} reference cell types${m.celltype_mapping_source === 'auto' ? ' with FastCCC\'s automatic name matching, reviewed manually' : m.celltype_mapping_source === 'user' ? ' with a curated mapping' : ' by name'}. `;
      txt += `For each interaction, the reference-derived null distribution of the query cell-type pair (rank-based, calibrated on ${P.calibration ? P.calibration.n_hk_genes + ' housekeeping genes, k = ' + P.calibration.k.toFixed(2) : 'housekeeping genes'}) defined a significance band; interactions significant in the query but not in the reference were labelled Up, those significant in the reference but not in the query Down. `;
    }
    if (HAS_BS && P.bootstrap) txt += `Robustness of each interaction was assessed by re-running the analysis on ${P.bootstrap.n} stratified subsamples of ${Math.round((P.bootstrap.fraction || 0.8) * 100)}% of the cells of every cell type and recording the fraction of subsamples in which it remained significant (bootstrap support). `;
    if (HAS_REP && P.replication) txt += `Replication was assessed by running FastCCC separately in each of ${P.replication.n_groups} samples and recording, for every interaction found in the full dataset, the fraction of samples in which it was significant on its own. `;
    if (HAS_Z) txt += `Standardised effect sizes were computed for the mean-based module as (score − null mean) / null standard deviation, where the null moments follow from FastCCC's permutation model (the per-cell-type mean of a gene is the mean of n independent draws from its distribution over all cells). `;
    if (P.genesets) txt += `Pathway enrichment of the ligand and receptor genes of selected interactions was assessed with a hypergeometric test against all genes of the ligand–receptor database (Benjamini–Hochberg correction), using the ${P.genesets.libraries.map(l => l.name.replace(/_/g, ' ')).join(', ')} gene sets from Enrichr. `;
    if (P.patterns && lastPatterns) txt += `Communication patterns were identified by non-negative matrix factorisation (k = ${lastPatterns.k}, computed in the page on the ${lastPatterns.n.toLocaleString()} selected interactions) of the cell type × ligand–receptor matrix of significant interactions, separately for outgoing and incoming signalling. `;
    txt += `Figures were exported from the FastCCC result explorer.`;
    el.innerHTML = `<p>${esc(txt)}</p>`;
    document.getElementById('copy-methods').addEventListener('click', () => { navigator.clipboard && navigator.clipboard.writeText(txt); document.getElementById('copy-methods').textContent = 'copied'; setTimeout(() => { document.getElementById('copy-methods').textContent = 'copy'; }, 1500); });
  }

  // ── figure builder ─────────────────────────────────────────────────────────
  const FIG_PANELS = [['chord', 'Communication flow (chord)'], ['heatmap', 'Sender × receiver'], ['dotplot', 'Dot plot'], ['roles', 'Signalling roles'], ['sankey', 'Sender → LR → receiver'], ['pathways', 'Pathway families'],
    ['scatter', 'Query vs reference'], ['drivers', 'Drivers'], ['disturb', 'Disturbance'], ['qc', 'Comparability'], ['diffmap', 'Condition difference'], ['patterns', 'Communication patterns'], ['enrichment', 'Enrichment'], ['signaltypes', 'Signal types'], ['modulegrid', 'Module grid'], ['umap', 'Embedding'], ['violin', 'Violins'], ['detail', 'Evidence']];
  function buildFigureBuilder() {
    const host = document.getElementById('figure-picks');
    const defaults = new Set(['chord', 'dotplot', 'roles', 'enrichment']);
    host.innerHTML = FIG_PANELS.filter(([id]) => document.getElementById(id) && document.getElementById(id).closest('.card').style.display !== 'none' && document.getElementById(id).closest('.row2, main') && getComputedStyle(document.getElementById(id).closest('.card')).display !== 'none')
      .map(([id, lab]) => `<label><input type="checkbox" data-f="${id}" ${defaults.has(id) ? 'checked' : ''}>${lab}</label>`).join('');
    document.getElementById('fig-build').addEventListener('click', () => {
      const picks = Array.from(host.querySelectorAll('input:checked')).map(cb => cb.dataset.f).filter(id => document.querySelector(`#${id} svg`));
      if (!picks.length) { document.getElementById('fig-note').textContent = 'nothing to compose: tick panels that are currently drawn'; return; }
      const { root, W, H } = composeSvg(picks.map(id => document.querySelector(`#${id} svg`)), true);
      downloadBlob(svgFile(root), 'fastccc_figure.svg');
      document.getElementById('fig-note').textContent = `composed ${picks.length} panel${picks.length === 1 ? '' : 's'} (${W} × ${Math.round(H)} px); open in Inkscape or Illustrator to arrange and edit text`;
    });
  }

  // ── orchestration ──────────────────────────────────────────────────────────
  function update() {
    refilter(); state.page = 0; renderTiles(); renderHeatmap(); renderChord(); renderDotplot(); renderRoles(); renderSankey(); renderPathways(); renderTable();
    if (KIND === 'reference') { renderScatter(); renderDrivers(); renderDisturbance(); }
    if (KIND === 'split') renderDiff();
    if (P.genesets) renderEnrichment();
    if (SIG) renderSignalTypes();
    if (KIND !== 'reference') { renderEffect(); renderStability(); }
    schedulePatterns(); if (P.replication) renderSamples();
    renderFindings(); fitSvgs();
    document.getElementById('pair-focus').textContent = state.pair ? `focused on ${CT[state.pair[0]].name} → ${CT[state.pair[1]].name} (click the cell again to clear)` : '';
  }
  buildFilters(); renderMethods(); renderShortlist(); if (P.cells) renderCells();
  if (!SIG) document.getElementById('signal-row').style.display = 'none';
  if (KIND === 'reference') document.getElementById('stability-row').style.display = 'none';
  if (!P.replication) document.getElementById('samples-card').style.display = 'none'; else renderSamples();
  if (KIND !== 'reference') document.getElementById('refqc-row').style.display = 'none'; else renderQC();
  renderMethodsText(); if (SIG) renderModuleGrid(); buildFigureBuilder();
  const hash = new URLSearchParams(location.hash.replace(/^#/, ''));
  setPreset(hash.has('highlights') ? 'highlights' : 'default');
  if (hash.has('pair')) { const [a, b] = hash.get('pair').split(',').map(Number); if (CT[a] && CT[b]) { state.pair = [a, b]; update(); } }
  if (hash.has('row')) { const i = +hash.get('row'); if (i >= 0 && i < N) { state.pinned.add(i); renderShortlist(); selectRow(i); } }
  renderDetail(); renderTable(); renderDotplot(); fitSvgs();
})();
