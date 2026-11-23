const status = document.querySelector('#status');
const content = document.querySelector('#explorer-content');
const search = document.querySelector('#search');
const list = document.querySelector('#cell-list');
const elements = Object.fromEntries([
  'cell-count', 'search-count', 'selected-id', 'cell-context', 'final-label',
  'pruned-label', 'pruning-detail', 'delta-next', 'delta-median', 'scores', 'markers', 'marker-description'
].map((id) => [id, document.getElementById(id)]));

function parseTsv(text) {
  const lines = text.trimEnd().split(/\r?\n/);
  const header = lines.shift().split('\t');
  return lines.filter(Boolean).map((line) => {
    const values = line.split('\t');
    if (values.length !== header.length) throw new Error('Invalid TSV row');
    return Object.fromEntries(header.map((key, index) => [key, values[index]]));
  });
}

function makeElement(tag, className, text) {
  const node = document.createElement(tag);
  if (className) node.className = className;
  if (text !== undefined) node.textContent = text;
  return node;
}

async function loadTable(name) {
  const response = await fetch(`./data/processed/${name}.tsv`);
  if (!response.ok) throw new Error(`Could not load ${name}.tsv (${response.status}).`);
  return parseTsv(await response.text());
}

function indexRows(rows, key) {
  const index = new Map();
  for (const row of rows) {
    if (!index.has(row[key])) index.set(row[key], []);
    index.get(row[key]).push(row);
  }
  return index;
}

function start(cells, scoreRows, markerRows, expressionRows, coordinateRows, referenceRows) {
  if (!cells.length) throw new Error('The cell table is empty.');
  const scoresByCell = indexRows(scoreRows, 'cell_id');
  const expressionByCell = indexRows(expressionRows, 'cell_id');
  const referenceExpression = new Map(referenceRows.map((row) => [`${row.label}\0${row.gene_id}`, row]));
  const markersByPair = new Map();
  for (const row of markerRows) {
    const key = `${row.label}\0${row.comparison}`;
    if (!markersByPair.has(key)) markersByPair.set(key, []);
    markersByPair.get(key).push(row);
  }
  let selected = cells[0];
  let selectedCandidate = null;
  let selectedBin = null;
  const coordinates = new Map(coordinateRows.map((row) => [row.cell_id, row]));
  if (coordinates.size !== cells.length || cells.some((cell) => !coordinates.has(cell.cell_id))) throw new Error('PCA coordinates do not match the query cells.');
  const palette = ['#316c58', '#bb7045', '#5b75a5', '#9c6690', '#9c853b', '#528a89', '#a7625f', '#7c7961'];
  const labels = [...new Set(cells.map((cell) => cell.singleR_label))].sort();
  const colors = new Map(labels.map((label, index) => [label, palette[index % palette.length]]));
  const labelFilter = document.querySelector('#label-filter');
  const prunedFilter = document.querySelector('#pruned-filter');
  const plot = document.querySelector('#pca-plot');
  const deltaPlot = document.querySelector('#delta-plot');
  const binWidth = 0.05;
  for (const label of labels) {
    const option = makeElement('option', '', label);
    option.value = label;
    labelFilter.append(option);
    const legend = makeElement('span', 'legend-item');
    const swatch = makeElement('span', 'legend-swatch');
    swatch.style.backgroundColor = colors.get(label);
    legend.append(swatch, document.createTextNode(label));
    document.querySelector('#plot-legend').append(legend);
  }
  elements['cell-count'].textContent = String(cells.length);

  function getFilteredCells(includeBin = true) {
    const term = search.value.trim().toLowerCase();
    return cells.filter((cell) => {
      const matchesSearch = `${cell.cell_id} ${cell.singleR_label}`.toLowerCase().includes(term);
      const matchesLabel = labelFilter.value === 'all' || cell.singleR_label === labelFilter.value;
      const matchesPruned = !prunedFilter.checked || cell.is_pruned === 'TRUE';
      const delta = Number(cell.delta_next);
      const matchesBin = !includeBin || selectedBin === null || (cell.delta_next !== '' && Math.floor(delta / binWidth) === selectedBin);
      return matchesSearch && matchesLabel && matchesPruned && matchesBin;
    });
  }

  function svgElement(tag, attributes = {}) {
    const node = document.createElementNS('http://www.w3.org/2000/svg', tag);
    for (const [key, value] of Object.entries(attributes)) node.setAttribute(key, value);
    return node;
  }

  function selectCell(cell, restoreFocus = false) {
    selected = cell;
    selectedCandidate = null;
    renderAll();
    if (restoreFocus) list.querySelector('.cell-option.active')?.focus();
  }

  function renderOverview() {
    const visible = getFilteredCells();
    const points = cells.map((cell) => coordinates.get(cell.cell_id));
    const xs = points.map((point) => Number(point.pc1));
    const ys = points.map((point) => Number(point.pc2));
    if (xs.some((value) => !Number.isFinite(value)) || ys.some((value) => !Number.isFinite(value))) throw new Error('PCA contains invalid coordinates.');
    const xMin = Math.min(...xs), xMax = Math.max(...xs), yMin = Math.min(...ys), yMax = Math.max(...ys);
    const mapX = (value) => 56 + (value - xMin) / (xMax - xMin || 1) * 558;
    const mapY = (value) => 360 - (value - yMin) / (yMax - yMin || 1) * 330;
    plot.replaceChildren();
    plot.append(svgElement('line', { x1: 56, y1: 360, x2: 614, y2: 360, class: 'axis-line' }), svgElement('line', { x1: 56, y1: 30, x2: 56, y2: 360, class: 'axis-line' }));
    const xLabel = svgElement('text', { x: 320, y: 397, 'text-anchor': 'middle', class: 'axis-label' }); xLabel.textContent = 'PC1';
    const yLabel = svgElement('text', { x: 16, y: 200, transform: 'rotate(-90 16 200)', 'text-anchor': 'middle', class: 'axis-label' }); yLabel.textContent = 'PC2';
    plot.append(xLabel, yLabel);
    for (const cell of visible) {
      const point = coordinates.get(cell.cell_id);
      const circle = svgElement('circle', { cx: mapX(Number(point.pc1)), cy: mapY(Number(point.pc2)), r: cell.cell_id === selected.cell_id ? 7 : 4, fill: colors.get(cell.singleR_label), class: cell.cell_id === selected.cell_id ? 'cell-point selected-point' : 'cell-point' });
      const title = svgElement('title'); title.textContent = `${cell.cell_id} · ${cell.singleR_label}${cell.is_pruned === 'TRUE' ? ' · pruned' : ''}`; circle.append(title);
      circle.addEventListener('click', () => selectCell(cell));
      plot.append(circle);
    }
    const selectionVisible = visible.some((cell) => cell.cell_id === selected.cell_id);
    document.querySelector('#plot-summary').textContent = `${visible.length} of ${cells.length} cells visible. Colours indicate final predicted labels; outlined point is selected.${selectionVisible ? '' : ' The selected cell is outside the current filters.'}`;
  }

  function renderDistribution() {
    const filtered = getFilteredCells(false).filter((cell) => cell.delta_next !== '' && Number.isFinite(Number(cell.delta_next)) && Number(cell.delta_next) >= 0);
    const maxDelta = Math.max(0, ...filtered.map((cell) => Number(cell.delta_next)));
    const countBins = Math.max(1, Math.floor(maxDelta / binWidth) + 1);
    const bins = Array.from({ length: countBins }, () => 0);
    for (const cell of filtered) bins[Math.floor(Number(cell.delta_next) / binWidth)]++;
    const maxCount = Math.max(1, ...bins);
    deltaPlot.replaceChildren();
    const width = 430 / countBins;
    bins.forEach((count, index) => {
      const height = count / maxCount * 175;
      const bar = svgElement('rect', { x: 46 + index * width, y: 200 - height, width: Math.max(1, width - 3), height, class: selectedBin === index ? 'histogram-bar active-bin' : 'histogram-bar', tabindex: '0', role: 'button', 'aria-label': `Delta ${ (index * binWidth).toFixed(2) } to ${ ((index + 1) * binWidth).toFixed(2) }: ${count} cells` });
      const title = svgElement('title'); title.textContent = `${count} cells: ${ (index * binWidth).toFixed(2) } ≤ delta.next < ${ ((index + 1) * binWidth).toFixed(2) }`; bar.append(title);
      const choose = (restoreFocus = false) => {
        selectedBin = selectedBin === index ? null : index;
        renderAll();
        if (restoreFocus) deltaPlot.querySelectorAll('.histogram-bar')[index]?.focus();
      };
      bar.addEventListener('click', () => choose());
      bar.addEventListener('keydown', (event) => { if (event.key === 'Enter' || event.key === ' ') { event.preventDefault(); choose(true); } });
      deltaPlot.append(bar);
    });
    deltaPlot.append(svgElement('line', { x1: 46, y1: 200, x2: 476, y2: 200, class: 'axis-line' }));
    const axis = svgElement('text', { x: 260, y: 237, 'text-anchor': 'middle', class: 'axis-label' }); axis.textContent = 'delta.next (bin width 0.05)'; deltaPlot.append(axis);
    const start = svgElement('text', { x: 46, y: 216, class: 'axis-label' }); start.textContent = '0'; deltaPlot.append(start);
    const end = svgElement('text', { x: 476, y: 216, 'text-anchor': 'end', class: 'axis-label' }); end.textContent = (countBins * binWidth).toFixed(2); deltaPlot.append(end);
    document.querySelector('#clear-bin').hidden = selectedBin === null;
    document.querySelector('#bin-summary').textContent = selectedBin === null ? `${filtered.length} cells with available delta.next.` : `Selected: ${ (selectedBin * binWidth).toFixed(2) } ≤ delta.next < ${ ((selectedBin + 1) * binWidth).toFixed(2) }. Histogram excludes the bin filter.`;
  }

  const steps = [
    ['Reference & query', 'Reference and query use different patients from the same study. The PCA shows query cells only; select any point or use the searchable list.', '#reference-overview'],
    ['Initial scores', 'SingleR compares a query cell with reference profiles and reports initial scores for each label. Click a candidate in the score chart to change the marker comparison.', '#score-view'],
    ['Marker evidence', 'SingleR selects genes distinguishing pairs of reference labels. Compare the selected query cell’s raw counts with reference summaries; the scales differ.', '#marker-view'],
    ['Fine-tuning & uncertainty', 'Fine-tuning may change the final label. Delta.next and pruning are distinct diagnostics. Iteration-by-iteration fine-tuning scores are not available in this export.', '#diagnostics']
  ];
  let currentStep = 0;
  function renderStep() {
    document.querySelector('#step-title').textContent = steps[currentStep][0];
    document.querySelector('#step-text').textContent = steps[currentStep][1];
    document.querySelector('#step-link').href = steps[currentStep][2];
    for (const button of document.querySelectorAll('[data-step]')) button.setAttribute('aria-pressed', String(Number(button.dataset.step) === currentStep));
    document.querySelector('#step-previous').disabled = currentStep === 0;
    document.querySelector('#step-next').disabled = currentStep === steps.length - 1;
  }
  for (const button of document.querySelectorAll('[data-step]')) button.addEventListener('click', () => { currentStep = Number(button.dataset.step); renderStep(); });
  document.querySelector('#step-previous').addEventListener('click', () => { currentStep--; renderStep(); });
  document.querySelector('#step-next').addEventListener('click', () => { currentStep++; renderStep(); });

  function renderList() {
    const matches = getFilteredCells();
    elements['search-count'].textContent = `${matches.length} matching ${matches.length === 1 ? 'cell' : 'cells'}`;
    list.replaceChildren();
    for (const cell of matches) {
      const button = makeElement('button', `cell-option${selected.cell_id === cell.cell_id ? ' active' : ''}`);
      button.type = 'button';
      button.setAttribute('aria-pressed', String(selected.cell_id === cell.cell_id));
      button.append(makeElement('span', 'cell-id', cell.cell_id), makeElement('span', 'cell-prediction', cell.singleR_label));
      button.addEventListener('click', () => selectCell(cell, true));
      list.append(button);
    }
    if (!matches.length) list.append(makeElement('p', 'empty', 'No cells match your search.'));
  }

  function renderEvidence() {
    elements['selected-id'].textContent = selected.cell_id;
    elements["cell-context"].textContent = "10x Genomics | Donor 2 query";
    elements['final-label'].textContent = selected.singleR_label;
    const pruned = selected.is_pruned === 'TRUE';
    elements['pruned-label'].textContent = selected.pruned_label || 'Not resolved (NA)';
    elements['pruning-detail'].textContent = pruned ? 'SingleR pruned this assignment' : 'SingleR retained this assignment';
    elements['delta-next'].textContent = selected.delta_next === '' ? 'Not available' : Number(selected.delta_next).toPrecision(4);
    elements['delta-median'].textContent = selected.delta_median === '' ? 'Not available' : Number(selected.delta_median).toPrecision(4);
    const scores = [...(scoresByCell.get(selected.cell_id) || [])].sort((a, b) => Number(a.rank) - Number(b.rank));
    elements.scores.replaceChildren();
    if (scores.length < 2) throw new Error(`Insufficient candidate scores for ${selected.cell_id}.`);
    for (const item of scores) {
      const row = makeElement('button', 'score-row');
      row.type = 'button';
      row.setAttribute('aria-label', `Compare ${item.label} markers with the leading candidate`);
      const activeCandidate = selectedCandidate && selectedCandidate !== scores[0].label ? selectedCandidate : scores[1].label;
      row.setAttribute('aria-pressed', String(item.label === activeCandidate));
      if (item.label === scores[0].label) {
        row.disabled = true;
        row.setAttribute('aria-label', `${item.label} is the leading initial candidate`);
      } else {
        row.addEventListener('click', () => { selectedCandidate = item.label; renderEvidence(); elements.scores.querySelector('[aria-pressed="true"]')?.focus(); });
      }
      const title = makeElement('div', 'score-title');
      title.append(makeElement('span', '', `${item.rank}. ${item.label}`), makeElement('strong', '', Number(item.score).toFixed(3)));
      const track = makeElement('div', 'score-track');
      const marker = makeElement('span', 'score-marker');
      marker.style.left = `${Math.max(0, Math.min(100, (Number(item.score) + 1) * 50))}%`;
      track.append(marker);
      row.append(title, track);
      elements.scores.append(row);
    }
    const top = scores[0];
    const next = scores.find((item) => item.label === selectedCandidate && item.label !== top.label) || scores[1];
    elements['marker-description'].textContent = `Reference-derived markers for ${top.label} versus ${next.label}, and the reverse comparison. Values are raw counts in this query cell.`;
    elements.markers.replaceChildren();
    const expression = new Map((expressionByCell.get(selected.cell_id) || []).map((row) => [row.gene_id, row.count]));
    for (const [label, comparison] of [[top.label, next.label], [next.label, top.label]]) {
      const group = makeElement('div', 'marker-group');
      group.append(makeElement('h4', '', `${label} vs ${comparison}`));
      const genes = (markersByPair.get(`${label}\0${comparison}`) || []).slice(0, 5);
      if (!genes.length) group.append(makeElement('p', 'subtle', 'No pairwise markers available.'));
      for (const gene of genes) {
        const row = makeElement('div', 'gene-row');
        const value = Number(expression.get(gene.gene_id));
        const max = Math.max(1, ...genes.map((item) => Number(expression.get(item.gene_id)) || 0));
        const bar = makeElement('span', 'gene-bar');
        const fill = makeElement('span', 'gene-bar-fill'); fill.style.width = `${Number.isFinite(value) ? value / max * 100 : 0}%`; bar.append(fill);
        row.append(makeElement('span', 'gene-name', gene.gene_id), bar, makeElement('span', 'gene-count', expression.has(gene.gene_id) ? `${expression.get(gene.gene_id)} counts` : 'Not available'));
        group.append(row);
        const referenceA = referenceExpression.get(`${label}\0${gene.gene_id}`);
        const referenceB = referenceExpression.get(`${comparison}\0${gene.gene_id}`);
        group.append(makeElement('p', 'reference-means', `Reference mean log expression: ${label} ${referenceA ? Number(referenceA.mean_log_expression).toFixed(2) : 'NA'} · ${comparison} ${referenceB ? Number(referenceB.mean_log_expression).toFixed(2) : 'NA'}`));
      }
      elements.markers.append(group);
    }
  }
  function renderAll() { renderList(); renderOverview(); renderDistribution(); renderEvidence(); }
  search.addEventListener('input', renderAll);
  labelFilter.addEventListener('change', renderAll);
  prunedFilter.addEventListener('change', renderAll);
  document.querySelector('#clear-bin').addEventListener('click', () => { selectedBin = null; renderAll(); });
  function resetFilters() { search.value = ''; labelFilter.value = 'all'; prunedFilter.checked = false; selectedBin = null; }
  document.querySelector('#reset-filters').addEventListener('click', () => { resetFilters(); renderAll(); });
  function selectExample(mode) {
    const available = cells.filter((cell) => cell.delta_next !== '' && Number.isFinite(Number(cell.delta_next)));
    const sorted = [...available].sort((a, b) => Number(a.delta_next) - Number(b.delta_next) || a.cell_id.localeCompare(b.cell_id));
    const example = mode === 'smallest' ? sorted[0] : mode === 'largest' ? sorted[sorted.length - 1] : [...cells].filter((cell) => cell.is_pruned === 'TRUE').sort((a, b) => a.cell_id.localeCompare(b.cell_id))[0];
    if (example) { resetFilters(); selectCell(example); }
  }
  document.querySelector('#example-clear').addEventListener('click', () => selectExample('largest'));
  document.querySelector('#example-ambiguous').addEventListener('click', () => selectExample('smallest'));
  document.querySelector('#example-pruned').addEventListener('click', () => selectExample('pruned'));
  renderStep();
  renderAll();
  status.hidden = true;
  content.hidden = false;
  document.querySelector('#visual-lab').hidden = false;
}

Promise.all(['cells', 'scores', 'markers', 'marker_expression', 'pca', 'reference_expression'].map(loadTable))
  .then((tables) => start(...tables))
  .catch((error) => {
    status.classList.add('error');
    status.textContent = `Unable to load the local example: ${error.message} Follow the complete R preparation workflow in README.md, then serve this directory over HTTP.`;
  });
