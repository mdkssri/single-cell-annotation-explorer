import { getAverageRanks, getQuantile, getSpearman } from './toy-math.js';

const genes = ['Gene A', 'Gene B', 'Gene C', 'Gene D', 'Gene E'];
const profiles = [
  { id: 'B reference 1', label: 'B-like', values: [9, 7, 2, 1, 4] },
  { id: 'B reference 2', label: 'B-like', values: [8, 6, 3, 2, 5] },
  { id: 'B reference 3', label: 'B-like', values: [7, 9, 1, 2, 4] },
  { id: 'T reference 1', label: 'T-like', values: [1, 2, 9, 7, 4] },
  { id: 'T reference 2', label: 'T-like', values: [2, 3, 8, 6, 5] },
  { id: 'T reference 3', label: 'T-like', values: [2, 1, 7, 9, 4] }
];
const presets = { b: [8, 7, 2, 1, 4], t: [1, 2, 8, 7, 4], mixed: [7, 2, 7, 2, 4] };
const state = { values: [...presets.b], enabled: genes.map(() => true) };
const geneContainer = document.querySelector('#toy-genes');
const referenceSelect = document.querySelector('#toy-reference');
const quantileInput = document.querySelector('#toy-quantile');
const presetSelect = document.querySelector('#toy-preset');
const warning = document.querySelector('#toy-warning');
const labelDetail = document.querySelector('#toy-label-detail');

function makeElement(tag, className, text) {
  const element = document.createElement(tag);
  if (className) element.className = className;
  if (text !== undefined) element.textContent = text;
  return element;
}

function makeSvgElement(tag, attributes = {}) {
  const element = document.createElementNS('http://www.w3.org/2000/svg', tag);
  for (const [name, value] of Object.entries(attributes)) element.setAttribute(name, value);
  return element;
}

function renderRankScatter(selected, query, referenceValues, correlation) {
  const plot = document.querySelector('#rank-scatter');
  const calculation = document.querySelector('#rank-calculation');
  plot.replaceChildren();
  const n = selected.length;
  if (n < 2) {
    calculation.textContent = 'Select at least two genes to see the rank calculation.';
    return;
  }
  const queryRanks = getAverageRanks(query);
  const referenceRanks = getAverageRanks(referenceValues);
  const meanQuery = queryRanks.reduce((sum, value) => sum + value, 0) / n;
  const meanReference = referenceRanks.reduce((sum, value) => sum + value, 0) / n;
  let cross = 0;
  let squareQuery = 0;
  let squareReference = 0;
  for (let index = 0; index < n; index++) {
    const x = queryRanks[index] - meanQuery;
    const y = referenceRanks[index] - meanReference;
    cross += x * y;
    squareQuery += x * x;
    squareReference += y * y;
  }
  const positionX = (rank) => 50 + (rank - 1) / (n - 1) * 260;
  const positionY = (rank) => 220 - (rank - 1) / (n - 1) * 175;
  plot.append(makeSvgElement('line', { x1: 50, y1: 220, x2: 310, y2: 220, class: 'axis-line' }));
  plot.append(makeSvgElement('line', { x1: 50, y1: 45, x2: 50, y2: 220, class: 'axis-line' }));
  plot.append(makeSvgElement('line', { x1: 50, y1: 220, x2: 310, y2: 45, class: 'math-diagonal' }));
  selected.forEach((geneIndex, index) => {
    const circle = makeSvgElement('circle', { cx: positionX(queryRanks[index]), cy: positionY(referenceRanks[index]), r: 6, class: 'math-point' });
    const title = makeSvgElement('title');
    title.textContent = `${genes[geneIndex]}: query rank ${queryRanks[index]}, reference rank ${referenceRanks[index]}`;
    circle.append(title);
    const label = makeSvgElement('text', { x: positionX(queryRanks[index]) + 8, y: positionY(referenceRanks[index]) - 7, class: 'math-point-label' });
    label.textContent = genes[geneIndex].replace('Gene ', '');
    plot.append(circle, label);
  });
  const xLabel = makeSvgElement('text', { x: 180, y: 255, 'text-anchor': 'middle', class: 'axis-label' }); xLabel.textContent = 'Query rank';
  const yLabel = makeSvgElement('text', { x: 15, y: 135, transform: 'rotate(-90 15 135)', 'text-anchor': 'middle', class: 'axis-label' }); yLabel.textContent = 'Reference rank';
  plot.append(xLabel, yLabel);
  calculation.textContent = correlation === null
    ? 'Undefined: one rank vector has zero variance, so the denominator is zero.'
    : `Σ cross-products = ${cross.toFixed(2)}; rank sums of squares = ${squareQuery.toFixed(2)} and ${squareReference.toFixed(2)}; ρ = ${correlation.toFixed(3)}.`;
}

function renderQuantileDetail(correlations, probability) {
  const chart = document.querySelector('#quantile-chart');
  const explanation = document.querySelector('#quantile-calculation');
  chart.replaceChildren();
  const group = correlations.filter((profile) => profile.label === labelDetail.value);
  if (group.some((profile) => profile.correlation === null)) {
    explanation.textContent = 'The quantile is unavailable because at least one reference correlation is undefined.';
    return;
  }
  const sorted = group.map((profile) => profile.correlation).sort((a, b) => a - b);
  const position = (sorted.length - 1) * probability;
  const lower = Math.floor(position);
  const upper = Math.ceil(position);
  const fraction = position - lower;
  sorted.forEach((value, index) => {
    const item = makeElement('div', `quantile-item${index === lower || index === upper ? ' quantile-used' : ''}`);
    item.append(makeElement('span', '', `x${index}`), makeElement('strong', '', value.toFixed(3)));
    chart.append(item);
  });
  const result = getQuantile(sorted, probability);
  explanation.textContent = `h = (${sorted.length} − 1) × ${probability.toFixed(1)} = ${position.toFixed(2)}; lower = x${lower}, upper = x${upper}; fraction = ${fraction.toFixed(2)}; score = ${result.toFixed(3)}.`;
}

for (const profile of profiles) {
  const option = makeElement('option', '', profile.id);
  option.value = profile.id;
  referenceSelect.append(option);
}

for (const [index, gene] of genes.entries()) {
  const row = makeElement('div', 'toy-gene');
  const toggleLabel = makeElement('label', 'toy-gene-label');
  const toggle = makeElement('input');
  toggle.type = 'checkbox';
  toggle.checked = true;
  toggle.setAttribute('aria-label', `Include ${gene}`);
  toggle.addEventListener('change', () => { state.enabled[index] = toggle.checked; render(); });
  toggleLabel.append(toggle, document.createTextNode(gene));
  const slider = makeElement('input');
  slider.type = 'range'; slider.min = '0'; slider.max = '10'; slider.step = '1'; slider.value = state.values[index];
  slider.setAttribute('aria-label', `${gene} query expression`);
  const output = makeElement('output', 'toy-gene-value', slider.value);
  slider.addEventListener('input', () => { state.values[index] = Number(slider.value); output.textContent = slider.value; render(); });
  row.append(toggleLabel, slider, output);
  geneContainer.append(row);
}

function render() {
  const selected = genes.map((_, index) => index).filter((index) => state.enabled[index]);
  const query = selected.map((index) => state.values[index]);
  const reference = profiles.find((profile) => profile.id === referenceSelect.value);
  const referenceValues = selected.map((index) => reference.values[index]);
  const correlations = profiles.map((profile) => ({ ...profile, correlation: getSpearman(query, selected.map((index) => profile.values[index])) }));
  const undefinedCorrelation = selected.length < 2 || correlations.some((profile) => profile.correlation === null);
  warning.hidden = !undefinedCorrelation;
  warning.textContent = selected.length < 2 ? 'Select at least two genes to calculate a correlation.' : 'Correlation is undefined when a selected profile has no rank variation. Change the gene selection or query values.';
  const ranksContainer = document.querySelector('#toy-ranks');
  ranksContainer.replaceChildren();
  if (selected.length) {
    const queryRanks = getAverageRanks(query);
    const referenceRanks = getAverageRanks(referenceValues);
    const header = makeElement('div', 'toy-rank-row toy-rank-header');
    for (const value of ['Gene', 'Query value → rank', 'Ref. value → rank']) header.append(makeElement('span', '', value));
    ranksContainer.append(header);
    selected.forEach((geneIndex, index) => {
      const row = makeElement('div', 'toy-rank-row');
      row.append(makeElement('strong', '', genes[geneIndex]), makeElement('span', '', `${query[index]} → ${queryRanks[index]}`), makeElement('span', '', `${referenceValues[index]} → ${referenceRanks[index]}`));
      ranksContainer.append(row);
    });
  }
  const chosen = correlations.find((profile) => profile.id === reference.id);
  document.querySelector('#toy-rank-result').textContent = `Spearman correlation: ${chosen.correlation === null ? 'undefined' : chosen.correlation.toFixed(3)}`;
  renderRankScatter(selected, query, referenceValues, chosen.correlation);
  renderQuantileDetail(correlations, Number(quantileInput.value));
  const correlationContainer = document.querySelector('#toy-correlations');
  const scoreContainer = document.querySelector('#toy-scores');
  correlationContainer.replaceChildren();
  scoreContainer.replaceChildren();
  const labelScores = [];
  for (const label of ['B-like', 'T-like']) {
    const group = correlations.filter((profile) => profile.label === label);
    const row = makeElement('div', 'toy-correlation-row');
    row.append(makeElement('strong', '', label), makeElement('span', '', group.map((profile) => profile.correlation === null ? 'undefined' : profile.correlation.toFixed(3)).join('  ·  ')));
    correlationContainer.append(row);
    const score = group.some((profile) => profile.correlation === null) ? null : getQuantile(group.map((profile) => profile.correlation), Number(quantileInput.value));
    labelScores.push({ label, score });
    const scoreRow = makeElement('div', 'toy-score');
    scoreRow.append(makeElement('span', '', `${label} quantile score`), makeElement('strong', '', score === null ? 'undefined' : score.toFixed(3)));
    scoreContainer.append(scoreRow);
    if (score !== null) {
      const track = makeElement('div', 'toy-score-track');
      const point = makeElement('span', 'toy-score-point');
      point.style.left = `${Math.max(0, Math.min(100, (score + 1) * 50))}%`;
      track.append(point);
      scoreContainer.append(track);
    }
  }
  if (labelScores.every((item) => item.score !== null)) {
    const difference = labelScores[0].score - labelScores[1].score;
    const text = Math.abs(difference) < 1e-12 ? 'Initial scores tie; no unique leader.' : `Initial leader: ${difference > 0 ? labelScores[0].label : labelScores[1].label}.`;
    scoreContainer.append(makeElement('p', 'toy-leader', `${text} Initial score gap: ${Math.abs(difference).toFixed(3)} (not delta.next).`));
  }
  document.querySelector('#toy-quantile-value').textContent = Number(quantileInput.value).toFixed(1);
}

function resetExample() {
  state.values = [...presets[presetSelect.value]];
  state.enabled = genes.map(() => true);
  for (const [index, row] of [...geneContainer.children].entries()) {
    row.querySelector('input[type="checkbox"]').checked = true;
    row.querySelector('input[type="range"]').value = state.values[index];
    row.querySelector('output').textContent = state.values[index];
  }
  quantileInput.value = '0.8';
  referenceSelect.selectedIndex = 0;
  render();
}

presetSelect.addEventListener('change', resetExample);
document.querySelector('#toy-reset').addEventListener('click', resetExample);
referenceSelect.addEventListener('change', render);
quantileInput.addEventListener('input', render);
labelDetail.addEventListener('change', render);
render();
