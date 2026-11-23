export function getAverageRanks(values) {
  const sorted = values.map((value, index) => ({ value, index })).sort((a, b) => a.value - b.value);
  const ranks = new Array(values.length);
  for (let start = 0; start < sorted.length;) {
    let end = start + 1;
    while (end < sorted.length && sorted[end].value === sorted[start].value) end++;
    const rank = (start + 1 + end) / 2;
    for (let index = start; index < end; index++) ranks[sorted[index].index] = rank;
    start = end;
  }
  return ranks;
}

export function getSpearman(left, right) {
  if (left.length !== right.length || left.length < 2) return null;
  const a = getAverageRanks(left);
  const b = getAverageRanks(right);
  const meanA = a.reduce((sum, value) => sum + value, 0) / a.length;
  const meanB = b.reduce((sum, value) => sum + value, 0) / b.length;
  let numerator = 0;
  let squareA = 0;
  let squareB = 0;
  for (let index = 0; index < a.length; index++) {
    const da = a[index] - meanA;
    const db = b[index] - meanB;
    numerator += da * db;
    squareA += da * da;
    squareB += db * db;
  }
  if (squareA === 0 || squareB === 0) return null;
  return numerator / Math.sqrt(squareA * squareB);
}

export function getQuantile(values, probability) {
  if (!values.length || !Number.isFinite(probability) || probability < 0 || probability > 1 || values.some((value) => !Number.isFinite(value))) return null;
  const sorted = [...values].sort((a, b) => a - b);
  const position = (sorted.length - 1) * probability;
  const lower = Math.floor(position);
  const fraction = position - lower;
  return sorted[lower] + fraction * (sorted[Math.min(lower + 1, sorted.length - 1)] - sorted[lower]);
}
