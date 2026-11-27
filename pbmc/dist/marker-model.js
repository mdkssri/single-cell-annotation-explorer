"use strict";

const MarkerModel = (() => {
  const mark = [
    { name: "T cell", genes: ["CD3D", "CD3E"] },
    { name: "B cell", genes: ["MS4A1", "CD79A"] },
    { name: "Monocyte", genes: ["LST1", "S100A9", "LYZ"] },
    { name: "Unassigned", genes: [] }
  ];
  const size = 2;
  const minimum = 1e-9;
  const defaultThreshold = .95;
  const convergence = 1e-5;
  const minimumAnchors = 5;
  const sum = values => values.reduce((a, b) => a + b, 0);
  const logKernel = (count, mean) => count * Math.log(mean) - (count + size) * Math.log(size + mean);

  function prepare(data) {
    const depth = data.totals.map(value => value / 10000);
    const genes = mark.flatMap(({ genes }, k) => genes.map(name => ({ name, k, column: data.genes.indexOf(name) })));
    if (genes.some(gene => gene.column < 0)) throw Error("A marker is missing from the expression panel.");
    for (const gene of genes) {
      gene.counts = data.cells.map(cell => cell[gene.column]);
      gene.upper = Math.max(minimum, ...gene.counts.map((count, i) => count / depth[i]));
    }
    return { genes, depth, count: depth.length };
  }

  function initialize(data) {
    const rates = data.genes.map(({ counts, upper }) => {
      const baseline = Math.max(minimum, sum(counts) / sum(data.depth) / 2);
      return { low: baseline, high: Math.max(baseline, Math.min(upper, baseline * 5)) };
    });
    return { rates, priors: mark.map(() => 1 / mark.length) };
  }

  function row(data, model, index) {
    const depth = data.depth[index];
    const values = model.priors.map((prior, k) => {
      const terms = data.genes.map((gene, j) => {
        const rate = gene.k === k ? model.rates[j].high : model.rates[j].low;
        const mean = depth * rate;
        const count = gene.counts[index];
        return { mean, count, score: logKernel(count, mean) };
      });
      const logPrior = Math.log(prior);
      const logScore = logPrior + sum(terms.map(term => term.score));
      return { terms, logPrior, logScore };
    });
    const largest = Math.max(...values.map(value => value.logScore));
    const total = largest + Math.log(sum(values.map(value => Math.exp(value.logScore - largest))));
    const probabilities = values.map(value => Math.exp(value.logScore - total));
    return { total, values, probabilities };
  }

  function expectation(data, model) {
    const rows = Array.from({ length: data.count }, (_, i) => row(data, model, i));
    return { rows, logLikelihood: sum(rows.map(entry => entry.total)) };
  }

  function fitRate(data, gene, weights) {
    if (sum(weights) <= 0) return minimum;
    function derivative(rate) {
      let value = 0;
      for (let i = 0; i < data.count; i++) {
        const mean = data.depth[i] * rate;
        value += weights[i] * (gene.counts[i] - mean) / (size + mean);
      }
      return value;
    }
    if (derivative(minimum) <= 0) return minimum;
    let low = minimum, high = gene.upper;
    for (let j = 0; j < 40; j++) {
      const middle = (high + low) / 2;
      if (derivative(middle) > 0) low = middle;
      else high = middle;
    }
    return (high + low) / 2;
  }

  function maximize(data, prediction) {
    const columns = mark.map((_, k) => prediction.rows.map(entry => entry.probabilities[k]));
    const priors = columns.map(column => sum(column) / data.count);
    const rates = data.genes.map(gene => {
      const highWeight = columns[gene.k];
      const lowWeight = highWeight.map(value => Math.max(0, 1 - value));
      let high = fitRate(data, gene, highWeight);
      let low = fitRate(data, gene, lowWeight);
      if (high < low) high = low = fitRate(data, gene, lowWeight.map(() => 1));
      return { high, low };
    });
    return { rates, priors };
  }

  function assignment(entry, threshold) {
    const order = entry.probabilities.map((p, k) => ({ p, k })).sort((a, b) => b.p - a.p);
    const [top, next] = order;
    return { top, next, kind: top.k === 3 ? "unassigned" : top.p >= threshold ? "anchor" : "ambiguous" };
  }

  function review(counts, prediction, threshold, pca, evaluate) {
    const used = new Set(mark.flatMap(item => item.genes));
    const includes = counts.genes.map((gene, i) => used.has(gene) ? -1 : i).filter(i => i >= 0);
    const cells = counts.cells.map(cell => includes.map(i => cell[i]));
    const { scores } = pca({ cells, totals: counts.totals }, 4);
    const byType = [[], [], []];
    const assignments = prediction.rows.map(entry => assignment(entry, threshold));
    assignments.forEach(({ top, kind }, i) => { if (kind === "anchor") byType[top.k].push(scores[i]); });
    const sizes = byType.map(group => group.length);
    if (sizes.some(n => n < minimumAnchors)) return { sizes, rows: null };
    const groups = byType.map(group => {
      const center = [0, 1, 2, 3].map(j => sum(group.map(point => point[j])) / group.length);
      const distances = group.map(point => sum(point.map((value, j) => (value - center[j]) ** 2)));
      const variance = Math.max(1e-5, sum(distances) / (4 * group.length));
      const limit = distances.sort((a, b) => a - b)[Math.ceil(.95 * distances.length) - 1];
      return { center, variance, prior: 1 / 3, limit };
    });
    const fit = { weights: [1, 1, 1, 1], groups, dimensions: 4 };
    return { sizes, rows: scores.map((point, i) => {
      if (assignments[i].kind !== "ambiguous") return null;
      const result = evaluate(point, fit);
      const supported = sum(result.values[result.winner].distance) <= groups[result.winner].limit;
      return { winner: result.winner, probability: result.probability, supported,
        corroborates: supported && result.probability >= .9 && result.winner === assignments[i].top.k };
    }) };
  }

  return { mark, size, defaultThreshold, convergence, minimumAnchors, prepare, initialize, expectation, maximize, assignment, review };
})();

if (typeof module !== "undefined") module.exports = MarkerModel;
