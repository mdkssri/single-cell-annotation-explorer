"use strict";

const ClusterModel = (() => {
  const sum = values => values.reduce((total, value) => total + value, 0);
  const tolerance = 1e-5;
  const normalize = (counts, total) => counts.map(count => Math.log1p(10000 * count / total));
  const dot = (a, b) => sum(a.map((value, i) => value * b[i]));

  function pca(data, dimensions = 4) {
    const rows = data.cells.map((counts, i) => normalize(counts, data.totals[i]));
    const width = rows[0].length;
    const means = new Float64Array(width);
    for (const row of rows) for (let j = 0; j < width; j++) means[j] += row[j] / rows.length;
    const centered = rows.map(row => row.map((value, j) => value - means[j]));
    const cov = new Float64Array(width * width);
    for (const row of centered) {
      for (let j = 0; j < width; j++) {
        const offset = j * width;
        for (let k = 0; k <= j; k++) cov[offset + k] += row[j] * row[k];
      }
    }
    const vectors = new Float64Array(width * width);
    for (let j = 0; j < width; j++) {
      vectors[j * width + j] = 1;
      for (let k = 0; k <= j; k++) {
        const value = cov[j * width + k] / (rows.length - 1);
        cov[j * width + k] = value;
        cov[k * width + j] = value;
      }
    }
    const totalVariance = sum(Array.from(means, (_, j) => cov[j * width + j]));
    const tolerance = Math.max(totalVariance, 1) * 1e-11;
    let converged = false;
    for (let sweep = 0; sweep < 40; sweep++) {
      let maximum = 0;
      for (let j = 0; j < width - 1; j++) {
        for (let k = j + 1; k < width; k++) {
          const off = cov[j * width + k];
          maximum = Math.max(maximum, Math.abs(off));
          if (Math.abs(off) <= tolerance) continue;
          const jj = j * width + j;
          const kk = k * width + k;
          const angle = (cov[kk] - cov[jj]) / (2 * off);
          const t = (angle >= 0 ? 1 : -1) / (Math.abs(angle) + Math.hypot(1, angle));
          const c = 1 / Math.hypot(1, t);
          const s = t * c;
          cov[jj] -= t * off;
          cov[kk] += t * off;
          cov[j * width + k] = cov[k * width + j] = 0;
          for (let h = 0; h < width; h++) {
            const hj = h * width + j;
            const hk = h * width + k;
            if (h !== j && h !== k) {
              const a = cov[hj], b = cov[hk];
              cov[j * width + h] = cov[hj] = c * a - s * b;
              cov[k * width + h] = cov[hk] = s * a + c * b;
            }
            const v = vectors[hj], w = vectors[hk];
            vectors[hj] = c * v - s * w;
            vectors[hk] = s * v + c * w;
          }
        }
      }
      if (maximum <= tolerance) { converged = true; break; }
    }
    if (!converged) throw Error("PCA did not converge. Reload and try again.");
    const axes = Array.from(means, (_, j) => j)
      .sort((a, b) => cov[b * width + b] - cov[a * width + a])
      .slice(0, dimensions);
    const loadings = axes.map(axis => {
      const v = Array.from(means, (_, j) => vectors[j * width + axis]);
      const largest = v.reduce((a, b) => Math.abs(a) >= Math.abs(b) ? a : b);
      return largest < 0 ? v.map(value => -value) : v;
    });
    return {
      scores: centered.map(row => loadings.map(axis => dot(row, axis))),
      fractions: axes.map(axis => cov[axis * width + axis] / totalVariance)
    };
  }

  function distanceTerms(point, center, weights) {
    return point.map((value, i) => weights[i] * (value - center[i]) ** 2);
  }

  function initialize(points, k, weights) {
    const dimensions = weights.filter(value => value > 0).length;
    const average = weights.map((_, j) => sum(points.map(point => point[j])) / points.length);
    const variance = Math.max(1e-8, sum(points.map(point => sum(distanceTerms(point, average, weights)))) / points.length / dimensions);
    const centers = [points.reduce((best, point) =>
      sum(distanceTerms(point, average, weights)) < sum(distanceTerms(best, average, weights)) ? point : best)];
    while (centers.length < k) {
      const nearest = point => Math.min(...centers.map(center => sum(distanceTerms(point, center, weights))));
      centers.push(points.reduce((best, point) => nearest(point) > nearest(best) ? point : best));
    }
    return {
      weights: [...weights], dimensions, floor: variance * 1e-4,
      groups: centers.map(center => ({ center: [...center], variance, prior: 1 / k }))
    };
  }

  function evaluate(point, model) {
    const values = model.groups.map(({ center, variance, prior }, index) => {
      const distance = distanceTerms(point, center, model.weights);
      const spread = -model.dimensions / 2 * Math.log(2 * Math.PI * variance);
      const penalty = -sum(distance) / (2 * variance);
      const logPrior = Math.log(prior);
      const logDensity = spread + penalty;
      return { index, distance, spread, penalty, logPrior, logDensity, weighted: logDensity + logPrior };
    });
    const max = Math.max(...values.map(value => value.weighted));
    const logTotal = max + Math.log(sum(values.map(value => Math.exp(value.weighted - max))));
    values.forEach(value => { value.probability = Math.exp(value.weighted - logTotal); });
    const winner = values.reduce((a, b) => a.probability >= b.probability ? a : b);
    return { values, logTotal, winner: winner.index, probability: winner.probability };
  }

  function expectation(points, model) {
    const rows = points.map(point => evaluate(point, model));
    return { rows, likelihood: sum(rows.map(row => row.logTotal)) };
  }

  function maximize(points, model, expectationResult, updates = []) {
    const groups = model.groups.map((group, k) => {
      const votes = expectationResult.rows.map(row => row.values[k].probability);
      const mass = sum(votes);
      if (mass === 0) {
        updates.push({ mass });
        return { ...group, prior: 0 };
      }
      const moments = model.weights.map((_, j) => sum(points.map((point, i) => votes[i] * point[j])));
      const center = moments.map(moment => moment / mass);
      const scatter = sum(points.map((point, i) => votes[i] * sum(distanceTerms(point, center, model.weights))));
      updates.push({ mass, moments, scatter });
      return { center, prior: mass / points.length, variance: Math.max(model.floor, scatter / (mass * model.dimensions)) };
    });
    return { ...model, groups };
  }

  function startFit(points, k, weights) {
    const model = initialize(points, k, weights);
    return { model, lastE: null, lastM: null, voteModel: null, next: "E", rounds: 0, history: [expectation(points, model).likelihood] };
  }

  function stepFit(points, fit) {
    if (fit.next === "E") {
      fit.lastE = expectation(points, fit.model);
      fit.voteModel = fit.model;
      fit.next = "M";
      return { previous: null, converged: false };
    }
    const previous = fit.model.groups;
    const updates = [];
    fit.model = maximize(points, fit.model, fit.lastE, updates);
    const value = expectation(points, fit.model).likelihood;
    const delta = value - fit.history.at(-1);
    fit.history.push(value);
    fit.next = "E";
    fit.rounds += 1;
    const converged = Math.abs(delta) < tolerance;
    fit.lastM = { updates, delta, converged };
    return { previous, converged };
  }

  return { normalize, pca, initialize, evaluate, expectation, maximize, startFit, stepFit, tolerance };
})();

if (typeof module !== "undefined") module.exports = ClusterModel;
