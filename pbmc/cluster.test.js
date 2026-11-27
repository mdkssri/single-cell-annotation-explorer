"use strict";

const assert = require("node:assert/strict");
const { readFileSync } = require("node:fs");
const { resolve } = require("node:path");
const { runInNewContext } = require("node:vm");
const { test } = require("node:test");
const { pca, initialize, evaluate, expectation, maximize, startFit, stepFit } = require("./dist/model");

const source = readFileSync(resolve(__dirname, "dist/pbmc.js"), "utf8");
const data = runInNewContext(source + "\nPBMC;");
const fitted = pca(data);

test("browser PCA matches independent R prcomp on the published teaching matrix", () => {
  assert.equal(fitted.scores.length, 256);
  assert.equal(data.genes.length, 160);
  const r = [.421986927239819, .173803362971707, .050532868177231, .045452694183480];
  assert.ok(fitted.fractions.every((value, i) => Math.abs(value - r[i]) < 1e-8));
});

test("each EM step preserves normalized probabilities and the M-step increases fit", () => {
  let model = initialize(fitted.scores, 4, [1, 1, 1, 1]);
  for (let i = 0; i < 12; i++) {
    const e = expectation(fitted.scores, model);
    assert.ok(e.rows.every(row => Math.abs(row.values.reduce((total, value) => total + value.probability, 0) - 1) < 1e-10));
    model = maximize(fitted.scores, model, e);
    assert.ok(expectation(fitted.scores, model).likelihood >= e.likelihood - 1e-8);
  }
});

test("weighted settings and cluster count change the fitted model without changing the PCA coordinates", () => {
  const initial = initialize(fitted.scores, 4, [1, 1, 1, 1]);
  const weighted = initialize(fitted.scores, 3, [1, .25, 0, 2]);
  const point = fitted.scores[31];
  const result = evaluate(point, weighted);
  assert.equal(initial.groups.length, 4);
  assert.equal(result.values.length, 3);
  assert.ok(result.values.every(group => group.distance[2] === 0));
  assert.ok(result.values.every(group => Number.isFinite(group.weighted)));
  assert.equal(result.winner, result.values.indexOf(result.values.reduce((a, b) => a.probability >= b.probability ? a : b)));
});

test("a cell experiment starts unscored and advances independently from a main fit", () => {
  const main = startFit(fitted.scores, 4, [1, 1, 1, 1]);
  stepFit(fitted.scores, main);
  const detail = startFit(fitted.scores, 3, [1, 1, 0, 1]);
  assert.equal(detail.lastE, null);
  stepFit(fitted.scores, detail);
  const priorVoteModel = detail.voteModel;
  const priorVotes = detail.lastE.rows[31];
  stepFit(fitted.scores, detail);
  assert.equal(detail.next, "E");
  assert.equal(detail.rounds, 1);
  assert.equal(detail.voteModel, priorVoteModel);
  assert.equal(detail.lastE.rows[31], priorVotes);
  assert.notEqual(detail.model, detail.voteModel);
  assert.equal(main.rounds, 0);
  assert.equal(main.model.groups.length, 4);
  const restarted = startFit(fitted.scores, 3, [1, 1, 0, 1]);
  assert.equal(restarted.lastE, null);
  stepFit(fitted.scores, restarted);
  assert.deepEqual(restarted.lastE.rows[31], priorVotes);
});

test("recorded M-step operands reproduce the updated parameters and stopping decision", () => {
  const fit = startFit(fitted.scores, 4, [.5, 1, 1, 1]);
  let converged = false;
  for (let i = 0; i < 100 && !converged; i++) {
    stepFit(fitted.scores, fit);
    ({ converged } = stepFit(fitted.scores, fit));
  }
  assert.equal(fit.rounds, 94);
  assert.equal(fit.lastM.converged, true);
  assert.equal(fit.lastM.delta, fit.history.at(-1) - fit.history.at(-2));
  const { mass, moments, scatter } = fit.lastM.updates[0];
  const group = fit.model.groups[0];
  assert.equal(group.prior, mass / fitted.scores.length);
  assert.deepEqual(group.center, moments.map(value => value / mass));
  assert.equal(group.variance, Math.max(fit.model.floor, scatter / (mass * fit.model.dimensions)));
});
