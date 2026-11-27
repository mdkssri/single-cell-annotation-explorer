"use strict";

const assert = require("node:assert/strict");
const { readFileSync } = require("node:fs");
const { join } = require("node:path");
const { test } = require("node:test");
const { runInNewContext } = require("node:vm");
const model = require("./dist/marker-model");
const { pca, evaluate } = require("./dist/model");

const source = readFileSync(join(__dirname, "dist/pbmc.js"), "utf8");
const { cells, totals, genes } = runInNewContext(source + "\nPBMC;");
const counts = { cells, totals, genes };

let fitted;
test("count-based EM converges with normalized responsibilities and ordered marker rates", () => {
  const data = model.prepare(counts);
  let parameters = model.initialize(data);
  let result = model.expectation(data, parameters);
  for (let i = 0; i < 100; i++) {
    parameters = model.maximize(data, result);
    assert.ok(parameters.rates.every(rate => rate.low > 0 && rate.high >= rate.low));
    const updated = model.expectation(data, parameters);
    const delta = updated.logLikelihood - result.logLikelihood;
    assert.ok(delta >= -1e-8);
    assert.ok(updated.rows.every(row =>
      Math.abs(row.probabilities.reduce((sum, value) => sum + value, 0) - 1) < 1e-10));
    result = updated;
    if (delta < model.convergence) {
      fitted = result;
      break;
    }
  }
  assert.ok(fitted);
});

test("threshold controls ambiguity and preserves the separate background category", () => {
  const groups = threshold => fitted.rows.map(row => model.assignment(row, threshold));
  const standard = groups(model.defaultThreshold);
  const strict = groups(.99);
  assert.ok(strict.filter(call => call.kind === "ambiguous").length > standard.filter(call => call.kind === "ambiguous").length);
  const isBackground = call => call.kind === "unassigned";
  assert.ok(standard.some(isBackground));
  assert.deepEqual(strict.map(isBackground), standard.map(isBackground));
  assert.ok([0, 1, 2].every(type => standard.some(call => call.kind === "anchor" && call.top.k === type)));
});

test("the optional review uses non-marker columns and evaluates only ambiguous cells", () => {
  const drop = new Set(model.mark.flatMap(type => type.genes));
  const check = model.review(counts, fitted, model.defaultThreshold, (input, dimensions) => {
    assert.deepEqual(Array.from(input.cells[0]), Array.from(cells[0]).filter((_, i) => !drop.has(genes[i])));
    return pca(input, dimensions);
  }, evaluate);
  assert.ok(check.sizes.every(n => n >= model.minimumAnchors));
  assert.ok(check.rows.some(row => row?.corroborates));
  assert.ok(check.rows.every((row, i) => Boolean(row) === (model.assignment(fitted.rows[i], model.defaultThreshold).kind === "ambiguous")));
});
