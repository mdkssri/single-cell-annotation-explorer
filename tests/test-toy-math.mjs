import assert from 'node:assert/strict';
import { getAverageRanks, getQuantile, getSpearman } from '../toy-math.js';

assert.deepEqual(getAverageRanks([5, 1, 5, 3]), [3.5, 1, 3.5, 2]);
assert.equal(getSpearman([1, 2, 3], [4, 5, 6]), 1);
assert.equal(getSpearman([1, 2, 3], [6, 5, 4]), -1);
assert.equal(getSpearman([1, 1, 1], [1, 2, 3]), null);
assert.equal(getSpearman([1], [1]), null);
assert.ok(Math.abs(getQuantile([0.1, 0.4, 0.9], 0.8) - 0.7) < 1e-12);
assert.equal(getQuantile([1, 2, 3], 0), 1);
assert.equal(getQuantile([1, 2, 3], 1), 3);
assert.equal(getQuantile([], 0.8), null);
assert.equal(getQuantile([1, 2], NaN), null);
assert.equal(getQuantile([1, 2], -0.1), null);
assert.ok(Math.abs(getSpearman([1, 1, 2], [1, 2, 3]) - Math.sqrt(3) / 2) < 1e-12);
assert.deepEqual(getAverageRanks([0, 0, 0]), [2, 2, 2]);
console.log('Toy rank, Spearman and quantile tests passed.');
