"use strict";

(() => {
  const { normalize, pca, evaluate, startFit, stepFit, tolerance } = ClusterModel;
  const get = id => document.getElementById(id);
  const palette = ["#146d70", "#8653b6", "#be572f", "#4764bd", "#966300", "#b44775", "#408423", "#62657a"];
  const idOf = index => "PBMC-" + String(index + 1).padStart(3, "0");
  const fmt = (value, digits = 3) => Number.isFinite(value) ? value.toFixed(digits).replace("-", "−") : "−∞";
  const percent = value => fmt(value * 100, 1) + "%";
  const state = { pca: null, model: null, lastE: null, next: "E", selected: 0, history: [], rounds: 0, baseline: null, bounds: null };
  let experiment;
  let timer = null;
  let experimentTimer = null;
  let animation = null;
  let plotGeometry;
  let stopReason = "";
  const note = message => { get("status").textContent = message; };
  const label = index => "Cluster " + (index + 1);
  const inputs = Array.from({ length: 4 }, (_, i) => i);
  function weightInputs(prefix) {
    return inputs.map(i => `<label for="${prefix}${i}">PC${i + 1}<input id="${prefix}${i}" type="number" value="1" min="0" max="4" step="0.25" required></label>`).join("");
  }
  get("weights").innerHTML = weightInputs("w");
  get("experiment-weights").innerHTML = weightInputs("experiment-w");
  get("cell-select").innerHTML = PBMC.cells.map((_, i) => `<option value="${i}">${idOf(i)}</option>`).join("");

  function requestedWeights(prefix) {
    const controls = inputs.map(i => get(prefix + i));
    if (controls.some(input => !input.validity.valid)) return null;
    const weights = controls.map(input => input.valueAsNumber);
    return weights.some(value => value > 0) ? weights : null;
  }

  function requested() {
    const weights = requestedWeights("w");
    const k = get("k");
    return weights && k.validity.valid ? { k: k.valueAsNumber, weights } : null;
  }

  function updateControls() {
    const config = requested();
    const model = state.model;
    const pending = model && config && (model.groups.length !== config.k || config.weights.some((w, i) => w !== model.weights[i]));
    get("initialize").disabled = !state.pca || !config;
    get("initialize").textContent = model ? "Apply settings and restart EM" : "Initialize mixture";
    get("settings-note").textContent = !config ? "Choose 2–8 clusters and keep at least one positive PC weight (0–4)."
      : pending ? "Settings edited. Apply them to start a new fit on the same PCA." : "1 uses the original PC scale. 0 omits that axis from the fit.";
    get("e").disabled = !model || state.next !== "E" || Boolean(timer);
    get("m").disabled = !model || state.next !== "M" || Boolean(timer);
    get("play").disabled = !model;
    get("pin").disabled = !model || Boolean(timer);
    get("pin").textContent = state.baseline ? "Replace pinned fit" : "Pin fit";
    get("play").textContent = timer ? "Pause EM" : "Play EM";
    get("phase").textContent = !model ? "Initialize the mixture to start."
      : state.next === "M" ? "E-step complete. Colors reflect the latest cluster probabilities."
      : state.rounds ? "M-step complete. Centers and ellipses moved. Run E to update colors."
      : "Initial centers shown. Run E to color cells.";
    get("likelihood").textContent = state.history.length ? fmt(state.history.at(-1), 2) : "—";
    get("progress").textContent = model ? `${state.rounds} completed M-steps${stopReason ? " · " + stopReason : ""}` : "One E-step and one M-step make a round.";
    const vote = state.lastE?.rows[state.selected];
    get("picked").textContent = vote ? `${label(vote.winner)} · ${percent(vote.probability)} at last E-step` : "Select a point or use the menu.";
    get("legend").innerHTML = model ? model.groups.map((group, i) =>
      `<span><i style="--cluster:${palette[i]}"></i>${label(i)} · ${percent(group.prior)}</span>`).join("") : "";
  }

  function canvasContext(canvas) {
    const { width, height } = canvas.getBoundingClientRect();
    if (!width || !height) return null;
    const ratio = Math.min(window.devicePixelRatio || 1, 2);
    const w = Math.round(width * ratio), h = Math.round(height * ratio);
    if (canvas.width !== w || canvas.height !== h) { canvas.width = w; canvas.height = h; }
    const ctx = canvas.getContext("2d");
    ctx.setTransform(ratio, 0, 0, ratio, 0, 0);
    ctx.clearRect(0, 0, width, height);
    return { ctx, width, height };
  }

  function renderPlot(fit = state, canvas = "plot", centers = fit?.model?.groups ?? []) {
    const view = canvasContext(get(canvas));
    if (!view || !state.pca) return;
    const { ctx, width, height } = view;
    const { scores } = state.pca;
    const { minX, maxX, minY, maxY } = state.bounds;
    const left = 48, right = width - 20, top = 22, bottom = height - 43;
    const unit = Math.min((right - left) / (maxX - minX), (bottom - top) / (maxY - minY));
    const cx = (left + right) / 2, cy = (top + bottom) / 2;
    const x = value => cx + (value - (minX + maxX) / 2) * unit;
    const y = value => cy - (value - (minY + maxY) / 2) * unit;
    const tracing = canvas === "experiment-plot";
    if (!tracing) plotGeometry = { x, y };
    ctx.font = "12px system-ui";
    ctx.strokeStyle = "#dce4df";
    ctx.fillStyle = "#4c625d";
    const span = Math.max(maxX - minX, maxY - minY);
    const tick = span > 64 ? 10 : span > 28 ? 5 : 2;
    for (let n = Math.floor(minX / tick) * tick; x(n) < right; n += tick) {
      if (x(n) < left) continue;
      ctx.beginPath(); ctx.moveTo(x(n), top); ctx.lineTo(x(n), bottom); ctx.stroke();
      ctx.textAlign = "center"; ctx.fillText(n, x(n), bottom + 16);
    }
    for (let n = Math.floor(minY / tick) * tick; y(n) > top; n += tick) {
      if (y(n) > bottom) continue;
      ctx.beginPath(); ctx.moveTo(left, y(n)); ctx.lineTo(right, y(n)); ctx.stroke();
      ctx.textAlign = "right"; ctx.fillText(n, left - 8, y(n) + 4);
    }
    ctx.textAlign = "center"; ctx.fillText("PC1", cx, height - 6);
    ctx.save(); ctx.translate(13, cy); ctx.rotate(-Math.PI / 2); ctx.fillText("PC2", 0, 0); ctx.restore();
    ctx.save();
    ctx.beginPath(); ctx.rect(left, top, right - left, bottom - top); ctx.clip();
    const weights = fit?.model?.weights ?? [];
    centers.forEach(({ center, variance }, i) => {
      if (weights[0] === 0 || weights[1] === 0) return;
      ctx.beginPath();
      ctx.ellipse(x(center[0]), y(center[1]), Math.sqrt(variance / weights[0]) * unit, Math.sqrt(variance / weights[1]) * unit, 0, 0, Math.PI * 2);
      ctx.fillStyle = palette[i]; ctx.globalAlpha = .08; ctx.fill();
      ctx.strokeStyle = palette[i]; ctx.globalAlpha = .65; ctx.lineWidth = 1.5; ctx.stroke();
    });
    if (tracing) {
      const point = scores[state.selected];
      centers.forEach(({ center }, i) => {
        ctx.strokeStyle = palette[i]; ctx.globalAlpha = .65; ctx.lineWidth = 1;
        ctx.setLineDash([3, 4]);
        ctx.beginPath(); ctx.moveTo(x(point[0]), y(point[1])); ctx.lineTo(x(center[0]), y(center[1])); ctx.stroke();
      });
      ctx.setLineDash([]);
    }
    function drawPoint(point, i) {
      const row = fit?.lastE?.rows[i];
      ctx.beginPath(); ctx.arc(x(point[0]), y(point[1]), i === state.selected ? 6 : tracing ? 2.7 : 3.6, 0, Math.PI * 2);
      ctx.fillStyle = row ? palette[row.winner] : "#69817a";
      ctx.globalAlpha = tracing && i !== state.selected ? .26 : row ? .32 + .68 * row.probability : .82;
      ctx.fill();
      if (i === state.selected) { ctx.globalAlpha = 1; ctx.strokeStyle = "#182e2a"; ctx.lineWidth = 2; ctx.stroke(); }
    }
    scores.forEach((point, i) => { if (i !== state.selected) drawPoint(point, i); });
    centers.forEach(({ center }, i) => {
      const px = x(center[0]), py = y(center[1]);
      ctx.globalAlpha = 1;
      ctx.beginPath(); ctx.rect(px - 7, py - 7, 14, 14);
      ctx.fillStyle = palette[i]; ctx.fill();
      ctx.lineWidth = 1.7; ctx.strokeStyle = "white"; ctx.stroke();
      ctx.fillStyle = "white"; ctx.textAlign = "center"; ctx.font = "bold 10px system-ui"; ctx.fillText(i + 1, px, py + 3.5);
    });
    drawPoint(scores[state.selected], state.selected);
    ctx.restore();
  }

  function renderTrend() {
    const view = canvasContext(get("trend"));
    if (!view || state.history.length === 0) return;
    const { ctx, width, height } = view;
    const min = Math.min(...state.history), max = Math.max(...state.history);
    const span = Math.max(.001, max - min);
    const points = state.history.map((value, i) =>
      [8 + i * (width - 16) / Math.max(1, state.history.length - 1), height - 8 - (value - min) * (height - 16) / span]);
    ctx.strokeStyle = "#26727f"; ctx.lineWidth = 2;
    ctx.beginPath(); points.forEach(([x, y], i) => i ? ctx.lineTo(x, y) : ctx.moveTo(x, y)); ctx.stroke();
    const [x, y] = points.at(-1);
    ctx.fillStyle = "#26727f"; ctx.beginPath(); ctx.arc(x, y, 3.5, 0, 2 * Math.PI); ctx.fill();
  }

  function render() {
    updateControls(); renderPlot(); renderTrend(); showMainFormulas();
  }

  function animateCenters(previous, fit = state, canvas = "plot") {
    cancelAnimationFrame(animation);
    if (matchMedia("(prefers-reduced-motion: reduce)").matches) return;
    const start = performance.now();
    const next = fit.model.groups;
    function frame(now) {
      const t = Math.min(1, (now - start) / 350);
      const p = 1 - (1 - t) ** 3;
      renderPlot(fit, canvas, next.map((group, i) => ({
        center: group.center.map((x, j) => previous[i].center[j] + p * (x - previous[i].center[j])),
        variance: previous[i].variance + p * (group.variance - previous[i].variance)
      })));
      if (t < 1) animation = requestAnimationFrame(frame);
    }
    animation = requestAnimationFrame(frame);
  }

  function pause() {
    clearInterval(timer);
    timer = null;
  }

  function runStep() {
    stopReason = "";
    const { previous, converged } = stepFit(state.pca.scores, state);
    render();
    if (previous) animateCenters(previous);
    return converged;
  }

  function select(index) {
    state.selected = index;
    get("cell-select").value = index;
    render();
  }

  get("pca").addEventListener("click", () => {
    get("pca").disabled = true;
    get("pca").textContent = "Computing PCA…";
    requestAnimationFrame(() => requestAnimationFrame(() => {
      const start = performance.now();
      try {
        state.pca = pca(PBMC);
        const values = state.pca.scores;
        const padded = index => {
          const vector = values.map(row => row[index]);
          const min = Math.min(...vector), max = Math.max(...vector);
          const pad = Math.max(.5, (max - min) * .10);
          return [min - pad, max + pad];
        };
        const [minX, maxX] = padded(0), [minY, maxY] = padded(1);
        state.bounds = { minX, maxX, minY, maxY };
        get("plot-empty").hidden = true;
        get("pca").textContent = "PCA computed";
        get("axes").textContent = `PC1 ${percent(state.pca.fractions[0])} · PC2 ${percent(state.pca.fractions[1])} · fit uses up to 4 PCs`;
        note(`PCA computed from ${PBMC.cells.length} cells in ${Math.round(performance.now() - start)} ms. Initialize the mixture to start EM.`);
        render();
      } catch (error) {
        get("pca").disabled = false;
        get("pca").textContent = "Run PCA";
        note(error.message);
      }
    }));
  });

  get("initialize").addEventListener("click", () => {
    const config = requested();
    if (!config || !state.pca) return;
    pause(); cancelAnimationFrame(animation);
    Object.assign(state, startFit(state.pca.scores, config.k, config.weights));
    stopReason = "";
    note(`${config.k} centers initialized. The underlying cell positions are ready for EM.`);
    render();
  });
  get("e").addEventListener("click", () => { runStep(); note("E-step computed every cell’s cluster probabilities."); });
  get("m").addEventListener("click", () => { runStep(); note("M-step updated the centers, spreads and mixing proportions."); });
  get("pin").addEventListener("click", () => {
    state.baseline = state.model;
    updateControls();
    note("Fit pinned. Change the settings, rerun EM, and inspect a cell to compare its scores.");
  });
  get("play").addEventListener("click", () => {
    if (timer) { pause(); stopReason = "paused"; updateControls(); return; }
    const target = state.rounds + 60;
    note("Playing E- and M-steps alternately.");
    timer = setInterval(() => {
      const converged = runStep();
      if (converged || state.rounds >= target) {
        pause();
        runStep();
        stopReason = converged ? "converged" : "iteration limit reached";
        updateControls();
        note(converged ? "EM converged. Inspect a cell to trace its prediction." : "Playback paused after 60 new rounds.");
      }
    }, 510);
    updateControls();
  });
  get("weights").addEventListener("input", updateControls);
  get("k").addEventListener("input", updateControls);
  get("cell-select").addEventListener("change", event => select(Number(event.target.value)));
  get("plot").addEventListener("click", event => {
    if (!plotGeometry) return;
    const box = get("plot").getBoundingClientRect();
    let nearest = -1, distance = 13 ** 2;
    state.pca.scores.forEach((point, index) => {
      const d = (plotGeometry.x(point[0]) - (event.clientX - box.left)) ** 2
        + (plotGeometry.y(point[1]) - (event.clientY - box.top)) ** 2;
      if (d < distance) { nearest = index; distance = d; }
    });
    if (nearest >= 0) { select(nearest); openInspector(); }
  });

  function showGenes() {
    const row = PBMC.cells[state.selected];
    const normalized = normalize(row, PBMC.totals[state.selected]);
    const genes = get("gene-mode").value === "all" ? PBMC.genes : PBMC.markers;
    const index = genes.map(gene => PBMC.genes.indexOf(gene));
    const max = Math.max(1, ...index.map(i => row[i]));
    get("gene-bars").innerHTML = index.map(i => `<div class="gene"><strong>${PBMC.genes[i]}</strong><div><i style="--fill:${row[i] / max * 100}%"></i></div><b>${row[i]}</b><span>log-normalized ${fmt(normalized[i], 2)}</span></div>`).join("");
  }

  function compareFit(title, model, point) {
    const row = evaluate(point, model);
    return `<span><b>${title}</b> · ${label(row.winner)} · ${percent(row.probability)}<br>K=${model.groups.length} · weights ${model.weights.join(", ")}</span>`;
  }

  function equation(name, math) {
    return `<div class="equation"><dt>${name}</dt><dd>${math}</dd></div>`;
  }

  function emFormulas(fit, cellIndex, clusterIndex, showDistance = true) {
    const row = fit.lastE?.rows[cellIndex];
    if (!row) return '<p class="muted">Run E to fill in the formulas with this cell’s current values.</p>';
    const value = row.values[clusterIndex];
    const model = fit.voteModel;
    const group = model.groups[clusterIndex];
    const distance = value.distance.reduce((a, b) => a + b, 0);
    const m = fit.lastM;
    const eRound = fit.next === "M" ? fit.rounds : fit.rounds - 1;
    const eEquations = [
      showDistance ? equation("Weighted squared distance · sum of PC terms", `${value.distance.map(term => fmt(term, 4)).join(" + ")} ≈ ${fmt(distance, 5)}`) : "",
      equation("Gaussian log density", `−(${model.dimensions}/2) × ln(2π × ${fmt(group.variance, 5)}) − ${fmt(distance, 5)} ÷ (2 × ${fmt(group.variance, 5)}) ≈ ${fmt(value.logDensity, 5)}`),
      equation("Log score including mixing proportion", `ln(${fmt(group.prior, 6)}) + (${fmt(value.logDensity, 5)}) ≈ ${fmt(value.weighted, 5)}`),
      equation("Log-sum-exp across all clusters", `ln(Σ exp(ℓₖ)) ≈ ${fmt(row.logTotal, 5)}`),
      equation("Probability: exp(log score − log-sum-exp)", `exp(${fmt(value.weighted, 5)} − (${fmt(row.logTotal, 5)})) ≈ ${percent(value.probability)}`)
    ].join("");
    let mEquations = "";
    if (m) {
      const stats = m.updates[clusterIndex];
      const updated = fit.model.groups[clusterIndex];
      const centers = stats.mass > 0 ? stats.moments.map((moment, j) => model.weights[j] > 0
        ? equation(`New PC${j + 1} center · weighted sum ÷ sum of responsibilities`, `${fmt(moment, 5)} ÷ ${fmt(stats.mass, 5)} ≈ ${fmt(updated.center[j], 5)}`) : "").join("") : "";
      mEquations = `<h4>M-step · round ${fit.rounds} · ${label(clusterIndex)}</h4><dl>
        ${equation("New mixing proportion · sum of responsibilities ÷ cells", `${fmt(stats.mass, 5)} ÷ ${PBMC.cells.length} ≈ ${fmt(updated.prior, 6)}`)}
        ${centers}
        ${stats.mass > 0 ? equation("New variance · weighted scatter ÷ (responsibilities × dimension), with floor",
          `max(${fmt(stats.scatter, 5)} ÷ (${fmt(stats.mass, 5)} × ${fit.model.dimensions}), ${fit.model.floor.toExponential(3)}) ≈ ${fmt(updated.variance, 5)}`) : ""}
        ${equation("Batch log-likelihood change · current minus previous",
          `${fmt(fit.history.at(-1), 10)} − (${fmt(fit.history.at(-2), 10)}) ≈ ${m.delta.toExponential(4)}`)}</dl>
        <p class="convergence ${m.converged ? "converged" : ""}">|Δ| ${m.converged ? "<" : "≥"} ${tolerance.toExponential(0)} · ${m.converged ? "Convergence threshold reached" : "Continue iterating"}</p>`;
    }
    return `<div class="formula-readout"><h4>E-step · ${idOf(cellIndex)} · ${label(clusterIndex)} · after ${eRound} M-steps</h4><dl>${eEquations}</dl>${mEquations}</div>`;
  }

  function showMainFormulas() {
    get("main-formulas").hidden = !state.model;
    if (!state.model) return;
    const select = get("main-formula-cluster");
    const k = state.model.groups.length;
    if (select.options.length !== k + 1) {
      select.innerHTML = '<option value="">Winning cluster</option>' + state.model.groups.map((_, i) =>
        `<option value="${i}">${label(i)}</option>`).join("");
    }
    const chosen = select.value === "" ? state.lastE?.rows[state.selected].winner ?? 0 : Number(select.value);
    get("main-formula-values").innerHTML = emFormulas(state, state.selected, chosen);
  }

  function showComponents() {
    const fit = experiment.fit;
    const result = fit.lastE?.rows[state.selected];
    if (!result) {
      get("cell-scores").innerHTML = '<p class="muted">Run the E-step to reveal the selected cell’s likelihoods and assignment.</p>';
      return;
    }
    const active = result.values[experiment.component];
    const point = state.pca.scores[state.selected];
    const model = fit.voteModel;
    const group = model.groups[active.index];
    const runnerUp = [...result.values].sort((a, b) => b.probability - a.probability)[1];
    const comparison = experiment.saved ?? state.baseline;
    const current = `<span><b>This E-step</b> · ${label(result.winner)} · ${percent(result.probability)}<br>K=${model.groups.length} · weights ${model.weights.join(", ")}</span>`;
    const previous = experiment.previous;
    get("cell-scores").innerHTML = `
      <div class="result" style="--cluster:${palette[result.winner]}"><strong>${label(result.winner)}</strong><span>${percent(result.probability)} model probability</span><span class="muted">Runner-up: ${label(runnerUp.index)} · ${percent(runnerUp.probability)}</span></div>
      ${previous ? `<div class="comparison"><span><b>Previous E-step</b> · ${label(previous.winner)} · ${percent(previous.probability)}</span>${current}</div>` : ""}
      ${comparison ? `<div class="comparison">${compareFit(experiment.saved ? "Previous run" : "Pinned main fit", comparison, point)}${current}<p>Cluster numbers are local to each fit.</p></div>` : ""}
      <div class="component-list">${result.values.map(value => `<button type="button" data-component="${value.index}" aria-pressed="${experiment.component === value.index}" style="--cluster:${palette[value.index]}"><span>${label(value.index)}</span><b>${percent(value.probability)}</b></button>`).join("")}</div>
      <details open><summary>${label(active.index)} likelihood calculation</summary>
        <div class="distance">${active.distance.map((term, i) => `<span>PC${i + 1}<small>${model.weights[i]} × (${fmt(point[i], 2)} − ${fmt(group.center[i], 2)})²</small><b>${fmt(term)}</b></span>`).join("")}</div>
        ${emFormulas(fit, state.selected, active.index, false)}
      </details>`;
  }

  function newExperiment(k, weights) {
    pauseExperiment();
    cancelAnimationFrame(animation);
    experiment = { cell: state.selected, stage: 0, fit: null, component: 0, saved: null, previous: null };
    inputs.forEach(i => { get("experiment-w" + i).value = weights[i]; });
    get("experiment-k").value = k;
  }

  function renderExperiment(message = "") {
    const { fit, stage } = experiment;
    const nextWeights = requestedWeights("experiment-w");
    const weights = fit?.model.weights ?? inputs.map(i => get("experiment-w" + i).valueAsNumber);
    get("experiment-status").setAttribute("aria-live", experimentTimer ? "off" : "polite");
    get("experiment-status").textContent = message;
    ["weight", "k", "fit"].forEach((name, i) => { get(`experiment-${name}-step`).hidden = stage !== i; });
    Array.from(get("experiment-steps").children).forEach((item, i) => {
      if (i === stage) item.setAttribute("aria-current", "step");
      else item.removeAttribute("aria-current");
    });
    get("experiment-next").disabled = !state.pca || !nextWeights;
    get("experiment-initialize").disabled = !get("experiment-k").validity.valid;
    get("coordinates").innerHTML = state.pca ? state.pca.scores[state.selected].map((value, j) =>
      `<span class="${weights[j] === 0 ? "inactive-axis" : ""}">PC${j + 1} <b>${fmt(value)}</b><small>${weights[j] === 0 ? "omitted" : Number.isFinite(weights[j]) ? "weight " + weights[j] : "enter weight"}</small></span>`).join("") : "<p class='muted'>Run PCA in the expression map to unlock this experiment.</p>";
    get("experiment-plot").parentElement.hidden = !state.pca;
    renderPlot(fit, "experiment-plot");
    if (!fit) return;
    get("experiment-e").disabled = fit.next !== "E" || Boolean(experimentTimer);
    get("experiment-m").disabled = fit.next !== "M" || Boolean(experimentTimer);
    get("experiment-auto").disabled = Boolean(fit.lastM?.converged) && !experimentTimer;
    get("experiment-auto").textContent = experimentTimer ? "Pause auto-step" : "Auto-step EM";
    get("experiment-round").textContent = `${fit.rounds} completed M-steps`;
    get("experiment-likelihood").textContent = `Batch log-likelihood ${fmt(fit.history.at(-1), 5)}${fit.lastM ? " · Δ " + fit.lastM.delta.toExponential(3) : ""}`;
    get("experiment-phase").textContent = fit.next === "M" ? "E-step: this cell has new probabilities. Run M to update the centers using all cells."
      : fit.rounds ? "M-step: the centers, spreads and proportions have updated. Displayed probabilities are from the previous E-step."
        : "The squares are the initial centers. Run E to score this cell.";
    showComponents();
  }

  function beginExperimentFit(message) {
    const weights = requestedWeights("experiment-w");
    const k = get("experiment-k");
    if (!state.pca || !weights || !k.validity.valid) return;
    pauseExperiment();
    cancelAnimationFrame(animation);
    if (experiment.fit?.voteModel) experiment.saved = experiment.fit.voteModel;
    experiment.fit = startFit(state.pca.scores, k.valueAsNumber, weights);
    experiment.previous = null;
    experiment.stage = 2;
    renderExperiment(message);
    get("experiment-e").focus({ preventScroll: true });
  }

  function pauseExperiment() {
    clearInterval(experimentTimer);
    experimentTimer = null;
  }

  function runExperimentStep(automatic = false) {
    const fit = experiment.fit;
    const eStep = fit.next === "E";
    if (eStep) experiment.previous = fit.lastE?.rows[state.selected] ?? null;
    const { previous, converged } = stepFit(state.pca.scores, fit);
    if (eStep) experiment.component = fit.lastE.rows[state.selected].winner;
    renderExperiment(converged ? `Convergence threshold reached at round ${fit.rounds}.`
      : eStep ? "This E-step has scored the selected cell." : "The M-step updated the mixture using all 256 cells.");
    if (previous) animateCenters(previous, fit, "experiment-plot");
    if (!automatic) get(eStep ? "experiment-m" : "experiment-e").focus({ preventScroll: true });
    return converged;
  }

  function openInspector() {
    pause(); cancelAnimationFrame(animation); updateControls();
    get("lab").hidden = true;
    const detail = get("cell-detail");
    detail.hidden = false;
    if (experiment?.cell !== state.selected) {
      const config = requested() ?? { k: 4, weights: [1, 1, 1, 1] };
      newExperiment(state.model?.groups.length ?? config.k, state.model?.weights ?? config.weights);
    }
    get("detail-title").textContent = idOf(state.selected);
    get("detail-context").textContent = `Donor 1 PBMC · ${PBMC.totals[state.selected].toLocaleString()} total UMI counts across the full library`;
    showGenes(); renderExperiment();
    detail.focus({ preventScroll: true }); detail.scrollIntoView({ block: "start" });
  }

  get("inspect").addEventListener("click", openInspector);
  get("main-formula-cluster").addEventListener("change", showMainFormulas);
  get("gene-mode").addEventListener("change", showGenes);
  get("experiment-weights").addEventListener("input", () => {
    renderExperiment(requestedWeights("experiment-w") ? "" : "Enter four weights from 0 to 4, with at least one positive value.");
  });
  get("experiment-next").addEventListener("click", () => {
    experiment.stage = 1; renderExperiment(); get("experiment-k").focus();
  });
  get("experiment-previous").addEventListener("click", () => {
    experiment.stage = 0; renderExperiment(); get("experiment-w0").focus();
  });
  get("experiment-k").addEventListener("input", () => {
    renderExperiment(get("experiment-k").validity.valid ? "" : "Use a whole cluster count from 2 to 8.");
  });
  get("experiment-initialize").addEventListener("click", () => beginExperimentFit("Initialized at round 0."));
  get("experiment-repeat").addEventListener("click", () => beginExperimentFit("Restarted from round 0 with the same settings."));
  get("experiment-edit").addEventListener("click", () => {
    pauseExperiment();
    cancelAnimationFrame(animation);
    const scored = experiment.fit.voteModel;
    if (scored) experiment.saved = scored;
    experiment.fit = null;
    experiment.stage = 0;
    renderExperiment(scored ? "Edit the settings and initialize a new run. The previous scored fit is saved for comparison."
      : "Edit the settings, then initialize a new run.");
    get("experiment-w0").focus();
  });
  get("experiment-reset").addEventListener("click", () => {
    newExperiment(4, [1, 1, 1, 1]);
    renderExperiment("Defaults restored: four clusters and all weights set to 1.");
    get("experiment-w0").focus();
  });
  get("experiment-e").addEventListener("click", () => runExperimentStep());
  get("experiment-m").addEventListener("click", () => runExperimentStep());
  get("experiment-auto").addEventListener("click", () => {
    if (experimentTimer) {
      pauseExperiment();
      renderExperiment("Auto-step paused. You can continue manually or resume playback.");
      return;
    }
    const limit = experiment.fit.rounds + 300;
    experimentTimer = setInterval(() => {
      const converged = runExperimentStep(true);
      const reachedLimit = experiment.fit.rounds >= limit;
      if (!converged && !reachedLimit) return;
      pauseExperiment();
      if (converged) runExperimentStep(true);
      renderExperiment(converged ? `Converged at round ${experiment.fit.rounds}. Final E-step shown.`
        : "Paused after 300 new rounds. Resume or use the manual controls.");
    }, 300);
    renderExperiment("Playing one E- or M-step at a time.");
  });
  get("cell-scores").addEventListener("click", event => {
    const target = event.target.closest("[data-component]");
    if (!target) return;
    experiment.component = Number(target.dataset.component);
    showComponents();
    get("cell-scores").querySelector('[aria-pressed="true"]').focus({ preventScroll: true });
  });
  get("back").addEventListener("click", () => {
    pauseExperiment();
    cancelAnimationFrame(animation);
    get("cell-detail").hidden = true; get("lab").hidden = false;
    render(); get("inspect").focus({ preventScroll: true });
    get("lab").scrollIntoView({ block: "start" });
  });
  new ResizeObserver(() => { renderPlot(); renderTrend(); }).observe(get("lab"));
  new ResizeObserver(() => renderPlot(experiment?.fit, "experiment-plot")).observe(get("experiment-plot"));
  render();
})();
