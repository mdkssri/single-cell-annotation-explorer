"use strict";

(() => {
  const { mark, size, defaultThreshold, convergence, minimumAnchors, prepare, initialize, expectation, maximize, assignment, review } = MarkerModel;
  const get = id => document.getElementById(id);
  const palette = ["#135cc9", "#a51a80", "#b75000", "#53605d"];
  const data = prepare(PBMC);
  const fmt = (value, places = 3) => Number.isFinite(value) ? value.toFixed(places) : "−∞";
  const fmtPositive = value => value < 1e-5 ? value.toExponential(3) : fmt(value, 5);
  const pct = value => (value * 100).toFixed(1) + "%";
  const idOf = i => "PBMC-" + String(i + 1).padStart(3, "0");
  const threshold = () => get("threshold").valueAsNumber / 100;
  const note = message => { get("status").textContent = message; };
  let points = null;
  let model = null;
  let posterior = null;
  let oldModel = null;
  let next = "E";
  let history = [];
  let iterations = 0;
  let selected = 0;
  let timer = null;
  let converged = false;
  let lastDelta = null;
  let positions = [];
  let calls = [];
  let mixed = [];
  let check = null;
  let activeType = null;

  get("cell").innerHTML = PBMC.cells.map((_, i) => `<option value="${i}">${idOf(i)}</option>`).join("");
  get("threshold").value = defaultThreshold * 100;
  get("legend").innerHTML = mark.map((item, k) =>
    `<span><i style="--cluster:${palette[k]}${k === 3 ? ";background:white;border:1px solid " + palette[k] : ""}"></i>${item.name}</span>`).join("");

  function plot() {
    const canvas = get("marker-plot");
    const { width, height } = canvas.getBoundingClientRect();
    if (!width || !height) return;
    const ratio = Math.min(window.devicePixelRatio || 1, 2);
    canvas.width = Math.round(width * ratio);
    canvas.height = Math.round(height * ratio);
    const ctx = canvas.getContext("2d");
    ctx.scale(ratio, ratio);
    positions = [];
    if (!points) return;
    const bounds = [0, 1].map(j => {
      const values = points.scores.map(point => point[j]);
      const min = Math.min(...values), max = Math.max(...values);
      return { mid: (min + max) / 2, span: 1.12 * (max - min) };
    });
    const unit = Math.min((width - 75) / bounds[0].span, (height - 56) / bounds[1].span);
    positions = points.scores.map(([x, y]) => [(width + 17) / 2 + (x - bounds[0].mid) * unit, (height - 14) / 2 - (y - bounds[1].mid) * unit]);
    ctx.strokeStyle = "#dae5dc";
    for (let column = 39; column < width - 15; column += 38) {
      ctx.beginPath(); ctx.moveTo(column, 14); ctx.lineTo(column, height - 26); ctx.stroke();
    }
    for (let row = 18; row < height - 23; row += 32) {
      ctx.beginPath(); ctx.moveTo(28, row); ctx.lineTo(width - 14, row); ctx.stroke();
    }
    ctx.fillStyle = "#4d6660"; ctx.font = "12px system-ui"; ctx.textAlign = "center";
    ctx.fillText("PC1", width / 2, height - 7);
    ctx.save(); ctx.translate(15, height / 2); ctx.rotate(-Math.PI / 2); ctx.fillText("PC2", 0, 0); ctx.restore();
    function point(i) {
      const [x, y] = positions[i];
      const result = calls[i];
      const radius = result?.kind === "ambiguous" ? (i === selected ? 9 : 6.2) : (i === selected ? 6.2 : 3.9);
      ctx.beginPath(); ctx.arc(x, y, radius, 0, Math.PI * 2);
      if (result?.kind === "unassigned") {
        ctx.fillStyle = "#fbfdfb"; ctx.fill();
        ctx.strokeStyle = palette[3]; ctx.lineWidth = 1.3; ctx.stroke();
      } else if (result?.kind === "ambiguous") {
        ctx.strokeStyle = "#fff"; ctx.lineWidth = 2.5; ctx.stroke();
        for (const [half, candidate] of [result.top, result.next].entries()) {
          const start = (half - .5) * Math.PI;
          ctx.beginPath(); ctx.arc(x, y, radius, start, start + Math.PI);
          ctx.closePath(); ctx.fillStyle = palette[candidate.k]; ctx.fill();
        }
        ctx.beginPath(); ctx.moveTo(x, y - radius); ctx.lineTo(x, y + radius);
        ctx.strokeStyle = "#fff"; ctx.lineWidth = 1.3; ctx.stroke();
      } else {
        ctx.fillStyle = result ? palette[result.top.k] : "#99aaa1";
        ctx.fill();
      }
      if (i === selected) {
        ctx.strokeStyle = "#163129"; ctx.lineWidth = 2;
        ctx.beginPath(); ctx.arc(x, y, radius + 3, 0, Math.PI * 2); ctx.stroke();
      }
    }
    positions.forEach((_, i) => { if (i !== selected && calls[i]?.kind !== "ambiguous") point(i); });
    mixed.forEach(i => { if (i !== selected) point(i); });
    point(selected);
  }

  function renderCell() {
    const cell = PBMC.cells[selected];
    const current = posterior?.rows[selected];
    const call = calls[selected];
    get("detail-title").textContent = idOf(selected);
    get("detail-status").textContent = !call ? "Awaiting the E-step" : call.kind === "unassigned" ? "Unassigned · background model leads"
      : call.kind === "anchor" ? `Anchor · ${mark[call.top.k].name}` : `Ambiguous · ${mark[call.top.k].name} / ${mark[call.next.k].name}`;
    get("marker-counts").innerHTML = mark.slice(0, 3).map(item =>
      `<div><strong>${item.name}</strong>${item.genes.map(gene => `<span>${gene} · <b>${cell[PBMC.genes.indexOf(gene)]}</b></span>`).join("")}</div>`).join("");
    get("probabilities").innerHTML = current ? current.probabilities.map((p, k) =>
      `<div class="mark-bar"><span>${mark[k].name}</span><div><i style="--amount:${p * 100}%;--color:${palette[k]}"></i></div><span>${pct(p)}</span></div>`).join("")
      : '<p class="muted">Run the E-step to calculate this cell’s probabilities.</p>';
    if (!current) {
      get("equations").innerHTML = '<p class="muted">The E-step uses the observed counts, current marker rates and type proportions.</p>';
      return;
    }
    const k = activeType ?? call.top.k;
    const geneIndex = k === 3 ? 0 : data.genes.findIndex(gene => gene.k === k);
    const gene = data.genes[geneIndex];
    const result = current.values[k];
    const term = result.terms[geneIndex];
    const oldRates = oldModel.rates[geneIndex];
    const newRates = model.rates[geneIndex];
    const rate = term.mean / data.depth[selected];
    const m = next === "E" ? `<p>New mixing proportion for ${mark[k].name}: <code>${fmt(model.priors[k] * data.count)} ÷ ${data.count} = ${fmt(model.priors[k], 5)}</code></p>
      <p>${gene.name} marker rate: <code>${fmt(oldRates.high, 5)} → ${fmt(newRates.high, 5)}</code></p>
      <p>${gene.name} background rate: <code>${fmt(oldRates.low, 5)} → ${fmt(newRates.low, 5)}</code></p>` : '<p class="muted">Run the M-step to see these rates and mixing proportions update.</p>';
    get("equations").innerHTML = `<label for="trace-type">Type to trace <select id="trace-type">${mark.map((item, i) =>
      `<option value="${i}"${i === k ? " selected" : ""}>${item.name}</option>`).join("")}</select></label>
      <div class="marker-equations"><section><h4>E-step · ${gene.name} example</h4>
      <p>Expected count: <code>${fmt(data.depth[selected])} × ${fmtPositive(rate)} ≈ ${fmtPositive(term.mean)}</code></p>
      <p>Count log-likelihood term: <code>${term.count} ln(${fmtPositive(term.mean)}) − (${term.count} + ${size}) ln(${size} + ${fmtPositive(term.mean)}) ≈ ${fmt(term.score, 5)}</code></p>
      <p>Log score across markers: <code>${fmt(result.terms.reduce((total, item) => total + item.score, 0), 5)} + (${fmt(result.logPrior, 5)}) ≈ ${fmt(result.logScore, 5)}</code></p>
      <p>Responsibility: <code>exp(${fmt(result.logScore, 5)} − (${fmt(current.total, 5)})) ≈ ${pct(current.probabilities[k])}</code></p></section>
      <section><h4>M-step · estimated parameters</h4>${m}<p class="muted">Rates are fit using all 256 cells, weighted by their responsibilities.</p></section></div>`;
  }

  function renderReview() {
    const panel = get("review-result");
    if (!check) {
      panel.textContent = "Run the first-round model, then review the selected ambiguous cell.";
    } else if (!check.rows) {
      panel.textContent = `Needs at least ${minimumAnchors} anchors per type at this threshold (T: ${check.sizes[0]}, B: ${check.sizes[1]}, monocyte: ${check.sizes[2]}).`;
    } else {
      const result = check.rows[selected];
      panel.textContent = !result ? `${idOf(selected)} is outside the ambiguous-cell review set.`
        : `${idOf(selected)} · ${mark[result.winner].name} ${pct(result.probability)} in the review model. ${result.corroborates ? "The other genes support the first-round candidate." : "Evidence remains inconclusive for the first-round candidate."}`;
    }
  }

  function render() {
    const running = timer !== null;
    get("e-step").disabled = !model || next !== "E" || running;
    get("m-step").disabled = !model || next !== "M" || running;
    get("auto").disabled = !model || (!running && converged);
    get("auto").textContent = running ? "Pause EM" : "Auto-step EM";
    get("reset").disabled = !model;
    get("initialize").disabled = Boolean(model);
    get("plot-empty").hidden = Boolean(points);
    get("phase").textContent = !model ? "Initialize to start from the raw count panel."
      : next === "M" ? "The latest E-step supplies the map colors."
      : iterations ? "M-step updated the rates. Colors show the previous E-step until you run E again."
        : "The model is initialized. Run E to reveal cell-type probabilities.";
    get("likelihood").textContent = history.length ? fmt(history.at(-1), 4) : "—";
    get("progress").textContent = model ? `${iterations} completed M-steps${lastDelta === null ? "" : " · Δ " + lastDelta.toExponential(2)}${converged ? " · converged" : ""}` : "";
    get("threshold-value").textContent = Math.round(threshold() * 100) + "%";
    const categories = { anchor: 0, ambiguous: 0, unassigned: 0 };
    calls = posterior?.rows.map(entry => assignment(entry, threshold())) ?? [];
    mixed = [];
    calls.forEach((call, i) => {
      categories[call.kind] += 1;
      if (call.kind === "ambiguous") mixed.push(i);
    });
    get("totals").innerHTML = [["anchor", "Anchors"], ["ambiguous", "Mixed color"], ["unassigned", "Unassigned"]].map(([key, label]) =>
      `<div><b>${posterior ? categories[key] : "—"}</b>${label}</div>`).join("");
    const mixedIndex = mixed.indexOf(selected);
    get("previous-mixed").disabled = !mixed.length;
    get("find").disabled = !mixed.length;
    get("mixed-position").textContent = mixedIndex >= 0 ? `${mixedIndex + 1} of ${mixed.length} · ${idOf(selected)}`
      : `${mixed.length} mixed-color ${mixed.length === 1 ? "cell" : "cells"}`;
    get("review").disabled = !posterior || running;
    plot(); renderCell(); renderReview();
  }

  function pause() {
    clearInterval(timer);
    timer = null;
  }

  function step() {
    if (next === "E") {
      oldModel = model;
      posterior = expectation(data, model);
      next = "M";
      check = null;
    } else {
      model = maximize(data, posterior);
      const likelihood = expectation(data, model).logLikelihood;
      lastDelta = likelihood - history.at(-1);
      history.push(likelihood);
      iterations += 1;
      next = "E";
      converged = Math.abs(lastDelta) < convergence;
    }
    render();
  }

  function select(i) {
    selected = i;
    activeType = null;
    get("cell").value = i;
    render();
  }

  get("initialize").addEventListener("click", () => {
    get("initialize").disabled = true;
    note("Computing the display PCA and initializing the marker model…");
    requestAnimationFrame(() => {
      try {
        points ??= ClusterModel.pca(PBMC, 4);
        model = initialize(data);
        history = [expectation(data, model).logLikelihood];
        get("axes").textContent = `PC1 ${pct(points.fractions[0])} · PC2 ${pct(points.fractions[1])} · display coordinates`;
        note("Ready. Run an E-step to see the first probabilities.");
      } catch (error) {
        note(error.message);
      }
      render();
    });
  });
  get("e-step").addEventListener("click", () => { step(); note("E-step refreshed the probabilities and point colors."); });
  get("m-step").addEventListener("click", () => { step(); note(converged ? "Convergence reached. Run E to view the final probabilities." : "M-step updated rates and proportions."); });
  get("auto").addEventListener("click", () => {
    if (timer !== null) { pause(); render(); note("Paused. Continue with a manual step or resume playback."); return; }
    const limit = iterations + 100;
    timer = setInterval(() => {
      step();
      if (!converged && iterations < limit) return;
      pause();
      if (next === "E") step();
      note(converged ? `Converged at round ${iterations}. Final E-step shown.` : "Paused after 100 rounds. Continue fitting with Auto-step EM.");
    }, 240);
    render();
    note("Alternating E- and M-steps.");
  });
  get("reset").addEventListener("click", () => {
    pause();
    model = posterior = oldModel = check = null;
    history = [];
    iterations = 0;
    converged = false;
    lastDelta = null;
    next = "E";
    activeType = null;
    get("threshold").value = defaultThreshold * 100;
    note("Model reset. Initialize to begin a new fit.");
    select(0);
  });
  get("threshold").addEventListener("input", () => { check = null; render(); });
  get("cell").addEventListener("change", event => select(Number(event.target.value)));
  function browseMixed(direction) {
    if (!mixed.length) return;
    const index = mixed.indexOf(selected);
    const start = direction === 1 ? 0 : mixed.length - 1;
    const nextIndex = index < 0 ? start : (index + direction + mixed.length) % mixed.length;
    select(mixed[nextIndex]);
  }
  get("previous-mixed").addEventListener("click", () => browseMixed(-1));
  get("find").addEventListener("click", () => browseMixed(1));
  get("marker-plot").addEventListener("click", event => {
    const box = get("marker-plot").getBoundingClientRect();
    let nearest = -1, best = 13 ** 2;
    positions.forEach(([x, y], i) => {
      const d = (event.clientX - box.left - x) ** 2 + (event.clientY - box.top - y) ** 2;
      if (d < best) { nearest = i; best = d; }
    });
    if (nearest >= 0) select(nearest);
  });
  get("equations").addEventListener("change", event => {
    if (event.target.id !== "trace-type") return;
    activeType = Number(event.target.value);
    renderCell();
    get("trace-type").focus({ preventScroll: true });
  });
  get("review").addEventListener("click", () => {
    check = review(PBMC, posterior, threshold(), ClusterModel.pca, ClusterModel.evaluate);
    renderReview();
  });
  new ResizeObserver(plot).observe(get("marker-plot"));
  render();
})();
