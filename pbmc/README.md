# PBMC Cluster Lab

An interactive BIOSTAT 615 project for data science students new to single-cell RNA sequencing. The lab uses a small set of real human peripheral blood mononuclear cells (PBMCs) to demonstrate how expression measurements become clusters and how marker evidence supports a cell-type annotation.

## Run locally

From the repository root, start a static file server. No JavaScript dependencies or build step are required.

```sh
python -m http.server 8000 --directory dist
```

Open [the clustering lab](http://localhost:8000/) or [the annotation lab](http://localhost:8000/two-round.html).

## Try the experiments

1. **Clustering lab.** Run PCA, set the number of clusters and the weights on the principal components, then run Gaussian-mixture EM. Inspect a cell to see how its cluster probabilities change. The cluster numbers identify groups fitted by the model.
2. **Annotation lab.** Click **Initialize model** and alternate **Run E-step** and **Run M-step**, or use **Auto-step EM**. Click **Next mixed cell** to step through ambiguous cells. Compare the raw marker counts and probabilities in the selected cell's detail panel. The two halves of each split dot identify its two leading candidates. Change the anchor threshold to see which cells switch between solid and split dots. Hollow dots show cells for which the background model ranks first.
3. **Optional review.** Use genes outside the marker panel to see whether they support the first candidate for the selected ambiguous cell. The main map continues to show the first-round annotations.

The annotation lab uses a simplified, fixed-dispersion negative-binomial model inspired by [CellAssign](https://pmc.ncbi.nlm.nih.gov/articles/PMC7485597/). It fits T-cell, B-cell, monocyte and background components using seven marker genes and each cell's library size. PCA supplies the display coordinates. The annotation fit uses the marker counts and does not load the provider's cell labels. The probabilities and threshold are for comparing evidence within this teaching model.

## Data and attribution

`dist/pbmc.js` is derived from [5k Human PBMCs (Donor 1) with Automated Cell Annotation](https://www.10xgenomics.com/datasets/5k_Human_Donor1_PBMC_3p_gem-x), published by **10x Genomics** under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). The teaching copy samples 256 cells with a fixed seed, retains raw counts for a 160-gene panel and the full-library total per cell, and omits original barcodes and provider cell-type labels. CC BY 4.0 applies to the included dataset derivative.

The teaching copy is bundled. To regenerate it, download the **Filtered feature barcode matrix (HDF5)** from the dataset page, install the R packages `Seurat`, `Matrix`, and `jsonlite`, and run this from the repository root:

```sh
Rscript prepare-data.R path/to/filtered_feature_bc_matrix.h5
```

## Tests

With Node.js installed, run:

```sh
node --test cluster.test.js marker.test.js
```

The tests verify PCA against an independent R calculation, EM likelihood and probability updates, the effects of clustering parameters and annotation thresholds, and the feature inputs to the optional review.

## AI-assisted development

Codex was used to implement and test changes after course feedback. The project evolved from a result display into a browser-based experiment with adjustable parameters, step-by-step EM, per-cell calculations, and a visual distinction between uncertain and unassigned cells. The final navigation update was tested by stepping through every ambiguous cell at the default threshold.
