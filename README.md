# Single-Cell Annotation Explorer

**[Live demo](https://mdkssri.github.io/single-cell-annotation-explorer/)**

Static SingleR teaching app.

**Local preview only:** From the repository root, run `python -m http.server 8771 --bind 127.0.0.1`, then open the address printed by the terminal and select `docs/`.

GitHub Pages: upload this folder contents, not the folder itself, then select main and /docs in Settings > Pages. Do not upload the original project folder or raw data.

## Results
Original gene x cell dimensions: 38606 x 5710; 38606 x 5989
Annotated barcodes matched: 5710 of 5710
Reference cells: 912
Query cells: 900
Shared genes: 38606
Reference labels: B cell 180; leukocyte 47; lymphocyte 104; macrophage 173; monocyte 180; mononuclear cell 48; T cell 180
Pruned: 48
Minimum delta.next: 1.11022302462516e-16
Highest initial label differs from final: 35
SingleR version: 2.12.0
Two computations and export validation: PASS

## Reproducibility
Use R 4.5.2, Bioconductor 3.22, SingleR 2.12.0, Seurat, rhdf5, SingleCellExperiment, scuttle, scrapper and matrixStats. See 10x-session.txt.
Download these inputs into a local input folder (not in GitHub):
donor1-counts.h5 | https://cf.10xgenomics.com/samples/cell-exp/9.0.0/5k_Human_Donor1_PBMC_3p_gem-x_5k_Human_Donor1_PBMC_3p_gem-x/5k_Human_Donor1_PBMC_3p_gem-x_5k_Human_Donor1_PBMC_3p_gem-x_count_sample_filtered_feature_bc_matrix.h5 | 9945e1b94d62375d4b845ac8cc2bbae9
donor2-counts.h5 | https://cf.10xgenomics.com/samples/cell-exp/9.0.0/5k_Human_Donor2_PBMC_3p_gem-x_5k_Human_Donor2_PBMC_3p_gem-x/5k_Human_Donor2_PBMC_3p_gem-x_5k_Human_Donor2_PBMC_3p_gem-x_count_sample_filtered_feature_bc_matrix.h5 | 009734dc3fa26732c1651041c698dfb0
donor1-cell-types.tar.gz | https://cf.10xgenomics.com/samples/cell-exp/9.0.0/5k_Human_Donor1_PBMC_3p_gem-x_5k_Human_Donor1_PBMC_3p_gem-x/5k_Human_Donor1_PBMC_3p_gem-x_5k_Human_Donor1_PBMC_3p_gem-x_count_cell_types.tar.gz | 62bab776d067007383a0aceeebbb4c5a
Run from this folder: Rscript scripts/10_prepare_10x.R . INPUT_DIR NEW_OUTPUT_DIR (the last directory must not exist).
The script runs SingleR twice and checks the six TSVs against R objects. For publication remove patient, sample and provided_annotation from cells.tsv, preserving the project cell IDs; retain the attribution file.
Fixed seed 6152026; reference types with fewer than 30 provided labels excluded, at most 180 sampled per retained type; 900 query cells sampled without using predictions.

## Data credits
Data: 10x Genomics (2024), 5k Human PBMCs (Donor 1) with Automated Cell Annotation and 5k Human PBMCs (Donor 2), Cell Ranger 9.0.0.
https://www.10xgenomics.com/datasets/5k_Human_Donor1_PBMC_3p_gem-x
https://www.10xgenomics.com/datasets/5k_Human_Donor2_PBMC_3p_gem-x
License: CC BY 4.0 https://creativecommons.org/licenses/by/4.0/
Changes: donor-separated reference and query; deterministic subsetting; SingleR scoring, query PCA, and reference-expression summaries. Public files replace original barcodes with project IDs and omit donor and sample columns.
Reference labels are 10x coarse automated labels, not ground truth. Donor 2 has no query truth labels. This teaching example does not establish generalization.

## Code license
Project code: MIT (see LICENSE). The 10x source datasets are separately attributed under CC BY 4.0; see docs/data/processed/ATTRIBUTION.txt.
