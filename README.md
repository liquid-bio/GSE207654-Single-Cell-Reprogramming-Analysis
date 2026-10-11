# GSE207654-Single-Cell-Reprogramming-Analysis

## 📜 Project Overview

This repository contains a computational reproduction of the single-cell RNA sequencing (scRNA-seq) analysis from the following study:

Qin, J., Zhang, J., Jiang, J., Zhang, B., Li, J., Lin, X., Wang, S., Zhu, M., Fan, Z., Lv, Y., He, L., Chen, L., Yue, W., Li, Y., & Pei, X. (2022). Direct chemical reprogramming of human cord blood erythroblasts to induced megakaryocytes that produce platelets. *Cell Stem Cell*, 29(8), 1229–1245.e7. https://doi.org/10.1016/j.stem.2022.07.004

**Objective**: To replicate the core computational findings of the original paper by independently implementing the full scRNA-seq analysis workflow on the publicly available dataset ([GSE207654](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207654)), which examines the reprogramming of erythroblasts into functional megakaryocytes using a small-molecule cocktail.

**Context**: This project was completed as an advanced coursework assignment, demonstrating proficiency in single-cell genomics analysis and reproducible computational research. The coursework report (English translation) is in [report/report.md](report/report.md).

## 🔬 Key Analysis Replicated

The analysis follows the methodological pipeline described in the original publication, covering:

- **Raw data processing**: loading and merging 10x Genomics data from four reprogramming time points (EB-D0, iMK-D3, iMK-D5, iMK-D7).
- **Quality control**: filtering cells by library size (< 2,500 UMIs) and mitochondrial fraction (> 15%), the thresholds of the paper; Scrublet doublet scores.
- **Normalization and integration**: log normalization, highly variable genes, cell-cycle score regression, scaling, and batch correction with BBKNN across time points.
- **Dimensionality reduction and clustering**: PCA, UMAP visualization and Leiden clustering.
- **Cell type annotation**: annotating clusters as EB, EB-like, iPEM, iMK-1 and iMK-2 from canonical marker genes (e.g., GYPA for erythroblasts, GP1BA for megakaryocytes) shown in the paper's Figures 5D–5I.
- **Biological interpretation**: differential expression between iMK-1 and iMK-2 and enrichment analysis with Enrichr (gseapy).
- **Doublet analysis**: an additional comparison of clustering with and without Scrublet doublet removal.

## 📁 Repository Structure

```
├── notebooks/code.ipynb   # full analysis
├── data/GSE207654/        # put the GEO files here (see tips.txt)
├── results/               # figures from the notebook
├── report/                # coursework report (English) and its figures
└── environment.yaml       # conda environment (exported on Windows)
```

## 📊 Results

| File | Content | Paper figure |
|---|---|---|
| `results/scatterPlot.pdf` | total counts vs detected genes, colored by mitochondrial fraction | – |
| `results/umap_non.pdf` | UMAP by time point, before batch correction | 5A |
| `results/umap_scanorama.pdf` | UMAP by time point, after batch correction (BBKNN) | 5B |
| `results/leiden_marker.pdf` | marker dot plot per Leiden cluster | 5E–5G |
| `results/celltype.pdf` | UMAP by annotated cell type | 5B |
| `results/celltype_marker.pdf` | marker dot plot per cell type | 5D |
| `results/celltype_sample_proportion.pdf` | sample composition of each cell type | 5C |

`results/celltype_marker2.pdf` (stacked violin plot) was saved blank; the same plot is Figure 11 of the report. The report shows all figures, including QC, cell-cycle correction, enrichment and the doublet analysis.

## 🚀 Quick Start

1. Clone the repository:

   ```bash
   git clone https://github.com/liquid-bio/GSE207654-Single-Cell-Reprogramming-Analysis.git
   cd GSE207654-Single-Cell-Reprogramming-Analysis
   ```

2. Set up the environment. Conda is recommended:

   ```bash
   conda env create -f environment.yaml
   conda activate scRNA
   ```

   The file was exported on Windows and contains Windows-only packages. On Linux or macOS, create the environment from the main packages instead: Python 3.9, scanpy 1.10.3, anndata 0.10.9, bbknn, leidenalg, scikit-image, gseapy 1.1.11, seaborn and jupyter.

3. Download the data from GEO (accession GSE207654). Put the 12 per-sample files (`GSM63044xx_*.barcodes.tsv.gz`, `*.genes.tsv.gz`, `*.matrix.mtx.gz`) in `data/GSE207654/`, as described in [data/GSE207654/tips.txt](data/GSE207654/tips.txt).

4. Run the notebook from the `notebooks/` folder:

   ```bash
   jupyter notebook notebooks/code.ipynb
   ```

   The notebook changes to `data/`, copies each sample into `data/outputs/`, and writes the merged object and figures to `data/1-DataPre/`.

## 📝 Notes

- **Corrupted GEO file**: the current GEO matrix of GSM6304416 (iMK_D7) is corrupted, both as a per-sample file and in `GSE207654_RAW.tar`. After 35,944,352 of the 35,951,943 declared entries it continues with an unrelated table, so `scipy.io.mmread` fails on it. Reading only the valid entries gives the same filtered sample as in this analysis (10,603 cells, 18,465 genes).
- **iMK-2 and doublets**: the paper does not mention doublet removal. The iMK_2 cluster in this analysis disappears after Scrublet doublet removal and has a high mean doublet score (0.32 vs 0.065 overall), consistent with doublets rather than the paper's CD53+ iMK-2 (report, Section 8).
- **Enrichment labels**: the notebook requests the GO and KEGG libraries in one `gp.enrichr` call. With gseapy 1.1.11 this returns the KEGG results under both labels, so the plots titled "GO" show KEGG pathways.
- **Batch correction**: the paper uses mutual nearest neighbors within Monocle 3; this analysis uses BBKNN, so cluster numbers and boundaries differ from the paper.
