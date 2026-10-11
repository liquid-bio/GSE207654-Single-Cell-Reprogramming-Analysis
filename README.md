# GSE207654-Single-Cell-Reprogramming-Analysis
## 📜 Project Overview  
This repository contains a complete computational reproduction of the single-cell RNA sequencing (scRNA-seq) analysis from the following study:  

Qin, J., Zhang, J., Jiang, J., Zhang, B., Li, J., Lin, X., Wang, S., Zhu, M., Fan, Z., Lv, Y., He, L., Chen, L., Yue, W., Li, Y., & Pei, X. (2022). Direct chemical reprogramming of human cord blood erythroblasts to induced megakaryocytes that produce platelets. Cell Stem Cell, 29(8), 1229–1245.e7. https://doi.org/10.1016/j.stem.2022.07.004

Objective: To replicate the core computational findings of the original paper by independently implementing the full scRNA-seq analysis workflow on the publicly available dataset (GSE207654), which examines the reprogramming of erythroblasts into functional megakaryocytes using a small-molecule cocktail.

Context: This project was completed as an advanced coursework assignment, demonstrating proficiency in single-cell genomics analysis and reproducible computational research. The coursework report (English translation) is in [report/report.md](report/report.md).

## 🔬 Key Analysis Replicated  
The analysis faithfully follows the methodological pipeline described in the original publication, covering:  

Raw Data Processing: Loading and merging 10x Genomics data from four reprogramming time points (EB-D0, iMK-D3, iMK-D5, iMK-D7).  

Rigorous Quality Control: Filtering cells based on gene counts, mitochondrial read percentage, and other standard metrics.  

Data Normalization & Integration: Applying log normalization, scaling, and batch correction using BBKNN to harmonize data across time points.  

Dimensionality Reduction & Clustering: Performing PCA, UMAP visualization, and Leiden clustering to identify cell populations.  

Cell Type Annotation: Annotating clusters based on canonical marker genes (e.g., GYPA for erythroblasts, GP1BA for megakaryocytes) as defined in the paper's Figures 5e-g.  

Biological Interpretation: Conducting differential expression and Gene Ontology (GO) enrichment analysis to validate the biological states of identified clusters.  


## 🚀 Quick Start
1. Clone the Repository  
bash  
git clone https://github.com/liquid-bio/GSE207654-Single-Cell-Reprogramming-Analysis.git  
cd GSE207654-Single-Cell-Reprogramming-Analysis   
2. Set Up the Computational Environment  
We recommend using Conda to create an isolated environment with all dependencies:  
bash  
conda env create -f environment.yaml  
conda activate scRNA  
3. Download the Data  
The raw sequencing data must be downloaded from the Gene Expression Omnibus (GEO):  
Accession: GSE207654  
Follow the instructions in data/GSE207654/tips.txt and put the required MTX and TSV files in data/GSE207654/.  
4. Execute the Analysis  
Open and run the Jupyter notebook to reproduce the entire analysis:  
bash  
jupyter notebook notebooks/code.ipynb  
