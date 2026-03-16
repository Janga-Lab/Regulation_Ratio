# Regulation_Ratio
**A singular multi-omic measurement of gene regulatory mechanisms**

Alexander Krohannon, Mansi Srivastava, Neel Sangani, Sarath Chandra Janga  
*Department of Biomedical Engineering and Informatics, Luddy School of Informatics, Computing, and Engineering, Indiana University Indianapolis*

Corresponding Authors: [amkrug@iu.edu](mailto:amkrug@iu.edu), [scjanga@iu.edu](mailto:scjanga@iu.edu)

---

## Overview

This repository contains the code and reference data used to compute and analyze the **Regulation Ratio (RR)** — a normalized, gene-level metric that quantifies the relative contributions of transcriptional and post-transcriptional regulatory mechanisms by integrating ATAC-seq, RNA-seq, and Protein Occupancy Profiling sequencing (POP-seq) data.

The framework was applied across three human cell lines — **HEK293**, **HepG2**, and **K562** — revealing three distinct regulatory classes: predominantly transcriptionally regulated, predominantly post-transcriptionally regulated, and balanced genes.

Raw and processed sequencing data are accessible at the Sequence Read Archive (SRA) under ID: **16057836**.

---

## Background

Gene expression in eukaryotic cells is governed by regulatory mechanisms operating at multiple levels. This project focuses on two major levels:

- **Transcriptional regulation** — quantified using ATAC-seq peaks mapping chromatin accessibility within an expanded gene boundary (5,000 bp upstream and 500 bp downstream), normalized by gene length in kilobases.
- **Post-transcriptional regulation** — quantified using POP-seq peaks mapping protein-RNA interaction sites within gene boundaries, normalized by gene length in kilobases.

The **Regulation Ratio** is defined as:

```
RR = (# Significant POP-seq peaks) / (# Significant ATAC-seq peaks)
```

Genes with RR > 1 are classified as **predominantly post-transcriptionally regulated**, RR < 1 as **predominantly transcriptionally regulated**, and RR = 1 as **balanced**.

---

## Key Findings

- Post-transcriptional regulation exhibited a **>4-fold higher** median density compared to transcriptional regulation in HEK293 cells, while HepG2 and K562 (both cancer-derived) showed the opposite, with ~2.5-fold higher transcriptional regulation.
- A **strong correlation** between transcriptional and post-transcriptional regulation was observed across all three cell lines (Pearson R² ≈ 0.60), suggesting coordinated rather than independent operation.
- **Transcriptionally regulated genes** tend to be significantly longer (median ~54,000 bp), while **post-transcriptionally regulated genes** have more isoforms and higher transcript abundance (median FPKM: 310 vs. 34 vs. 6.75 for post-transcriptional, transcriptional, and balanced, respectively).
- **55.8% of genes** maintained identical regulatory classification across all three cell lines, with enriched functional pathways showing consistent directionality — suggesting regulatory strategy is a fundamental property of gene function, not just cellular context.
- Proliferation-associated pathways (Myc Targets, E2F Targets, G2M Checkpoint, mTORC1 Signaling) showed consistent post-transcriptional enrichment across all three cell lines, while developmental and tissue-identity pathways (Hedgehog Signaling, Epithelial Mesenchymal Transition, Notch Signaling) were consistently transcriptionally enriched.

---

## Repository Structure

```
Regulation_Ratio/
├── Regulation_Count.py       # Computes the Regulation Ratio profile for a cell line
├── Profile_Analyzer.py       # Analyzes and visualizes a computed profile
├── Gene_Information.txt      # Ensembl gene annotations (GRCh38)
├── Gene_Type_Data.txt        # Ensembl transcript type annotations
├── Protein_Information.txt   # UniProt protein annotations for DepMap integration
├── gene_regulation.yml       # Conda environment specification
├── HEK293/                   # Precomputed profiles and results for HEK293
├── HepG2/                    # Precomputed profiles and results for HepG2
└── K562/                     # Precomputed profiles and results for K562
```

---

## Installation

It is recommended to use the provided Conda environment:

```bash
conda env create -f gene_regulation.yml
conda activate gene_regulation
```

The environment includes all required dependencies: `numpy`, `pandas`, `scipy`, `matplotlib`, `seaborn`, `statannot`, `cliffs_delta`, and more.

---

## Usage

### Step 1 — Generate a Regulation Ratio Profile (`Regulation_Count.py`)

This script takes ATAC-seq and POP-seq peak files for a given cell line and computes per-gene regulatory metrics, including the Regulation Ratio and regulatory group classification. Optionally, isoform counts, RNA-seq expression, and proteomics data can be incorporated.

**Required inputs:**

| Argument | Description |
|---|---|
| `cell_type` | Cell line name (e.g., `HEK293`) |
| `atac_file` | Path to ATAC-seq peak file (tab-separated BED-like format) |
| `pop_file` | Path to POP-seq peak file (tab-separated BED-like format) |

**Optional inputs:**

| Argument | Description |
|---|---|
| `--iso_file` | Path to StringTie GTF file for isoform counting |
| `--rna_file` | Path to Ballgown CSV for transcript abundance |
| `--prot_file` | Path to DepMap proteomics CSV (default: `DepMap_Proteomics.csv`) |
| `--upstream` | Bases upstream of gene boundary for ATAC peak detection (default: 5000) |
| `--downstream` | Bases downstream of gene boundary for ATAC peak detection (default: 500) |

**Example:**

```bash
python Regulation_Count.py HEK293 HEK293_ATAC_peaks.bed HEK293_POP_peaks.bed \
    --iso_file HEK293_isoforms.gtf \
    --rna_file HEK293_rna.csv \
    --prot_file DepMap_Proteomics.csv
```

**Output:** A CSV file named `<cell_type>_Piranha_Profile.csv` (or `<cell_type>_Piranha_Profile_<upstream>.csv` if `--upstream` is specified), containing per-gene regulatory metrics and group classifications.

**Peak file format:** Tab-separated with at least three columns — chromosome, peak start, and peak end (no `chr` prefix in the chromosome column). This matches standard output from MACS3 (ATAC-seq) and Piranha (POP-seq).

---

### Step 2 — Analyze and Visualize a Profile (`Profile_Analyzer.py`)

This script takes a profile CSV generated by `Regulation_Count.py` and produces a suite of statistical comparisons and publication-quality figures.

**Required inputs:**

| Argument | Description |
|---|---|
| `profile` | Path to the profile CSV file |
| `output` | A label used in output file names and plot titles |

**Example:**

```bash
python Profile_Analyzer.py HEK293_Piranha_Profile.csv HEK293
```

**Outputs:**

- `<output>_Regulation_Density.png` — Box plot comparing transcriptional vs. post-transcriptional regulatory densities across all genes
- `<output>_Regulation_Ratio.png` — Scatter plot of regulatory densities colored by regulatory group (Regulation Ratio correlation)
- `<output>_Gene_Length.png` — Violin plot of gene lengths across regulatory groups
- `<output>_Isoform_Count.png` — Violin plot of isoform counts across regulatory groups
- `<output>_Transcript_Abundance.png` — Violin plot of transcript abundance (FPKM) across regulatory groups
- `<output>_Protein_Expression.png` — Violin plot of protein expression levels across regulatory groups
- `<output>_Categorical.png` — Bar chart of transcript type distribution and enrichment across regulatory groups
- `<output> Transcript Type Enrichment.csv` — Fisher's exact test results for transcript type enrichment

> **Note:** `Profile_Analyzer.py` expects `Gene_Type_Data.txt` to be present in the working directory for transcript type enrichment analysis.

---

## Input Data Preparation

The peak files used as input for `Regulation_Count.py` should be generated using the following tools, consistent with the methods in the associated publication:

- **ATAC-seq peaks:** Called with [MACS3](https://github.com/macs3-project/MACS) using a shift of 50 bp and extension of 100 bp, with FDR correction (p < 0.05). Biological replicates should be pooled prior to peak calling.
- **POP-seq peaks:** Called with [Piranha](http://smithlabresearch.org/software/piranha/) (v1.2.1) using a bin size of 20 bp, no normalization, and logarithmic covariates, with FDR correction (p < 0.05). Biological replicates should be pooled prior to peak calling.
- All data should be aligned to **GRCh38/hg38** using [HISAT2](http://daehwankimlab.github.io/hisat2/).

---

## Citation

If you use this code or the Regulation Ratio framework in your work, please cite:

> Krohannon A, Srivastava M, Sangani N, Janga SC. *Regulation Ratio: A singular multi-omic measurement of gene regulatory mechanisms.* (2025).

---

## Funding

This work was supported by the National Institute of General Medical Sciences under Award Number **R01GM123314**, a Lilly Research Award Program (LRAP) grant, and the IUI Luddy School of Informatics, Computing and Engineering.

---

## License

This project is licensed under the GPL-3.0 License. See [LICENSE](LICENSE) for details.
