# TReATAD + BV2 + BioHyV – Reproducibility Repository

This repository contains code for the analyses and figures in:

> **[Paper Title]**  
> Authors: [List of Authors]  
> Journal / Preprint / Year: [Journal name, Year or bioRxiv link]  
> DOI: [Insert DOI once minted]

This is the minimal, paper-focused version of the project, containing only the scripts and data pointers necessary to reproduce the main results.

---

## 📄 Structure of this repository

treatad_bv2_biohyv_paper/
├── data/ # (Optional) small example or subset data
├── scripts/ # Main analysis scripts in execution order
│ ├── 01_download.R
│ ├── 02_preprocess.R
│ ├── 03_analysis.R
│ └── 04_figures.R
├── outputs/ # Generated outputs: figures, tables, etc.
│ ├── figures/
│ ├── tables/
│ └── session_info.txt
├── renv.lock # lockfile capturing the R package environment
├── LICENSE
├── CITATION.cff
├── README.md
└── DATA_AVAILABILITY.md


---

## 🧰 Requirements & Setup

### System dependencies

You may need the following non-R system libraries (platform names examples):

- GCC / g++ / gfortran and `make` (for building packages)
- `libcurl` / `libssl` / `openssl`
- `libxml2`
- Compression: `zlib`, `xz`, `bzip2`
- Graphics / rendering: `libpng`, `libjpeg`, `freetype`, `harfbuzz`, `fribidi`
- (Optional, for full runs) SLURM client tools if running on HPC cluster

Example install (Ubuntu / Debian):
```bash
sudo apt-get update
sudo apt-get install -y build-essential gfortran git \
  libcurl4-openssl-dev libssl-dev libxml2-dev \
  zlib1g-dev libbz2-dev liblzma-dev \
  libpng-dev libjpeg-dev \
  freetype-dev libharfbuzz-dev libfribidi-dev
```

### R environment

Clone the repository and restore package versions via renv:

```r
install.packages("renv")
renv::restore()
```

This ensures that you are running with the same package versions used for the analyses.

### 🔍 Data availability

Detailed data access instructions and dataset identifiers are in DATA_AVAILABILITY.md. In summary:

Some input data are available via Synapse (or another controlled repository).

Use the commands shown there (e.g. synapser::synGet(...)) to download required data.

For users without full data access, a small example dataset is included under data/ so that scripts can be run in demonstration mode.

---

## ▶️ How to reproduce the analyses

The scripts should be run in a specific order. Each script contains descriptive comments about its inputs and outputs.

A suggested workflow:

```bash
Rscript scripts/01_download.R
Rscript scripts/02_preprocess.R
Rscript scripts/03_analysis.R
Rscript scripts/04_figures.R
```

After these run, you can inspect outputs in:

outputs/figures/ — figures for the manuscript

outputs/tables/ — intermediate tables / results

outputs/session_info.txt — package and system session details

---

## 🪪 Citation & release

Please cite the manuscript associated with this repository.

A CITATION.cff file is included for machine-readable citations.

After publication, a release (e.g. v1.0-paper) will be created, and a DOI minted via Zenodo. The DOI badge will appear above here.

--- 

## 📜 License

This code is released under the [MIT License](). Unless otherwise noted, reuse is permitted under those terms.