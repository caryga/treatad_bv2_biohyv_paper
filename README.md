# TREAT-AD -- Bioinformatics Hypothesis Validation in BV2 cells

This repository contains code for the analyses and figures in:

> Integrated Phenotypic and Proteomic Screening Identifies Top-Tier Alzheimer's Disease Therapeutic Targets
> Authors: Cary GA, Li Q, Wiley JC, Paisie CA, Du Y, Zoeller EL, Duong D, Fu H, Seyfried NT, Levey AI, Betarbet R, Carter GW for The Emory-Sage-SGC-JAX TREAT-AD Center
> Journal: _in progress_  
> DOI: _in progress_ 

This is the minimal, paper-focused version of the project, containing the scripts and data pointers necessary to reproduce the main results.

---

## 🧰 Requirements & Setup

### System dependencies

You may need the following non-R system libraries (platform names examples):

- GCC / g++ / gfortran and `make` (for building packages)
- `libcurl` / `libssl` / `openssl`
- `libxml2`
- Compression: `zlib`, `xz`, `bzip2`
- Graphics / rendering: `libpng`, `libjpeg`, `freetype`, `harfbuzz`, `fribidi`

### R environment

Clone the repository and restore package versions via renv:

```r
install.packages("renv")
renv::restore()
```

This ensures that you are running with the same package versions used for the analyses.

### 🔍 Data availability

Detailed data access instructions and dataset identifiers are in DATA_AVAILABILITY.md. In summary:

Most input data are available via Synapse.

Use the commands shown there (e.g. `get` from the `synapseclient` python package) to download required data.

---

## ▶️ How to reproduce the analyses

1) To analyze the cell phenotype data, run `bv2_cell_assays_lmm.R`
2) To analyze the proteomics data, run `bv2_proteomics.R`
3) To collect microglial datasets, run `microglial_cell_sets.R`
4) To generate the figures, using the output from the scripts in steps #1-3, follow the `bv2_biohyv_manuscript_updated.Rmd` notebook

---

## 🪪 Citation & release

Please cite the manuscript associated with this repository.

--- 

## 📜 License

This code is released under the [MIT License](). Unless otherwise noted, reuse is permitted under those terms.
