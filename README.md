# Mixed Stratification

**Created by:** Georgi Borros
**Date:** 5 July 2026

## Overview

This repository contains the code and data used to replicate a comparison of stratification methods: medoid-based stratification with genetic algorithm (GA) sample allocation, medoid-based stratification with proportional allocation, and the `SamplingStrata` method.

---

## Code Replication

| Folder/File | Description |
|---|---|
| `medoids_mixed_gasampsi/` | Code used to run the medoids for stratification with the GA allocation method, for each dataset |
| `Proportional/` | Code used to run the medoids for stratification with proportional allocation |
| `SamplingStrata/` | Code used to run the SamplingStrata method for each dataset |

---

## Functions

### `core_functions/`

Key functions used to perform the methods:

- **`function_calc_variance`** — Calculates measures of variation: the design effect and the coefficient of variation.
- **`function_sample_size`** — Calculates proportional allocation.
- **`medmix_gasampsi`** — The genetic algorithm function used by `medoids_mixed_gasampsi`.

---

## Data

- **`multivar_datasets/`** — The datasets used for replication.
- **`dataset_generation2026.R`** — Script used to generate the simulated datasets.

---

## Coefficients of Variation (CVs)

- **`cv_extract`** — Extracts the CVs from the results after running the medoids scripts.
- **`CVs/`** — Contains the CVs corresponding to the minimum design effect from the medoids methods, for input into the SamplingStrata runs.

---

## Working Folders

- **`RAW26`** — Once results are run, they should be copied to the RAW26 folder [not uploaded here due to space constraints], where key results are extracted to *WIP26* using the extraction functions in this folder. 
- **`WIP26`** — Contains extracted results from *RAW26*. 
---

## Key Note

All code was run on the UCT HPC cluster, parallelised across 40 cores to get a sense of variability across 40 runs. To reproduce locally, this code can be minimally edited by adjusting the number of cores and parallel runs (`iter`).
