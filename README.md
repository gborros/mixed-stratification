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
| `SamplingStrata/` | Code used to run the SamplingStrata method for the GHS application |

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

## Running the Code

### HPC

All code was run on the UCT HPC cluster (Ada). Simulated data results were parallelised across 40 cores to get a sense of variability across 40 runs. Application results were run over 10 runs. To reproduce locally, this code can be minimally edited by adjusting the number of cores and parallel runs (`iter`).

### `local_testing/`

For a quicker, non-parallelised sanity check, `local_testing/` contains scripts that let a user run each method for a single specified seed and inspect the outcome directly, without needing HPC access. This is useful for debugging a method or checking output structure before committing to a full 40-core HPC run.

---

## SamplingStrata Folder

`SamplingStrata/` contains the scripts used to run the `SamplingStrata` method for the GHS application, across both the continuous and atomic variants, using different CV inputs. Naming follows `ss_ghs_<method>_<CV input>[_fpc]`, plus a set of Bethel allocation scripts that use MED-GA-derived stratifications and CVs.

| Script | Description |
|---|---|
| `ss_ghs_cont_1` | Continuous SamplingStrata method, with manually specified 0.05 CVs as input |
| `ss_ghs_atomic_1` | Atomic SamplingStrata method, with manually specified 0.05 CVs as input |
| `ss_ghs_cont_2` | Continuous SamplingStrata method, with MED-GA CVs as input |
| `ss_ghs_cont_2_fpc` | Continuous SamplingStrata method, with MED-GA (FPC-adjusted) CVs as input |
| `ss_ghs_atomic_2` | Atomic SamplingStrata method, with MED-GA CVs as input |
| `ss_ghs_atomic_2_fpc` | Atomic SamplingStrata method, with MED-GA (FPC-adjusted) CVs as input |
| `bethel_ghs_cont` | Bethel allocation method, using MED-GA stratification and CVs (MED-GA inputs derived from `ss_ghs_cont_1` inputs) |
| `bethel_ghs_atomic` | Same as `bethel_ghs_cont`, but using atomic-derived inputs |
| `bethel_ghs_cont_fpc` | Bethel allocation method, using MED-GA stratification and FPC-adjusted CVs (MED-GA inputs derived from `ss_ghs_cont_1` inputs) |
| `bethel_ghs_atomic_fpc` | Same as `bethel_ghs_cont_fpc`, but using atomic-derived inputs |

---

## `medoids_mixed_gasampsi/` Folder

Contains the scripts used to run the MED-GA method for the GHS application, using strata numbers and/or sample size inputs sourced from the `SamplingStrata` and Bethel runs above.

| Script | Description |
|---|---|
| `medmix_ghs_atomic` | MED-GA method, using `ss_ghs_atomic_1` inputs |
| `medmix_ghs_cont` | MED-GA method, using `ss_ghs_cont_1` inputs |
| `medmix_ghs_bethel` | MED-GA method, using Bethel sample size inputs and `ss_ghs_cont_1` strata number inputs |
| `medmix_ghs_bethel_fpc` | MED-GA method, using Bethel sample size inputs (where Bethel was supplied with FPC-adjusted CVs) and `ss_ghs_cont_1` strata number inputs |
| `medmix_ghs_bethel_atomic` | MED-GA method, using Bethel sample size inputs and `ss_ghs_atomic_1` strata number inputs |
| `medmix_ghs_bethel_fpc_atomic` | MED-GA method, using Bethel sample size inputs (where Bethel was supplied with FPC-adjusted CVs) and `ss_ghs_atomic_1` strata number inputs |
| `run_medmix_gasampsiMIX1` – `run_medmix_gasampsiMIX6` | Scripts to run the MED-GA method on the simulated datasets |

---
