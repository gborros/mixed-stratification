# Mixed Stratification #
## Created by: Georgi Borros ##
## 5 July 2026 ##

## Code replication ##
medoids_mixed_gasampsi: folder containing the code used to run the medoids for stratification with GA allocation method for each dataset
Proportional: code used to run the medoids for stratification with proportional allocation
SamplingStrata: code used to run the SamplingStrata method fo each dataset

## Functions ##
core_functions: the key functions used to perform the methods 
  - function_calc_variance: used to calculate measures of variation. In this case it is the design effect and coefficient of variation.
  - function_sample_size: function used to calculate proportional allocation.
  - medmix_gasampsi: the genetic algorithm function used by medoids_mixed_gasampsi.

## Data ##
multivar_datasets: the datasets used. 
dataset_generation2026.R: the script used to create the simulated datasets.

## CVs ##
cv_extract: extracts the CVs from the results after running the medoids' scripts. 
CVs folder: contains the CVs corresponding to the minimum design effect from the medoids' methods, for input into the SamplingStrata runs. 

### Key note ###
All code was run on the UCT HPC cluster, parallelised across 40 cores for replication. This code can be minimally edited for local reproduction by changing the number of cores and parallel runs (iter). 
