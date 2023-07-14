# Introduction
All scripts that were used to generate main text and supplementary figures are contained within this project’s git repository: https://github.com/GarciaLab/OptogeneticDissection.git

This document contains instructions for how to use the scripts in this repository to replicate the results shown in Figures 2-4 of the main text. 

For questions, please contact us at:
- nlammers@berkeley.edu (Nicholas Lammers)
- jiaxi.zhao@berkeley.edu  (Jiaxi Zhao)

# General usage notes
- **Software:** All scripts were written and run using Matlab version 2021b. For best results, use this version of the software. 
- **Installation:** Beyond the installation on Matlab itself, our software does not require any installation steps. 
- **Demo:** The datasets used for the actual analysis presented throughout our work are relatively small, and thus computationally tractable, even on a laptop. As a result, we do not include a separate “demo” dataset, since the user should be able to process our full dataset with relative ease.
- **Replication:** This document includes instructions for replicating the analyses required to generate the principal results depicted in Figures 2-4 of the main text. 
- **Data:** Our git repository includes “raw” MS2 trace datasets in the “data” folder. These datasets are themselves an intermediate result, as they are generated using a previously published pipeline for estimating fluorescent spot intensities from the actual raw image files produced by our microscope. 

# Figure 1
## Overview
All scripts related to data processing for the results presented in Figures 1 and related supplemental figures are contained within the root folder.


## Figure 1G
To generate the figure, please run `figure_optogenetics_import_export.m`.

# Figure 2
## Overview
All scripts related to data processing for the results presented in Figures 2 and related supplemental figures are contained within the root folder.

## Figure 2A and E, Figure S3 and S6: Wildtype dynamics
To generate the figure, please run `figure_repression_wildtype.m`.

## Figure 2B: Single embryo dynamics under different illumination conditions
To generate the figure, please run `figure_repression_modular_illumination_single_embryo.m`.

## Figure 2D, Figure S5: Dynamics under different illumination conditions
To generate the figure, please run `figure_repression_modular_illumination.m`.

## Figure S1: Eve rescue comparison
To generate the figure, please run `supp_figure_eve_rescue.m`.

## Figure S2: Nuclei position calibration
To generate the figure, please run `supp_figure_stripe_position_correction.m`.


# Figure 3
## Overview
All scripts related to data processing for the results presented in Figures 3 and related supplemental figures are contained within the root folder.

## Figure 3 and Figure S7
To generate the figure, please run `figure_reactivation_response.m`.


# Figure 4
## Overview

All scripts related to data processing and inference for the results presented in Figures 4 are contained within the subfolder entitled “burst_analysis_scripts”. Scripts with the prefix “cpHMM” pertain to the binary inference results presented in Figure 4 C. Scripts with the prefix “io” pertain to the non-steady-state results presented in Figure 4D-H. In each case, the scripts are numbered in the order in which they should be run. 

**Before running these scripts, change your working directory to be:**
`.\OptogeneticDissection\burst_analysis_scripts\`

## Figure 4C: cpHMM inference
1. **Results generation:** To generate results that underlie Figure 4C, run the following scripts in order (no changes to parameters should be required):
    - **`cpHMM_inference01_build.m`**: this combines, cleans, and interpolates raw MS2 trace data collected from multiple embryos.
    - **`cpHMM_inference02_run_binary.m`**: this calls the core cpHMM parameter inference functions.
    - **`cpHMM_inference03_compile_binary.m`**: this aggregates across multiple inference replicates to calculate burst parameter averages and uncertainty measures for different Knirps groups of interest. 

2. **Figure generation:** 
    - Navigate out to the main OptogeneticsDissection folder and open **`figure_cpHMM_analysis_binary.m`**.  
    - Note that the “dateString” parameter is set to match the datetime corresponding to the results shown in the paper. If you have re-run the inference results, these will be saved under a new datetime, so be sure to update this parameter accordingly. 
    - The script is organized as a notebook. Execute each section in turn. 
        - The second action generates p value estimates
        - The final section will generate the bar plots from Figure 4C

## Figure 4D-H: cpHMM inference
3. **Results generation:** Navigate to the “burst_analysis_scripts” subfolder.  To generate results that underlie Figure 4 E-H, run the following scripts in order (no changes to parameters should be required):
    - **`io_simulation01_build.m`:** this combines, cleans, and interpolates raw MS2 trace data collected from multiple embryos to build datasets that can be used for parameter sweeps and MCMC inference.
    - **`io_simulation02_run_sweep.m`:** This uses parameter sweeps to identify molecular parameters from Equations 1 and S2 that best our experimental data. 
        - Note that this script draws upon cpHMM inference results. The “DateStr” parameter indicates which cpHMM results folder will be read for this. Change accordingly if you have re-run cpHMM inference and wish to use the new results.

    - **`io_simulation03_conduct_fits_v2.m`:** This script uses MCMC inference to refine the parameter estimates from above and to estimate uncertainty bounds.

4. **Figure generation:** 
    - Navigate out to the main OptogeneticsDissection folder and open **`figure_io_sweep_results_v2.m`**.  
    - Note that the “DateStr” parameter is set to match the datetime corresponding to the results shown in the paper. If you have re-run the inference results, these will be saved under a new datetime, so be sure to update this parameter accordingly. 
    - The script is organized as a notebook. Execute each section in turn to generate the panels and results that correspond to Figure 4 D-H in the main text.




