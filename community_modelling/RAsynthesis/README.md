# Case study: Rosmarinic acid synthesis

The notebooks and scripts in this directory were used to replicate results from [Balancing the non-linear rosmarinic acid biosynthetic pathway by modular co-culture engineering](https://doi.org/10.1016/j.ymben.2019.03.002) using static and dynamic metabolic modelling methods.

## Content overview

All genome-scale metabolic models generated for the strains in the case study can be found in the [GEMs/](../RAsynthesis/GEMs/) folder.

Code for attempts that were abandonned are located in the [archive/](../RAsynthesis/archive/) folder -- should be ignored.

Scripts and functions that are used in the final analysis can be found in the [functions](../RAsynthesis/functions/) folder.
Most central are [create_GEMS.py](create_GEMS.py), [dfba_comets.py](functions/dfba_comets.py) and [plot_results.py](functions/plot_results.py).

### Functions overview

| File                          | Description             |
|-------------------------------|-------------------------|
| [create_GEMS.py](create_GEMS.py) | Functions to create all GEMs for the case study |
| [functions/all_exp_pairwise.py](functions/all_exp_pairwise.py) | Script for running 36 experiments with all inoculation and glucose-to-xylose -ratio combinations for CAL11, SAL11, and MAM3 tri-culture |
| [functions/data_analysis.py](functions/data_analysis.py) | Functions to processing the experimental data retreievd from case study publication|
| [functions/dfba_comets.py](functions/dfba_comets.py) | Functions to run dFBA simulations for consortia |
| [functions/inoc_ratio_coculture_script.py](functions/inoc_ratio_coculture_script.py) | Script for running 7 experiments with all inoculation ratio combinations for RAU2:RAD4 co-culture |
| [functions/inoc_ratio_glc_script.py](functions/inoc_ratio_glc_script.py) | Script for running 16 experiments with all inoculation ratio combinations for CAL2, SAL9, and MAM2 tri-culture |
| [functions/modify_GEM.py](functions/modify_GEM.py) | Function to add flux coupling constraint in COBRApy |
| [functions/plot_results.py](functions/plot_results.py) | Functions to make all plots in thesis |
| [functions/steadiercom.py](functions/steadiercom.py) | An updated version of reframed.community.simulate.py, extending steadiercom to include variance analysis. |
| [functions/updated_steadycom.py](functions/updated_steadycom.py) | Extension of SteadyCom function from reframed |

### Results overview

| Notebook                                 | Description                                |
|------------------------------------------|--------------------------------------------|
| [results_coculture.ipynb](../RAsynthesis/results_coculture.ipynb)                 | RAD2 and RAU4 coculture simulation             |
| [results_glc_tri.ipynb](../RAsynthesis/results_glc_tri.ipynb)                       | CAL2, SAL9, and MAM2 triculture simulation       |
| [results_mixed_sub_tri.ipynb](../RAsynthesis/results_mixed_sub_tri.ipynb) | CAL11, SAL11, and MAM3 triculture simulation     |