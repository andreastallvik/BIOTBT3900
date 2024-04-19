# Case study: Rosmarinic acid synthesis

The notebooks and scripts in this directory were used to replicate results from [Balancing the non-linear rosmarinic acid biosynthetic pathway by modular co-culture engineering](https://doi.org/10.1016/j.ymben.2019.03.002) using static and dynamic metabolic modelling methods.

## Content overview

GEMs for the BiGG universal model, and the iHK1487 and iML1515 models are located in the [../GEMs/](../GEMs/) folder.

| Folder                        | Description             |
|-------------------------------|-------------------------|
|[GEMs/](../RAsynthesis/GEMs/)  | All genome-scale metabolic models generated for the strains in the case study |
| [archive/](../RAsynthesis/archive/) | Old code that should please be ignored |
| [exp_data/](../RAsynthesis/exp_data/) | Experimental retrieved from the case study |
| [functions/](../RAsynthesis/functions/) | Scripts and functions that are used in the analysis |
| [results/](../RAsynthesis/results/) | All generated plots in the thesis |

### Functions overview

| File                          | Description             |
|-------------------------------|-------------------------|
| [GEM_eval.py](GEM_eval.py) | Functions used to test validity of GEMs |
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

Consortia simulations are conducted in the following notebooks:

| Notebook                                 | Description                                |
|------------------------------------------|--------------------------------------------|
| [experimental_data.ipynb](experimental_data.ipynb) | Calculations on experimental data | 
| [results_coculture.ipynb](../RAsynthesis/results_coculture.ipynb)                 | RAD2 and RAU4 co-culture simulations             |
| [results_glc_tri.ipynb](../RAsynthesis/results_glc_tri.ipynb)                       | CAL2, SAL9, and MAM2 tri-culture simulations       |
| [results_mixed_sub_tri.ipynb](../RAsynthesis/results_mixed_sub_tri.ipynb) | CAL11, SAL11, and MAM3 tri-culture simulations     |

The most relevant functions used to generate reults are found in [create_GEMS.py](create_GEMS.py), [dfba_comets.py](functions/dfba_comets.py) and [plot_results.py](functions/plot_results.py).
