# Case-study: Butanol production
Dynamic community metabolic modelling of a co-culture of *Thermoanaerobacterium thermosaccharolyticum* M5 and *Clostridium acetobutylicum* NJ4 from the publication [Consolidated bioprocessing performance of a two-species microbial consortium for butanol production from lignocellulosic biomass
](https://doi.org/10.1002/bit.27464).

## Project workflow

The main steps in this project were:

1. Automatic reconstruction of GEMs (CareveMe)
1. Model validation and curation (COBRApy)
1. Static modelling of individual strains (COBRApy)
1. Dynamic modelling of individual strains (COMETS)
1. Dynamic modelling of community (COMETS)

## Content overview

The GEM for the BiGG universal model, is located in the [../GEMs/](../GEMs/) folder.
The medium used for model reconstruction and simulation is in the [medium.tsv](medium.tsv) file.
The main results are presented in the notebooks [m5_sim.ipynb](m5_sim.ipynb), [nj4_sim.ipynb](nj4_sim.ipynb), and [community_sim.ipynb](community_sim.ipynb).

| Folder | Description |
|--------|-------------|
|[GEMs](../CBP_butanol/GEMs/) | Genome-scale metabolic models constructed for the case study strains, and GEMs from closely related strains |
|[escher_maps/](../CBP_butanol/escher_maps/) | .json files to make escher map visualisations|
|[exp_data/](../CBP_butanol/exp_data/) | Experimental data retrieved from the case study |
|[figures/](../CBP_butanol/figures/) | All figures generated for the thesis|
|[grid_search_results/](../CBP_butanol/grid_search_results/) | Raw results from grid search over fermentation conditions |
|[utils/](../CBP_butanol/utils/) | Functions and scripts for simulations and analysis |


### Overview of utils/

The [utils/](../CBP_butanol/utils/) folder contains functions and scripts to run simulation and analysis, here follows a overview of the content:

| File | Description |
|------|-------------|
|[comets_functions.py](utils/comets_functions.py) | Functions to run single and sequential dFBA simulations |
|[data_processing.py](utils/data_processing.py) | Functions to clean the experimental data retrieved from case study publication |
|[fermentation_optimisation.py](utils/fermentation_optimisation.py) | Script to run grid search over inoculation time, inoculation ratio, and xylan concentration |
|[flux_coupling.py](utils/flux_coupling.py) | Function to add flux coupling constraints |
|[kinetic_params.py](utils/kinetic_params.py) | Dictionary of kinetic parameter values |
|[m5_modifications.py](utils/m5_modifications.py) | Modifications for curation of M5 model |
|[model_validation.py](utils/model_validation.py) | Helper functions for model validation |
|[nj4_modifications.py](utils/nj4_modifications.py) | Modifications for curation of NJ4 model |
|[read_excel.py](utils/read_excel.py) | Functions to read GEMs from Excel files |
|[static_sim.py](utils/static_sim.py) | Helper functions to run and plot FBA and FVA simulations |
|[unit_conversion.py](utils/unit_conversion.py) | Functions to perform unit conversions |


### Notebooks overview

| Notebook | Description |
|----------|-------------|
|[community_sim.ipynb](community_sim.ipynb) | Results for M5 and NJ4 community simulation|
|[data.ipynb](data.ipynb) | Inspection and calculations on experimental data from publication |
|[m5_sim.ipynb](m5_sim.ipynb) | Results from the M5 monoculture simulation |
|[m5_validation.ipynb](m5_validation.ipynb) | Used to validate and curate the M5 GEM |
|[nj4_sim.ipynb](nj4_sim.ipynb) | Results from the NJ4 monoculture simulation |
|[nj4_validation.ipynb](nj4_validation.ipynb) | Used to validate and curate the NJ4 GEM |