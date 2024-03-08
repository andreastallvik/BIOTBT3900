# Case-study: Butanol production
Community metabolic modelling of a co-culture of *Thermoanaerobacterium thermosaccharolyticum* M5 and *Clostridium acetobutylicum* NJ4 from the case-study [Consolidated bioprocessing performance of a two-species microbial consortium for butanol production from lignocellulosic biomass
](https://doi.org/10.1002/bit.27464).

## Workflow:
1. Automatic reconstruction of GEMs for M5 and NJ4 (carveme)
1. Model validation and curation (cobrapy)
1. Static modelling of individual strains (cobrapy)
1. Dynamic modelling of individual strains (COMETS)
1. Dynamic modelling of community (COMETS)

## Contents
- /utils folder contains python scrips with helper functions and utility for running simulations 
- notebooks for validation and simulation of each strain
- notebooks for community simulation
