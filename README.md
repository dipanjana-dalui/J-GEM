# GEM (Julia version)
[G]illespie [E]co-evolutionary [M]odel

one sentence blurb

## Coded by
Julia Version
- Dipanjana Dalui (ddalui2@unl.edu)
  
MATLAB versions (archived elsewhere)
- John DeLong (jdelong2@unl.edu)

## Installation
For now, J-GEM is not a package, and you cannot install it on Julia using Pkg.add.
You will have to either clone, fork or download all the scripts from the GitHub account. 

#### Prerequisites
Julia v1.1x 

VSCode (optional, recommended)

## Documentation
An quick reference user guide can be found [here](). 

For detailed documentation, please see [add link to methods paper]

## Citation
TBA

## License
TBA

## Quick Start Guide
Clone the repository to your local machine.
To set up any GEM analysis, you will be interacting with 2 files:
1. **_xx_setup.jl_** ->  This file will walk you through the various blocks that you will configure for your own analysis.  
2. **_xx_bdTerms.jl_** -> a script of functions for the birth and deaths of each state. This is your book-keeping function, arguably the most important part of the algorithm setup. For births, you will take into account all processes that are leading to an increase in the state's abundance. Similarly, for death, make sure you account for all processes that can lead to individuals leaving the state. 
The prefix **_xx_** indicates the example case. 

## Questions/Comments
- email ddalui2@unl.edu or jdelong2@unl.edu
