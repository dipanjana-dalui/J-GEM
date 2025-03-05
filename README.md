# USER-GUIDE-TO-J-GEM
 Tutorial and user guide code with examples

#### Files you will need: 

`bdLM_GEM_main.jl` - the main file where you will set up the simulation, and also the core simulation process 
(perhaps break this into two and made the core GEM loop a function to make it more modular?)

Below is the list of auxillary function files you would need to load to your working environment (all prefixed with bd_LM_):
- `WhoIsNext.jl`
- `DrawNewTraits.jl`
- `PickEvent.jl`
- `InitiatePop.jl`
- `PickIndiv.jl`
- `CaclAvgFreq.jl`

Each function file comes with a scract file section that is commented out, but can be used to check and debug each function. 
