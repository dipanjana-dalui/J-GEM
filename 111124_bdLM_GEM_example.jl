# birth-death logistic model GEM example


## PART OF THE MODULE
using DifferentialEquations
using Plots
using DataFrames
using Distributions
using Random

include("AppendState.jl")
include("CalcAvgFreq.jl")
include("DrawNewTrait.jl")
include("InitiatePopulation.jl")
include("MedianCi.jl")
include("PickEvent.jl")
include("PickIndividual.jl")


Random.seed!(42)  # use only when debugging 

# Load required scripts: We want to make sure that functions that will be called inside the GEM 
# are correctly loaded. Double check if you have asynchronous version numbers across dependent files. 
GEM_type = Array{Int64}([1 2 3])
num_rep = 100

