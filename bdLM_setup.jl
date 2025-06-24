#= Julia-GEM setup file

This is your set up file.

Anything that is preceeded by #=*** X X ***=# is for user to change

Anything that is preceeded by #### need not be changed


=#

include("Packages.jl")
include("AuxillaryFunc.jl")



#=*** set seed ***=#
Random.seed!(42)  # use only when debugging 

#=***************************************************
**         DEFINE MODEL & BIRTH-DEATH FUNC         **
***************************************************=#
#= 
Use this space to write down your model and birth/death terms
dN_i(t)/dt = F(r_i,N_i(t),a_ij)

r_i = intrinsic birth sp i - intrinsic death sp i

birth terms Sp_i (N_i): 
birth_N_i = b_i*N_i

death terms Sp_i (N_i):
death_N_i = d_i*N_i + (competition_ij) + (predation_ij)

=#

#=***************************************************
**                INITIAL CONDITION                **
***************************************************=#

# define initial poopulation abundances 
R0 = Vector{Float64}([10.0]) #we picked the initial population 
# R0 = Array{Float64}([10.0, 20.0]) #for two species initial population

# bd-logistic parameters distribution  
b_max_mu = 4.0
b_max_sigma = 0.01

b_s_mu = 0.0012
b_s_sigma = 0.0

d_s_mu = 1e-5
d_s_sigma = 0.0

d_min_mu = 1
d_min_sigma = 0.0

b_max = rand(LogNormal(log(b_max_mu), b_max_sigma), 1)  # max birth
d_min = rand(LogNormal(log(d_min_mu), d_min_sigma), 1) # min death
b_s = rand(LogNormal(log(b_s_mu), b_s_sigma), 1) # density dependence of birth
d_s = rand(LogNormal(log(d_s_mu), d_s_sigma), 1) # density dependence of death

params = [vec(b_max)[1], vec(d_min)[1], vec(b_s)[1], vec(d_s)[1]]
# par_names = ("b_max", "d_min", "b_s", "d_s")

# calculate initial constant 
r_max = b_max-d_min
K = floor(vec((b_max - d_min)/(b_s + d_s))[1])

no_species = length(R0) ## also, no_species = size(state_par_match, 1) 
no_params = length(params)  ## also, size(state_par_match, 2)


#=***************************************************
**                  DESIGN CHOICES                 **
***************************************************=#
# Major decision time: here you will decide on how many GEM vers you want to run
# in total. You can also choose to run the GEM vers separately.

# For example, I will create a GEM ver array of Integer 1, 2, 3: these are our
# 3 vers of GEMS.

GEM_ver = Vector{Int64}([1, 2, 3])
num_rep = 5
t_max = 5.0 # we will keep it low for checking purpose
min_time_step_to_store = 0.1

#= fix later: decide on the proper dimensions of h2 and cv
so you can streamline it to scale to #species and #param 
cv_vect = SMatrix{no_params,no_species}(0.2, 0, 0, 0)
size(cv_vect)
h2_vect = SMatrix{2,no_species}(0, 0.25) 
size(h2_vect)
or, 
h2_vect = [0.2 0 0 0] # row corresponds to state, 
                      #col corresponds to trait heri
=#

h2_vect = [0.0] #1xno_states; col parameter later decides which h2 value
            # will get multiplied to parent traits

cv_vect = [0.2, 0, 0, 0]
 
#=
******************************************************
**          SANITY CHECK: DETERMINISTIC ODE         **
******************************************************
## Okay let's test the parameters in a ODE simulations
# du := the derivative value
# u := integer state value 
# p := parameters (r_max, K)
# t := time
using DifferentialEquations
using Plots; gr() 

r_max = 2
K = 500
f(u,p,t) = r_max * u * (1 - u / K)
u0 = 10.0
tspan = (0.0,10)
prob = ODEProblem(f,u0,tspan)
sol = solve(prob)

plot(sol, linewidth=3,
     title="solution to the logistic equation",
     xaxis="Time (t)", yaxis="R(t)",
     label="K=500") 

=#
######################################################
# If the deterministic model looks good to you, you can 
# now move on to the GEM. We need to done one more thing 
# before we are ready to call the GEM simulation function
# - we have to make some containers to store the simulation
# outputs. 

#=***************************************************
**               PAR & GENOTYPE MATCH              **
***************************************************=#
state_par_match = Array{Int64}([1 1 1 1]) #no_col = params, no_row = state
state_geno_match = Array{Int64}([0 0 0 0])
geno_par_match = Array{Int64}([0 0 0 0])

# might be good practice to write down the order of of the parameters 
# par_names = ("b_max", "d_min", "b_s", "d_s")
# geno_names = ("g_b_max", "g_d_min", "g_b_s", "g_d_s")
which_par_quant = state_par_match - geno_par_match
no_columns = no_params + 1 + size(state_geno_match, 2) 


######################################################
##     STORAGE CONTAINERS FOR SIMULATION OUTPUT     ##
######################################################
stand_time = range(0, t_max, step = min_time_step_to_store);
stand_time = collect(stand_time);
num_time_steps = length(stand_time);

#pop_stand_out = fill(NaN, no_species, num_time_steps, num_rep)
#x_stand_out = fill(NaN, no_columns-1,num_time_steps, no_species,num_rep) #trait
#x_var_stand_out = fill(NaN, no_columns-1,num_time_steps, no_species, num_rep) # trait variance 

pop_stand_out_all = fill(NaN, no_species, num_time_steps, num_rep, length(GEM_ver));
x_stand_out_all = fill(NaN, no_columns-1,num_time_steps, no_species,num_rep, length(GEM_ver));
x_var_stand_out_all = fill(NaN, no_columns-1,num_time_steps, no_species,num_rep, length(GEM_ver));

#=***************************************************
**           CALL GEM SIMULATION FUNCTION          **
***************************************************=#
include("bdLM_v2_GEM_main.jl")

GEM_run = GEM_sim(GEM_ver::Vector{Int64}, 
    t_max::Float64,
    min_time_step_to_store::Float64,
    no_species::Int64,
    num_rep::Int64,
    no_columns::Int64,
    no_params::Int64,
    R0::Vector{Float64},
    which_par_quant::Matrix{Int64},
    state_geno_match::Matrix{Int64},
    state_par_match::Matrix{Int64},
    params::Vector{Float64},
    cv_vect::Vector{Float64},
    h2_vect::Vector{Float64},
    pop_stand_out_all::AbstractArray{Float64, 4} # Output array to be populated
    )

# The GEM_sim function returns the population time series, and the parameter mean/variances
pop_time_series_df = GEM_run[1]
#CSV.write("Pop_Time_Series.csv",pop_time_series_df)

trait_mean_df = GEM_run[2]
#CSV.write(trait_mean_df)

trait_var_df = GEM_run[3]
#CSV.write(trait_var_df)



# Jun 24: everything up to this point works.
# Things to change/discuss with John about changing:
#         - (i)  parallet the replicates
#         - (ii) add a separate birth death function, or 
#                make it part of the setup file?
# OTHER DESIGN CHOICES
# Modules
# struct 
#
#
#
#
#
#