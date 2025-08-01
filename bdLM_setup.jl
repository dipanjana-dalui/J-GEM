#= Julia-GEM setup file
This is your set up file.
Anything that is preceeded by #=*** X X ***=# is for user to change
Anything that is preceeded by #### need not be changed

=#

include("functions/Packages.jl")
include("functions/AuxiliaryFunc.jl")

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
N0 = Vector{Int64}([10.0]) #we picked the initial population 
# N0 = Array{Float64}([10.0, 20.0]) #for two species initial population

# bd-logistic parameters distribution  
b_max_mu = 4.0
b_max_sigma = 0.0

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

param_init = [vec(b_max)[1], vec(d_min)[1], vec(b_s)[1], vec(d_s)[1]]
# par_names = ("b_max", "d_min", "b_s", "d_s")

# calculate initial constant 
r_max = b_max-d_min
K = floor(vec((b_max - d_min)/(b_s + d_s))[1])

no_species = length(N0) ## also, no_species = size(state_par_match, 1) 
no_params = length(param_init)  ## also, size(state_par_match, 2)

include("bdLM_bdTerms.jl")

#=***************************************************
**                  DESIGN CHOICES                 **
***************************************************=#
# Major decision time: here you will decide on how many GEM vers you want to run
# in total. You can also choose to run the GEM vers separately.

# For example, I will create a GEM ver array of Integer 1, 2, 3: these are our
# 3 vers of GEMS.

GEM_ver = Vector{String}(["ver1","ver2"])
num_rep = 10
t_max = 15.0 # we will keep it low for checking purpose
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

h2_vect = [0.0;
            0.2]
              ## rows: GEM versions, cols: state ID
h2_mat = reshape(h2_vect, length(GEM_ver), no_species)

cv_vect = [0.0;
            0.2]
cv_mat = reshape(cv_vect, length(GEM_ver), no_species)
#cv_vect = collect(transpose([0.2 0 0 0 0 0 0;
#                             0 0 0 0 0 0 0])) #### THIS IS WRONG???
                             ## row = state
                             ## col = GEM version
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

par_names = ["b_max", "d_min", "b_s", "d_s"]
geno_names = ["g_1", "g_2", "g_3", "g_4"]

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
#include("functions/bdLM_v2_GEM_main.jl")

GEM_run = GEM_sim(GEM_ver, 
                  t_max,
                  no_species,
                  num_rep,
                  no_columns,
                  no_params,
                  N0,
                  which_par_quant,
                  state_geno_match,
                  state_par_match,
                  param_init,
                  cv_mat,
                  h2_mat,
                  par_names,
                  geno_names,
                  pop_stand_out_all
                )


# The GEM_sim function returns the population time series, and the parameter mean/variances
pop_time_series_df = GEM_run[1]
#CSV.write("Pop_Time_Series.csv",pop_time_series_df)

trait_mean_df = GEM_run[2]
#CSV.write(trait_mean_df)

trait_var_df = GEM_run[3]
#CSV.write(trait_var_df)

Pop_Plot(pop_time_series_df, 1)

Trait_Plot(trait_mean_df, trait_var_df, 1, "b_max")
