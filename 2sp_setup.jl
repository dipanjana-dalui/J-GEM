#= Julia-GEM setup file
This is your set up file.
Anything that is preceeded by #=*** X X ***=# is for user to change
Anything that is preceeded by #### need not be changed
=#

######################################################
##     AUXILIARY FUNCTIONS AND REQUIRED PACKAGES    ##
######################################################
include("functions/Packages.jl")

include("functions/AuxiliaryFunc.jl")

#=*** set seed ***=#
Random.seed!(42)  # use only when debugging 

#=***************************************************
**         DEFINE MODEL & BIRTH-DEATH FUNC         **
***************************************************=#
#= 
Use this space to write down your model and birth/death terms
dH_1(t)/dt = r1*H1*(1-alpha*H1) - sp*H1*P1
dP_1(t)/dt = sp*f*H1*P1 - m1*P1 

list of parameters
r := growth rate of host 
r_i = intrinsic birth sp i - intrinsic death sp i
alpha := 1/K = self regulation
sp := space clearace rate/attack rate 
f := conversion efficiency  
m1 := mortality

birth terms H1: 
birth_H_1 = b_i*H_1

death terms H1:
death_H_1 =  + 

birth term for P1:
birth_P = sp*f*H1*P1

death term for P1:
death_P = m1*P1

=#

#=***************************************************
**                INITIAL CONDITION                **
***************************************************=#

# define initial poopulation abundances 
N0 = Vector{Int64}([5, 1]) #we picked the initial population 
# R0 = Array{Float64}([10.0, 20.0]) #for two species initial population

# bd-logistic parameters distribution  
b_max_mu = 0.8
b_max_sigma = 0.0

d_min_mu = 0.4
d_min_sigma = 0.0

#=
b_s_mu = 1e-2
b_s_sigma = 0.0

d_s_mu = 1e-5
d_s_sigma = 0.0
=#

scr_mu = 0.005
scr_sigma = 0

fec_mu = 0.05
fec_sigma = 0

m_mu = 0.01
m_sigma = 0

b_max = rand(LogNormal(log(b_max_mu), b_max_sigma), 1)  # max birth
d_min = rand(LogNormal(log(d_min_mu), d_min_sigma), 1) # min death
#b_s = rand(LogNormal(log(b_s_mu), b_s_sigma), 1) # density dependence of birth
#d_s = rand(LogNormal(log(d_s_mu), d_s_sigma), 1) # density dependence of death

scr = rand(LogNormal(log(scr_mu), scr_sigma), 1)
fec = rand(LogNormal(log(fec_mu), fec_sigma), 1)
m = rand(LogNormal(log(m_mu), m_sigma), 1)

param_init = [vec(b_max)[1], vec(d_min)[1], vec(scr)[1], vec(fec)[1], vec(m)[1]]

# calculate initial constant 
r_max = b_max-d_min
#K = floor(vec((b_max - d_min)/(b_s + d_s))[1])

no_species = length(N0) ## also, no_species = size(state_par_match, 1) 
no_param = length(param_init)  ## also, size(state_par_match, 2)


#=***************************************************
**                  DESIGN CHOICES                 **
***************************************************=#
# Major decision time: here you will decide on how many GEM vers you want to run
# in total. You can also choose to run the GEM vers separately.

# For example, I will create a GEM ver array of Integer 1, 2, 3: these are our
# 3 vers of GEMS.

num_rep = 3
GEM_ver = Vector{String}(["ver1", "ver2"])
t_max = 10.0 # we will keep it low for checking purpose
min_time_step_to_store = 0.5

#= fix later: decide on the proper dimensions of h2 and cv
so you can streamline it to scale to #species and #param 
cv_vect = SMatrix{no_param,no_species}(0.2, 0, 0, 0)
size(cv_vect)
h2_vect = SMatrix{2,no_species}(0, 0.25) 
size(h2_vect)
or, †
h2_vect = [0.2 0 0 0] # row corresponds to state, 
                      #col corresponds to trait heri
=#


### CHECK WITH MANUSCRIPT and MATLAB version
h2_vect = [0.0 0.0; 0.1 0.1] ## rows: GEM versions, cols: state ID

cv_vect = [0.0 0.0; 0.2 0.2] ## rows: GEM versions, cols: state ID

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
state_par_match = Array{Int64}([1 1 0 0 0; 0 0 1 1 1]) #no_col = param_init, no_row = state
state_geno_match = Array{Int64}([0 0 0 0; 0 0 0 0])
geno_par_match = Array{Int64}([0 0 0 0 0; 0 0 0 0 0])

no_columns = no_param + 1 + size(state_geno_match, 2) 
which_par_quant = state_par_match - geno_par_match

par_names = ["b_max", "d_min", "scr", "fec", "m"]
geno_names = ["g_1", "g_2", "g_3", "g_4"]
######################################################
##     STORAGE CONTAINERS FOR SIMULATION OUTPUT     ##
######################################################
stand_time = range(0, t_max, step = min_time_step_to_store);
stand_time = collect(stand_time);
num_time_steps = length(stand_time);

pop_stand_out_all = fill(NaN, no_species, num_time_steps, num_rep, length(GEM_ver));
x_stand_out_all = fill(NaN, no_columns-1,num_time_steps, no_species,num_rep, length(GEM_ver));
x_var_stand_out_all = fill(NaN, no_columns-1,num_time_steps, no_species,num_rep, length(GEM_ver));

#=***************************************************
**           CALL GEM SIMULATION FUNCTION          **
***************************************************=#
#include("functions/2sp_v2_GEM_main.jl")


GEM_run = GEM_sim(GEM_ver, 
                  t_max,
                  no_species,
                  num_rep,
                  no_columns,
                  no_param,
                  N0,
                  which_par_quant,
                  state_geno_match,
                  state_par_match,
                  param_init,
                  cv_vect,
                  h2_vect,
                  par_names,
                  geno_names,
                  pop_stand_out_all
                )

    

# The GEM_sim function returns the population time series, and the parameter mean/variances
pop_time_series_df = GEM_run[1]

#CSV.write("Pop_Time_Series.csv",pop_time_series_df)

trait_mean_df = GEM_run[2]
#CSV.write(trait_mean_df)

trait_mean_df[trait_mean_df.GEM_ver.==2,:]
trait_mean_df[trait_mean_df.GEM_ver.==2,:]



trait_var_df = GEM_run[3]
#CSV.write(trait_var_df)



Pop_Plot(pop_time_series_df, 2)

Trait_Plot(trait_mean_df, trait_var_df, 
          STATEID::INT, GEM_VER::INT, "trait_name")