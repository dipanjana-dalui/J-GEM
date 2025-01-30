# birth-death logistic model GEM example


## PART OF THE MODULE
using DifferentialEquations
using Plots
using DataFrames
using Distributions
using Random
using Statistics
using StaticArrays
using BenchmarkTools
using StatsBase

#include("AppendState.jl")

include("CalcAvgFreq.jl")
include("DrawNewTrait.jl")
include("InitiatePopulation.jl")
include("MedianCi.jl")
include("PickEvent.jl")
include("PickIndividual.jl")


Random.seed!(42)  # use only when debugging 

#= *** DESIGN CHOICES *** =#


# Major decision time: here you will decide on how many GEM types you want to run
# in total. You can also choose to run the GEM types separately.

# For example, I will create a GEM type array of Integer 1, 2, 3: these are our
# 3 types of GEMS.

GEM_type = Array{Int64}([1 2 3])
GEM_ver = Array{Int64}([1 2 3]) ## another upper level loop is needed

num_rep = 5

h2_levels = Array{Float64}([0, 0.25, 0.5])
cv_levels = Array{Float64}([0, 0.2])

#= *** INITIAL CHOICES *** =#

# define initial poopulation 
R0 = Array{Float64}([10.0]) #we picked the initial population 
# R0 = Array{Float64}([10.0, 20.0]) #for two species initial population

# bd-logistic parameters distribution  
b_max_mu = 4.0
b_max_sigma = 0.01

b_s_mu = 0.0012
b_s_sigma = 0.0

d_s_mu = 0.001
d_s_sigma = 0.0

d_min_mu = 1
d_min_sigma = 0.0

b_max = rand(LogNormal(log(b_max_mu), b_max_sigma), 1)  # max birth
d_min = rand(LogNormal(log(d_min_mu), d_min_sigma), 1) # min death
b_s = rand(LogNormal(log(b_s_mu), b_s_sigma), 1) # density dependence of birth
d_s = rand(LogNormal(log(d_s_mu), d_s_sigma), 1) # density dependence of death

params = [vec(b_max)[1], vec(d_min)[1], vec(b_s)[1], vec(d_s)[1]]

# calculate initial constant 
r_max = b_max-d_min
K = floor(vec((b_max - d_min)/(b_s + d_s))[1])

no_species = length(R0) ## also, no_species = size(state_par_match, 1) 
no_params = length(params)  ## also, size(state_par_match, 2)

state_par_match = SMatrix{no_species, no_params}(1, 1, 1, 1)
state_geno_match = SMatrix{no_species, no_params}(0, 0, 0, 0)
geno_par_match = SMatrix{no_species, no_params}(0, 0, 0, 0)

which_par_quant = state_par_match - geno_par_match
no_columns = no_params + 1 + size(state_geno_match, 2) 

t_max = 20.0 # we will keep it low for checking purpose
min_time_step_to_store = 1

#= fix later
cv_vect = SMatrix{no_params,no_species}(0.2, 0, 0, 0)
size(cv_vect)
h2_vect = SMatrix{2,no_species}(0, 0.25) 
size(h2_vect)
or, 
h2_vect = [0.2 0 0 0] # row corresponds to state, 
                      #col corresponds to trait heri
=#
h2 = [0.2] #1xno_states; col parameter later decides which h2 value
            # will get multiplied to parent traits 
######################################################
######################################################
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
     title="solution to the logistic growth equation",
     xaxis="Time (t)", yaxis="R(t)",
     label="K=500") 

# All looks good
######################################################

for j = i:length(GEM_type)
    t = 0
    stand_time = range( 0, t_max, step = min_time_step_to_store)
    num_time_steps = length(stand_time)

    ## memory preallocation should happen here for the entire simulations
    ## these are all the data structure we are storing
    pop_stand = zeros(no_species, num_time_steps, num_rep )
    x_stand = fill(NaN, no_columns-1,num_time_steps, no_species,num_rep) #trait
    x_var_stand = fill(NaN, no_columns-1,num_time_steps, no_species, num_rep) # trait variance 

    pop_data_out = fill(NaN, 3, num_time_steps, no_species)

    for i = 1:num_rep
        t = 0
        #save a copy inside the loop to update
        R = R0
        ## all structures made here are temp and exist only inside the
        ## replication loop
        #init_comm_mat =  MMatrix{R0, no_columns}(fill(NaN, comm_mat_rows, no_columns))
        init_comm_mat =  Array{Float64}(fill(NaN, Int(sum(R0)), no_columns))

        pop_slice = Array{Float64}(fill(0.0, no_species, num_time_steps))
        x_slice = fill(NaN, no_columns-1, num_time_steps, no_species)
        x_var_slice = fill(NaN, no_columns-1, num_time_steps, no_species)

        # store details of initial state 
        pop_slice[:,1] = R0 #first col/time step gets the initial pop 
        

        #draw initial population: 10x9 Matrix{Float64}
        x_dist_init = InitiatePop(R0, which_par_quant, state_geno_match, state_par_match,
            init_comm_mat, params, cv_vect) ## sum(R0) number of rows and no_col num of cols
            
        for ii = 1:no_species #8×21×1 Array{Float64, 3}:[:, :, 1] 
            x_slice[:, 1, ii] = CalcAverageFreqs(ii, no_columns, no_params, x_dist_init)
            x_var_slice[1:no_params,1,ii] = var(x_dist_init[x_dist_init[:,1] .== ii,2:no_params+1],dims=1)
        end

        # count up each individual for all states
        for jj = 1:length(R)
                x = init_comm_mat[:,1] #extract first col
                #typeof(x)
                R[jj] = count(.==(jj), x)
                #R[jj] = count(x -> (x.==jj), init_comm_mat)
        end # this is suppose to give the count of individuals in each population


        x_dist = x_dist_init
        time_step_index = 2
        time_step = stand_time[time_step_index]

        while t < t_max && sum(R) >= 1 # <--- MAKE SURE THIS IS THE RIGHT R VALUE - 
                                  # perhaps we do need the the jj loop l153-159
            
            # Find who is next
            FindWhoNext = WhoIsNext(x_dist,no_species,no_columns,no_params,y0,
                                          state_par_match,state_geno_match)
            params_next = FindWhoNext[1]
            genotypes_next = FindWhoNext[2]
            whosnext = FindWhoNext[3]
            # R birth
            b_max = params[1] # max birth
            d_min = params[2] # min death
            b_s = params[3] # density dependence of birth
            d_s = params[4] 

            birth_R =  (b_max - b_s*R[1])*R[1]
            death_R =  (d_min - d_s*R[1])*R[1]

            terms = [birth_R, death_R]

            [c_sum, row, col] = PickEvent(terms, no_species)

            #c_sum = [10 20]
            #row = 1 #event
            #col = 1 #state
            ## row = 1 && col == 1 - sp1 birth
            ## row = 2 && col == 2 - sp1 death
            ## row = 1 && col == 1 - sp2 birth
            ## row = 2 && col == 2 - sp2 death 
            if row == 1
                parent_traits = x_dist[Int(whosnext[col]), 2:no_columns] 
                        
                new_trait = DrawNewTraits(x_dist,parent_traits,h2_vect,no_params,no_columns,col)
                # new_trait = offspring_traits, offspring_genotypes
                # ([4.676604793653037 1.0 0.0011999999999999995 0.0010000000000000002], [NaN NaN NaN NaN])

                new_trait_row = hcat(col, new_trait)
                x_dist = vcat(x_dist, new_trait_row) 
            elseif row == 2 # death 
                # delete the individual by making a new copy of the matrix
                # without the row 
                x_dist = x_dist[1:size(x_dist, 1) .!= Int(whosnext[col]), :]
            end

            ## UPDATE ABUNDANCES 
            #test_x_dist = hcat(x_dist, )
            for jj in 1:no_species
                R[jj] = sum(x_dist[:,1].== jj)
            end

            if t > time_step
                pop_slice[1:no_species,time_step_index] .= R # assign current values to sliced standard times
                
                for ii in 1:no_species
                    x_slice[:,time_step_index,ii] = CalcAverageFreqs(ii,no_columns,no_params,x_dist)
                    x_var_slice[1:no_params,time_step _index,ii] = var(x_dist[x_dist[:,1]==ii,2:no_params+1],1)
                end
                time_step_index = time_step_index + 1 # advance to next standardized time
                time_step = stand_times(time_step _index)
                =#
            end

            #last thing to do before exiting the loop:
            # ADVANCE TIME 
            time_advance = exp(-1/c_sum[end])/(c_sum[end])
            if isnan(time_advance) == 0
                t = t + time_advance
            else
                break
            end

        end #end of while loop running the core GEM

        # store the last value of the replicate``
        pop_slice[1:no_species, time_step_index] = R
        for ii = 1:no_species
              x_slice[:,time_step_index,ii] = CalcAverageFreqs(ii,no_columns,no_params,x_dist);
              x_var_slice[1:no_params,time_step_index,ii] = var(x_dist[x_dist[:,1].== ii,2:no_params+1],dims=1);
        end

        ## save the pop and param/geno values for the whole replicate in premade containers  c
        pop_stand[:, :, i] = pop_slice
        x_stand[:, :, :, i] = x_slice 
        x_var_stand[:, :, :, i] = x_var_slice ##

    end ## end of for loop for single simulation replicate

    upper_ci_level = 75
    lower_ci_level = 25

    # pop_stand is no spp x no time step x no replicates
    #=MATLAB
    for ii = 1:no_species # dimension 3 in sequence is number of species
        pop_data_out[1,:,ii] = prctile(pop_stand[ii,:,:],lower_ci_level,3);
        pop_data_out[2,:,ii] = prctile(pop_stand[ii,:,:],50,3);
        pop_data_out[3,:,ii] = prctile(pop_stand[ii,:,:],upper_ci_level,3);
    end
    =#
    .
    x_data_out = MediansCI(upper_ci_level,lower_ci_level,x_stand) # trait
    x_var_data_out = MediansCI(upper_ci_level,lower_ci_level,x_var_stand) # variance in trait 
    
       



end
