# First job is to test for type stability in all our 
# auxiliary functions

# we will first create any object type necessary for the test

# this creates a 1x4 matrix that accepts Float64 as its elements
using DifferentialEquations
using Plots
using DataFrames
using Distributions
using Statistics
using Random
using StaticArrays
using BenchmarkTools

# pick an initial resource density
R0 = 10 #we picked the initial population 
typeof(R0) #Float64 


# bd logistic parameters for resource 
b_max_mu = 4.0
b_max_sigma = 0.01

b_s_mu = 0.0012
b_s_sigma = 0.0

d_s_mu = 0.001
d_s_sigma = 0.0

d_min_mu = 1
d_min_sigma = 0.0

# now, to actually draw the value. We have every parameter set like this as 
# this is the most flexible way to turn on selection on any of parameters. Of
# course you can just set a Float64 type value for any parameter that you 
# think is biologically not evolvable.  
b_max = rand(LogNormal(log(b_max_mu), b_max_sigma), 1)
d_min = rand(LogNormal(log(d_min_mu), d_min_sigma), 1)
b_s = rand(LogNormal(log(b_s_mu), b_s_sigma), 1)
d_s = rand(LogNormal(log(d_s_mu), d_s_sigma), 1)
 

# remember to vectorise before taking floor or round
#a = ((b_max - d_min)/(b_s + d_s)) #1494.019278317911
#typeof(a)
#a_val = a[1,1] # accessing the element
#a_val_vec = vec(a)[1] #flattens/vectorizes the array

#K_vec = vec((b_max - d_min)/(b_s + d_s))[1]
#K = floor(K_vec)
r_max = b_max-d_min
K = floor(vec((b_max - d_min)/(b_s + d_s))[1])


# here you have to be careful with the number of states 
# and the number of parameters. We will have as many rows 
# number of states, and 1+parameters number of columns.
# You can define a static array as follows- {2,3} is the size.
# first arg is rows = number of states 
# second arg is cols = number of param
# b = SMatrix{2,3}(1, 2, 3, 4, 5, 6)

state_par_match = SMatrix{1, 4}(1, 1, 1, 1) ## make static array?
typeof(state_par_match)

state_geno_match = SMatrix{1, 4}(0, 0, 0, 0)
#@btime state_geno_match = SMatrix{1, 5}(0, 0, 0, 0, 0)

geno_par_match = SMatrix{1, 4}(0, 0, 0, 0)

# you can check time difference if you run the following two lines
#@btime geno_par_match = SMatrix{1, 5}(0, 0, 0, 0, 0)
#@btime geno_par_match = Array{Int64}([0 0 0 0 0])

which_par_quant = state_par_match - geno_par_match

#no_species = size(state_par_match, 1)
# Here we will make two param data-structures. We will use the
# mutable version of static arrays MVector
# n_par = size(state_par_match, 2)+1  
n_par = size(state_par_match, 2)
params = MVector{n_par}(fill(NaN, n_par))

#params[1] = NaN
params[1] =  vec(b_max)[1]; # max birth
params[2] = vec(d_min)[1]; # min death
params[3] = vec(b_s)[1]; # density dependence of birth
params[4] = vec(d_s)[1]; # density dependence of death
params
# params has one element more than the number of params, hence
# we take care to save no_params seperately
no_species = size(state_par_match, 1)
no_params = size(state_par_match, 2)
no_columns = no_params + 1 + size(state_geno_match, 2) 

## the _match onjects are arrays, size(array, 2) gives me the number of cols


t_max = 20.0 # we will keep it low for checking purpose
min_time_step_to_store = 1

GEM_type = 1
num_GEM_var = 1 #let's look at an example where b_max is under selection
num_rep = 5 # number of GEM simulations/version

params .* cv_vect
cv_vect = [0.2, 0, 0, 0]
h2_vect = [0 h] 

##############################################
##
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


##############################################
### We have all our predefintions, now we have to stich together
### the algorithm 

#You will set up your parallelizing here before you start your 
# loop over the number of replicates. 
# You might have higher loops that should be placed before/outside 
# the number of replicates for loop. These could correspond to GEM Version 
# and GEM type. Defintion of version and type?

for j = 1:length(GEM_type) #only 1 in this case 
    ## add parallelizing here
    ##
    ##
    for i = 1:num_rep
            t = 0
            # standardize time steps for storing time series
            stand_time = range( 0, t_max, step = min_time_step_to_store) #define t_max and  
            # min_time_step in block 1 or 2. 
            num_time_steps = length(stand_time)
      ##
            # The main thing to remember is some of these will be dynamics
            #   MAIN CHECK HERE: StaticArrays MIGHT NOT BE IDEAL SINCE 
            #   SOME OF THESE STRUCTURES GROW BIGGER (StaticArrays BEST USED FOR 
            #   FOR SMALL IMMUTABLE OBJECTS)

            init_comm_mat =  MMatrix{R0, no_columns}(fill(NaN, R0, no_columns))
            init_comm_mat =  Array{Float64}(fill(NaN, R0, no_columns))
            typeof(init_comm_mat) #check that you have the right type
            # here we are creating an initial community matrix
            # the rows should be the total community size. If yoou have
            # more than one state, remember to add all of them up.
            # Number of cols is 1 more than the total number of 
            # parameters and genotypes.

            # now we will create arrays that store the data for each replication
            pop_slice = Array{Float64}(fill(0.0, no_species, num_time_steps))
            # we will keep it a M-Matrix too remain general, but for one state M-vector also works
            typeof(pop_slice) #-- check to see elements are Float64

            x_slice = fill(NaN, no_columns-1, num_time_steps, no_species)
            # since we only have 1 species/state here,  we essentially have a matrix
            # but defining it this way makes it flexible for expansion.
            # Here we are creating a matrix whose rows are total number of param+genotype
            # and the cols are standardised time steps we want to save our data at.
            # The last arg is the number of species - each gets one slice/matrix
            x_var_slice = fill(NaN, no_columns-1, num_time_steps, no_species)

            # let's assign starting values 
            pop_slice[:,1] = [R0]
            # note that if you have more states, you will replace 0's in the
            # first col with all the starting values. 
            #  pop_slice[:,1] = [Y0] <- Y0 is a vector of starting values
            R = R0 #save a copy inside the loop to update
            
            x_dist_init = InitiatePop(y0, which_param_quant, state_geno_match, state_par_match,
            init_comm_mat, params, cv_vect)

            # We will sample initial values and frequencies
            
            for ii = 1:no_species
                  #x_slice = fill(NaN, no_columns-1, num_time_steps, no_species)
                  x_slice[:, 1, ii] = CalcAveragesFreqs(ii, no_columns, no_params, x_dist_init)
                  x_var_slice[1:no_params,1,ii] = 
                  var(x_dist_init[x_dist_init[:,1] .== ii,2:no_params+1],dims=1)
            end

            ## 8
            
            #=
            A = [1, 6]
            x_dist_init[:,1]  
            temp = fill(2.0, 8)   
            test_arr = vcat(x_dist_init[:,1], temp)
            for jj = 1:2
            A[jj] = count(==(jj), test_arr)
            end
            A
            x=test_arr2[:,1]
            test_arr2 = hcat(test_arr, fill(3, length(test_arr)))
            count(x->(x==2), test_arr2)
            =#

            # setting up init values for the loop
            for jj = 1:length(R)
                  x = init_comm_mat[:,1] #extract first col
                  R[jj] = count(x -> (x==jj), init_comm_mat)
            end

            x_dist = x_dist_init
            time_step_index = 2
            time_step = stand_time[time_step_index]

            while t < t_max && sum(R)>0
                  # FIND WHO IS NEXT
                  ## FLAG 1
                  [params_next,genotypes_next,whosnext] = FindWhosNext(x_dist,no_species,no_columns,no_params,y0,state_parameter _match,state_genotype _match);

                  # R birth
                  b_max = params[1] # max birth
                  d_max = params[2] # min death
                  b_s = params[3] # density dependence of birth
                  d_s = params[4] 

                  birth_R =  (b_max - b_s*R[1])*R[1]
                  death_R =  (d_min - d_s*R[1])*R[1]

                  terms = [birth_R, death_R]

                  [event_index, c_sum, row, col] = PickEvent(terms, no_species) ## FLAG 2

                  
            
            