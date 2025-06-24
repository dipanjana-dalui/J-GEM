##############################################
#		FUNCTION INITIATE POPULATION         #
##############################################
function InitiatePop(y0, which_par_quant, state_geno_match, state_par_match,
	init_comm_mat, params, cv_vect)
	y0 = R0
	end_row = cumsum(y0)
	starting_row = [1; 1 .+ end_row[1:length(end_row)-1]]

	#traits_to_assign = which_par_quant .* y0 # CHECK
	
	# you want each row multiply each row with each corresponding state
	# Matlab should do y0'. Julia's default is a col vector
	gts_to_assign = state_geno_match .* y0
	num_gts = size(state_geno_match, 2)
	for qq = 1:length(y0)
		init_comm_mat[Int(starting_row[qq]):Int(end_row[qq]), 1] .= qq
		for zz = 1:length(params)
			#zz = 1
			traits_to_assign = state_par_match .* y0
			# SINCE BOTH PARAM AND CV_VECT HAVE 4 ROWS AND 1 COL - JULIA IS BY DEFAULT A COL VEC 
			# WE WILL SWAP THE INDICES TO BE (zz, qq)
			temp = PickIndiv(params[zz,qq],cv_vect[zz, qq]*params[zz, qq],Int(traits_to_assign[qq,zz])) 

			if !isempty(temp) #when temp has a non-zero value
				init_comm_mat[Int(starting_row[qq]):Int(end_row[qq]), 1+zz] .= temp
			end
		end

		# ALL?
		if all(gts_to_assign[qq,:] .> 0)
			genotype = zeros(y0[qq], num_gts) 
			# makes a matrix with communuty size no. of rows, and number of geno cols
			
			# go row by row and randomly choose genotype to set to 1
			for yy = 1:y0[qq]
				temp = rand(1:num_gts) # try this alt temp = rand(1:num_gts, 1)
				genotype[yy, temp] = 1
			end			
		
		init_comm_mat[Int(starting_row[qq]):Int(end_row[qq]), 2 + length(params):1+length(params)+num_gts] = genotype
		end
	end
	return init_comm_mat
end

##############################################
#					SCRATCH			         #
##############################################

#= 
using DifferentialEquations
using Plots
using DataFrames
using Distributions
using Statistics
using Random
using StaticArrays
using BenchmarkTools
## scratch
y0 = R0 = 10
end_row = cumsum(y0)
#end_row # [10 25 50] 
#length(end_row)-1 # 3 -1 
#1:length(end_row)-1 # 1:2
#end_row[1:length(end_row)-1] # extract the elememts [1 3 6]
#1 .+ end_row[1:length(end_row)-1] # add 1 to each element 
#[1; 1 .+ end_row[1:length(end_row)-1]] # makes new array with [1 11 26]

state_par_match = SMatrix{1, 4}(1, 1, 1, 1) ## make static array?
state_geno_match = SMatrix{1, 4}(0, 0, 0, 0)
geno_par_match = SMatrix{1, 4}(0, 0, 0, 0)
which_param_quant = state_geno_match - geno_par_match
no_species = size(state_par_match, 1)
no_params = size(state_par_match, 2)
no_columns = no_params + 1 + size(state_geno_match, 2)  
init_comm_mat =  MMatrix{R0, no_columns}(fill(NaN, R0, no_columns))

n_par = size(state_par_match, 2) 
params = MVector{n_par}(fill(NaN, n_par))

params[1] = 4 ; # max birth
params[2] = 1; # min death
params[3] = 0.001; # density dependence of birth
params[4] = 0.001; # density dependence of death
params
cv_vect = [0.2, 0, 0, 0]
##

InitiatePop(R0, which_par_quant, state_geno_match, state_par_match,
	init_comm_mat, params, cv_vect)
=#
