##############################################
#		FUNCTION INITIATE POPULATION         #
##############################################
function InitiatePop(N0::Vector{Int64}, which_par_quant::Matrix{Int64}, state_geno_match::Matrix{Int64}, 
	state_par_match::Matrix{Int64}, init_comm_mat::Matrix{Float64}, params::Vector{Float64}, 
	cv_vect::Matrix{Float64}, j::Int64)
	
	# some conversions
    state_par_match = zero_to_nan(state_par_match)
    state_geno_match = zero_to_nan(state_geno_match)
    #geno_par_match = zero_to_nan(geno_par_match)
	which_par_quant = zero_to_nan(which_par_quant)

	y0 = N0
	end_row = cumsum(y0)
	starting_row = [1; 1 .+ end_row[1:length(end_row)-1]]
	n_sp = length(y0)
	traits_to_assign = which_par_quant .* y0 # it accounts for the no. indiv / spp
	params_to_pick = collect(transpose(state_par_match .* repeat(params', n_sp, 1)))

	gts_to_assign = state_geno_match .* y0
	num_gts = size(state_geno_match, 2)
	for qq = 1:n_sp
		init_comm_mat[Int(starting_row[qq]):Int(end_row[qq]), 1] .= qq
		for zz = 1:length(params)
			if !isnan(params_to_pick[zz,qq]) 
				temp = PickIndiv(params_to_pick[zz,qq],cv_vect[j, qq]*params_to_pick[zz, qq],Int(traits_to_assign[qq,zz])) 
			else 
				#temp = zeros(Int(y0[qq]),1)
				temp = fill(NaN, Int(y0[qq]),1)
			end 
			init_comm_mat[Int(starting_row[qq]):Int(end_row[qq]), 1+zz] .= temp
		end
	
		# ALL?
		#=if all(gts_to_assign[qq,:] .> 0)
			 	# makes a matrix with communuty size no. of rows, and number of geno cols
			genotype = zeros(Int(y0[1]), num_gts) 
			# go row by row and randomly choose genotype to set to 1
			for yy = 1:y0[1]
				yy = 1
				temp2 = rand(1:num_gts) # try this alt temp = rand(1:num_gts, 1)
				genotype[yy, temp2] = 1
			end			
		end
		init_comm_mat[Int(starting_row[qq]):Int(end_row[qq]), 2 + length(params):1+length(params)+num_gts] = genotype
		=#
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
