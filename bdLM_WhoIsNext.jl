## scratch
#
using Distributions
using Random
using StaticArrays


R0 = 10
y0 = R0 = 10
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
params[4] = 0.001;
params
cv_vect = [0.2, 0, 0, 0]


x_dist_init = InitiatePop(y0, which_param_quant, state_geno_match, state_par_match,
init_comm_mat, params, cv_vect)
x_dist = x_dist_init

#WhoIsNext(x_dist, no_species, no_columns, no_params, y0, state_par_match, state_geno_match)

##############################################
#				 FUNCTION
##############################################
#function WhoIsNext(x_dist, no_species, no_columns, no_params, y0,
#	state_par_match, state_geno_match)

	params_next = zeros(no_species, no_params)
	genotype_next = zeros(no_species, size(state_geno_match, 2))
	whosnext = fill(NaN,length(R0)) 
	for zz = 1:no_species #loop through all states that are present in ODE form
		zz = 1
		ind_in_state = findall(x_dist[:,1] .== zz) # this finds the index for the zzth state
		# this gives the range of index values b/c we are searching for the index all the
		# individuals with species index zz.
		which_params = findall(state_par_match[1, :] .!= 0) # finding the indices of the non zero elements
		which_genotype = 1:size(state_geno_match, 2) #1:ncol of genotypes (NOT just the non zero) -- WHY?
		#while !isempty(ind_in_state)
		if !isempty(ind_in_state)
			which_row = rand(1:length(ind_in_state)) # generate a random interger between 1 and the index
			whosnext[zz] = ind_in_state[which_row]
			params_next[zz, which_params] = x_dist[Int(whosnext[zz]), 1 .+ which_params]
			genotype_next[zz,which_genotype] = x_dist[Int(whosnext[zz]), 2 .+ no_params:no_columns]
		end
	end
	return params_next, genotypes_next, whosnext
end
