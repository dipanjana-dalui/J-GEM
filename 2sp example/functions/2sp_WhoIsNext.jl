
##############################################
#			FUNCTION WHO IS NEXT?	         #
##############################################
function WhoIsNext(x_dist::Matrix{Float64}, no_species::Int64, no_columns::Int64, no_param::Int64, 
	N0::Vector{Int64}, state_par_match::Matrix{Int64}, state_geno_match::Matrix{Int64})

	state_par_match = zero_to_nan(state_par_match)
	state_geno_match = zero_to_nan(state_geno_match)
	
	params_next = fill(NaN,no_species, no_param)
	genotype_next = fill(NaN, no_species, size(state_geno_match, 2))
	whosnext = fill(NaN,length(N0)) 

	for zz = 1:no_species
		ind_in_state = findall(x_dist[:,1] .== zz) # this finds the index for the zzth state
		which_params = findall(.!isnan.(state_par_match[zz, :])) # finding the indices of the non NaN elements	
		which_genotype = 1:size(state_geno_match, 2) #1:ncol of genotypes (NOT just the non zero) -- WHY?
		#while !isempty(ind_in_state)
		if !isempty(ind_in_state)
			which_row = rand(1:length(ind_in_state)) 
			whosnext[zz] = ind_in_state[which_row]
			params_next[zz, which_params] = x_dist[Int(whosnext[zz]), 1 .+ which_params]
			genotype_next[zz,which_genotype] = x_dist[Int(whosnext[zz]), 2 .+ no_param:no_columns]
		#else
			#params_next[zz, which_params] = 0 #params_next[zz, which_params] 
			#genotypes_next[zz,which_genotypes] = 0 # if pop is gone, set genotypes to 0
		end
	end
	return params_next, genotype_next, whosnext
end
