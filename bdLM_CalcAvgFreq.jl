##############################################
#		FUNCTION CALC AVG FREQUENCIES        #
##############################################

function CalcAveragesFreqs(ii, no_columns, no_params, x_dist)
	# x_dist = x_dist_init
	# 2:no_params+1
	# ii=1
	# x_dist[x_dist[:, 1] .== ii, 2:no_params+1]
	# mapslices(median ∘ skipmissing, x_dist[x_dist[:, 1] .== ii, 2:no_params+1], dims=1)
	# Step 1: Quantitative trait medians (ignoring NaN values)
	# takes the mean of the parameter columns 
	qt_means = mapslices(median ∘ skipmissing, x_dist[x_dist[:, 1] .== ii, 2:no_params+1], dims=1)
	## if this part breaks, change the dim and check. Dims = 1 gives you colsums

	# Step 2: Discrete trait frequencies
	# take the sum of the genotype cols 
	gt_freqs = sum(x_dist[x_dist[:, 1] .== ii,2+no_params:no_columns], dims=1) ./ sum(x_dist[:, 1] .== ii)
	# Step 3: Combine the results
	means_freqs = hcat(qt_means, gt_freqs)
	return means_freqs

end


#=
ii = 1
no_params = size(state_par_match, 2)
no_columns = no_params + 1 + size(state_geno_match, 2) 
x_dist_init
=#

# CalcAveragesFreqs(ii, no_columns, no_params, x_dist_init)