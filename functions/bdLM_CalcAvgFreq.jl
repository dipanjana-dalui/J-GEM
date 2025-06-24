##############################################
#		FUNCTION CALC AVG FREQUENCIES        #
##############################################
function CalcAverageFreqs(ii, no_columns, no_params, x_dist)
	# Step 1: Quantitative trait medians (ignoring NaN values)
	# takes the mean of the parameter columns 
	qt_means = mapslices(median ∘ skipmissing, x_dist[x_dist[:, 1] .== ii, 2:no_params+1], dims=1)
	
	# Step 2: Discrete trait frequencies
	# take the sum of the genotype cols 
	gt_freqs = sum(x_dist[x_dist[:, 1] .== ii,2+no_params:no_columns], dims=1) ./ sum(x_dist[:, 1] .== ii)
	# Step 3: Combine the results
	mean_freqs = hcat(qt_means, gt_freqs)
	return mean_freqs

end

##############################################
#					SCRATCH			         #
##############################################
#=
ii = 1 #number of species
no_columns = 9
no_params = 4
x_dist = x_dist_init

CalcAverageFreqs(ii, no_columns, no_params, x_dist_init)

=#
