##############################################
#		FUNCTION CALC AVG FREQUENCIES        #
##############################################
function CalcMedian(ii::Int64, no_columns::Int64, no_param::Int64, x_dist::Matrix{Float64})
	# Step 1: Quantitative trait medians (ignoring NaN values)
	qt_medians = mapslices(median ∘ skipmissing, x_dist[x_dist[:, 1] .== ii, 2:no_param+1], dims=1)
	# Step 2: Discrete trait frequencies
	# take the sum of the genotype cols 
	gt_freqs = sum(x_dist[x_dist[:, 1] .== ii,2+no_param:no_columns], dims=1) ./ sum(x_dist[:, 1] .== ii)
	# Step 3: Combine the results
	median_freqs = hcat(qt_medians, gt_freqs)
	return median_freqs
end

function CalcMean(ii::Int64, no_columns::Int64, no_param::Int64, x_dist::Matrix{Float64})
	# Step 1: Quantitative trait medians (ignoring NaN values)
	qt_means = mapslices(mean ∘ skipmissing, x_dist[x_dist[:, 1] .== ii, 2:no_param+1], dims=1)
	# Step 2: Discrete trait frequencies
	# take the sum of the genotype cols 
	gt_freqs = sum(x_dist[x_dist[:, 1] .== ii,2+no_param:no_columns], dims=1) ./ sum(x_dist[:, 1] .== ii)
	# Step 3: Combine the results
	mean_freqs = hcat(qt_means, gt_freqs)
	return mean_freqs
end

function CalcVar(ii::Int64, no_param::Int64, x_dist::Matrix{Float64})
	temp_size = size(x_dist[x_dist[:,1] .== ii,2:no_param+1])
	if temp_size[1] == 1 #has only one row = only one individual
		trait_var = Matrix{Float64}(fill(0, temp_size[1], length(2:no_param+1)))
	else
		trait_var = mapslices(var ∘ skipmissing, x_dist[x_dist[:,1] .== ii,2:no_param+1],dims=1)
	end
	return trait_var
end
