function DrawNewTraits(x_dist, parent_traits, h2, no_params, no_columns, col)
	
	# QUANTITATIVE TRAITS
	
	pop_mean = mean(x_dist[x_dist[:,1]==col, 2:no_params+1])
	pop_stdev = std(x_dist[x_dist[:,1]==col, 2:no_params+1])
	exp_offspring_traits = h2[col] .* parent_traits[1:no_params] .+ [1-h2[col]] .* pop_mean
	sigma_offspring_traits = sqrt(1-(h2[col])^2)*pop_stdev
	MU = log(exp_offspring_traits .^2 ./ sqrt(sigma_offspring _traits . ^2+exp_offspring _traits . ^2)); # mean for lognormal
	SIGMA = sqrt(log(sigma_offspring _traits .^2 ./exp_offspring_traits .^2 + 1)); # std for lognormal
	offspring_traits = lognrnd(MU,SIGMA,1,length(SIGMA)); # pull traits
	offspring_genotypes = parent_traits[no_params+1:no_columns-1]
	traits_out = [offspring_traits offspring_genotypes] ## this does hcat
	## for vcat, do [offspring_traits ; offspring_genotypes]
	return traits_out
	
end
