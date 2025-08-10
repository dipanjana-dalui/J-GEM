## Function to convert 0s to NaN in matrix
function zero_to_nan(Bool_mat::Matrix{Int64})
	return [x==0 ? NaN : x for x in Bool_mat]
end

## function to set a threshold for zero
function deci_threshold(df::DataFrame)
	#trait = trait_mean_df.b_max
	#round.(trait_mean_df.b_max, digits=4)

	df_mod = mapcols(col -> round.(col, digits=6), df)
	return df_mod
end
