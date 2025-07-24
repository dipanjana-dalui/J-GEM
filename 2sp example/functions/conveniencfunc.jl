## Function to convert 0s to NaN in matrix
function zero_to_nan(Bool_mat::Matrix{Int64})
	return [x==0 ? NaN : x for x in Bool_mat]
end