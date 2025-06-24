#=
function AppendStates(num_similar_states, mat)
	temp = []
	for i in 1:size(mat, 1)
		push!(temp, repeat(mat[i, :], num_similar_states[i], 1))
		end
	return vcat(temp...)
end
=#

#NOT USING THIS FUNCTION THIS TIME