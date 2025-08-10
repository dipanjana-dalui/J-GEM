##############################################
#		  	  FUNCTION PICK EVENTS           #
##############################################

function PickEvent(terms::Vector{Float64}, no_species::Int64)
	terms = reshape(terms, 1, length(terms)) #reshape(ele_to_reshape, new_row, new_col)
	c_sum = cumsum(terms, dims=2)  
	pie_slices = c_sum ./ c_sum[end] #generated weighted slices b/w 0-1
	r_num = rand()
	
	less_than = r_num .< pie_slices  #BitMatrix
	less_than = collect(less_than)  #Matrix{Bool}
	event_mat = reshape(less_than, :, no_species) #reshaped Matrix{Bool}
	
	row = -1
	col = -1
	for r = 1:size(event_mat,1) #want to go by row
		col_ind_in_row_r = findfirst(==(1),event_mat[r,:])
		if col_ind_in_row_r != nothing
			row = r
			col = col_ind_in_row_r
			break
		end
	end

	return c_sum, row, col 
	

end
