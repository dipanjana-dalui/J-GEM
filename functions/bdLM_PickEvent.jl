##############################################
#		  	  FUNCTION PICK EVENTS           #
##############################################

function PickEvent(terms, no_species)
	terms = reshape(terms, 1, length(terms)) #reshape(ele_to_reshape, new_row, new_col)
	c_sum = cumsum(terms, dims=2)  
	pie_slices = c_sum ./ c_sum[end] #generated weighted slices b/w 0-1
	r_num = rand()
	
	less_than = r_num .< pie_slices  #BitMatrix
	less_than = collect(less_than)  #Matrix{Bool}
	event_mat = reshape(less_than, :, no_species) #reshaped Matrix{Bool}
	index = findfirst(==(1), event_mat) #CartesianIndex
	#event_index = Tuple(findfirst(less_than))
	#event_id = event_index[2]

	#bd_event_picked = terms[event_id]
	#bd_event_mat = reshape(terms, :, no_species)
	#index = findfirst(==(bd_event_picked), bd_event_mat)
	row = index[1]
    col = index[2]
	return c_sum, row, col 
	#return event_index, c_sum, row, col 

end

##############################################
#					SCRATCH			         #
##############################################
#=
R = [10 11]
no_species = length(R)
terms = [10.0; 10.0; 10; 10]
length(terms)
reshape(terms, 1,length(terms))
=#