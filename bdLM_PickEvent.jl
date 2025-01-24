R = [10]
no_species = length(R)
terms = [10.0; 20.0; 30.0]
length(terms)
reshape(terms, 1,length(terms))


function PickEvent(terms, no_species)
	terms = reshape(terms, 1, length(terms)) #reshape(ele_to_reshape, new_row, new_col)
	c_sum = cumsum(terms, dims=2)  
	pie_slices = c_sum ./ c_sum[end] #generated weighted slices b/w 0-1
	#r_num = rand()
	r_num = 0.85
	
	less_than = r_num .< pie_slices
	less_than = collect(less_than)
	
	
	event_index = findfirst(less_than)
	typeof(event_index)
	Tuple(event_index)

	CS_elements = reshape(1:length(terms), :, no_species)
    a = [0, 1]
	CS_elements .== a
	tuple(findall(CS_elements .== a)[1])
	tuple(findall(CS_elements .== event_index)[1])
	
	return event_index, c_sum, row, col 

end

	#=
	
    # Generate random number and compare to pie slices
    r_num = rand()
    less_than = r_num .< pie_slices
    
    # Find the event index (first index where less_than is true)
    event_index = findfirst(less_than)
    
    # Create CS_elements and reshape it
    CS_elements = 1:length(terms)
    CS_elements = reshape(CS_elements, :, no_species)
    
    # Find the row and column corresponding to the event_index
    row, col = tuple(findall(CS_elements .== event_index)[1])
	
	
	return #will depend on what we want out of this code
end


#=
# Helper function to convert a linear index to Cartesian (row, col)
function ind2sub(dims, index)
		nrows, ncols = dims
		row = div(index - 1, ncols) + 1
		col = (index - 1) % ncols + 1
		return row, col
end
	