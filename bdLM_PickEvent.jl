R = [10]
no_species = length(R)
terms = [20; 30]
length(terms)
reshape(terms, 1,length(terms))


function PickEvent(terms, no_species)
	terms = reshape(terms, 1, length(terms)) #reshape(ele_to_reshape, new_row, new_col)
	c_sum = cumsum(terms, dims = 2)  
	pie_slices = c_sum ./ c_sum[end] #generated weighted slices b/w 0-1
	#r_num = rand()
	r_num = 0.85
	
	less_than = r_num .< pie_slices
	typeof(less_than)
	
	event_index = findfirst(less_than)
	typeof(event_index)

	CS_elements = reshape(1:length(terms), :, no_species)
    
	#row, col = ind2sub(size(CS_elements), event_index)
	CartesianIndex(CS_elements)
## ASK John
# we can just assign the row and col from the cartesian coordinates, right?
	return event_index, c_sum, row, col 

end
#=
#function PickEvent(terms, no_species)
	terms = reshape(terms, 1, length(terms)) #make everything one row
	c_sum = cumsum(terms, dims = 2)   
	pie_slices = c_sum ./ c_sum[end] #generated weighted slices b/w 0-1
	r_num = rand()
	less_than = r_num .< pie_slices
	event_index = findfirst(less_than)

	# Create CS_elements and reshape it
    CS_elements = 1:length(terms)
    CS_elements = reshape(CS_elements, :, no_species)
	index = findfirst(CS_elements .== event_index)
	row, col = Tuple(findall(CS_elements .== event_index)[1])
	########
	
	terms = reshape(terms, 1, length(terms))
	c_sum = cumsum(terms, dims = 2)   
	
	pie_slices = c_sum ./ c_sum[end]
	r_num = rand()
	less_than = r_num .< pie_slices
	
	# Find the event index (first index where less_than is true)
	event_index = findfirst(less_than)
	# Create CS_elements and reshape it
	CS_elements = 1:length(terms)
	CS_elements = reshape(CS_elements, :, no_species)
	
	# Find the linear index where CS_elements == event_index
	linear_index = findfirst(CS_elements .== event_index)
		
	# Convert the linear index to row and column
	row, col = ind2sub(size(CS_elements), linear_index)
	
	return event_index, c_sum, row, col
	end
	


	##
	
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
	