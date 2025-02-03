##############################################
#		  FUNCTION MEDIANS AND CI            #
##############################################
sequence = x_stand

dimensions = size(sequence)
dimensions[4]

dim 1 is (no_column-1) - number of parameters/genotype
dim 2 is number of time steps 
dim 3 is number of species
dim 4 is number of replicates

we want to make a container with following Dims
new_dim 1 is 95% CI and median (=3)
new_dim 2 is number of time steps 
new_dim 3 is number of parameters (no_columns-1)
new_dim 4 is number of species 

function MedianCI(upper_ci_level, lower_ci_level, sequence)
    dimensions = size(sequence) #8x21x1x5
    A = fill(NaN, 3, dimensions[2], dimensions[1], dimensions[3])
    # 3x21x8x1 

    for ii = 1:dimensions[3] #loop through species 
        for i = 1:dimensions[1] #loop through parameters
            A[1, :, i, ii] =  percentile.(eachcol(sequence[1, :, 1, :]), 25)
            
            A[1, :, i, ii] =  percentile.(eachcol(mymatrix), lower_ci_level)
            A[2, :, i, ii] =  percentile.(eachcol(mymatrix), 50)
            A[3, :, i, ii] =  percentile.(eachcol(mymatrix), upper_ci_level)
