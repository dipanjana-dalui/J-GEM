##############################################
#		  FUNCTION PICK individuals          #
##############################################
## This function picks individuals with mean parameter, and std parameter * cv 

function PickIndiv(x, std, N)
	# this picks N traits with mean x and std deviation stand
	MU = log(x .^2 ./ sqrt(std .^2 + x .^2))
	SIGMA = sqrt(log(std .^2 ./x .^2 + 1))
	sample = LogNormal(MU, SIGMA)
	x_out = rand(sample, N)
	return x_out
end
