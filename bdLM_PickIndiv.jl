#using Distributions

function PickIndiv(x, std, N)
	# this picks N traits with mean x and std deviation stand
	MU = log(x .^2 ./ sqrt(std .^2 + x .^2))
	SIGMA = sqrt(log(std .^2 ./x .^2 + 1))
	sample = LogNormal(MU, SIGMA)
	x_out = rand(sample, N)
	#x_out = lognormal(MU, SIGMA, N, length(SIGMA)) # matlab ver
	# length(SIGMA) -- to loop over various values of sigma??
	return x_out
end

#=
using Distributions
using Plots
check
x = 1
std = 0.1
N = 10
L=PickIndiv(x, std, N)
histogram(L)
=#
