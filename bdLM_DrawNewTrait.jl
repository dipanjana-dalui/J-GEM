##############################################
#		  FUNCTION DRAW NEW TRAITS           #
##############################################

function DrawNewTraits(x_dist, parent_traits, h2_vect, no_params, no_columns, col)
	
	# QUANTITATIVE TRAITS
	pop_mean = mapslices(mean ∘ skipmissing, x_dist[x_dist[:, 1] .== col, 2:no_params+1], dims=1)
	# we want the means to be for each parameter
	pop_var = mapslices(var ∘ skipmissing, x_dist[x_dist[:, 1] .== col, 2:no_params+1], dims=1)
	pop_stdev = sqrt.(pop_var)

	## CHECK WITH h2 structure, may not need col (state ID)
	exp_offspring_traits = h2_vect[col] .* reshape(parent_traits[1:no_params], 1, 4) .+ (1-h2_vect[col]) .* pop_mean
	# 1×4 Matrix{Float64}:
 	# 3.7684  1.0  0.0012  0.001
	sigma_offspring_traits = sqrt(1-(h2_vect[col])^2)*pop_stdev
	# 1×4 Matrix{Float64}:
	# 0.836441  0.0  0.0  0.0
	MU = log.(exp_offspring_traits.^2 ./ sqrt.(sigma_offspring_traits.^2 .+ exp_offspring_traits.^2)) # mean for lognormal
	# 1×4 Matrix{Float64}:
	# 1.3026  0.0  -6.72543  -6.90776
	
	SIGMA = sqrt.(log.(sigma_offspring_traits.^2 ./ exp_offspring_traits.^2 .+ 1)) # std for lognormal
	# 1×4 Matrix{Float64}:
 	# 0.219299  0.0  0.0  0.0
	
	####### 
	offspring_traits = Matrix{Float64}(undef, 1, length(SIGMA))
	for i in 1:length(SIGMA)
		offspring_traits[1,i] = rand(LogNormal(MU[i], SIGMA[i])) # pull traits
	end
	# offspring_traits
	# 1×4 Matrix{Float64}:
 	# 3.72261  1.0  0.0012  0.001

	offspring_genotypes = parent_traits[no_params+1:no_columns-1]
	offspring_genotypes = reshape(offspring_genotypes, 1, length(SIGMA))
	
	## CHECK WITH John
	#traits_out = [offspring_traits offspring_genotypes]
	## this does hcat
	## for vcat, do [offspring_traits ; offspring_genotypes]
	return offspring_traits, offspring_genotypes
	
end


##############################################
#					SCRATCH			         #
##############################################
#=

R0 = 10
y0 = R0 = 10
state_par_match = SMatrix{1, 4}(1, 1, 1, 1) ## make static array?
state_geno_match = SMatrix{1, 4}(0, 0, 0, 0)
geno_par_match = SMatrix{1, 4}(0, 0, 0, 0)
which_param_quant = state_geno_match - geno_par_match
no_species = size(state_par_match, 1)
no_params = size(state_par_match, 2)
no_columns = no_params + 1 + size(state_geno_match, 2)

x_dist
col = 1 ## confusing - is this the index or the value??
# col right now is 1 and acts as the species id
parent_traits = x_dist[Int(whosnext[1]), 2:no_columns]
no_columns
no_params

=#