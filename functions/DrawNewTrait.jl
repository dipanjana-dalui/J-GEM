##############################################
#		  FUNCTION DRAW NEW TRAITS           #
##############################################
function DrawNewTraits(x_dist::Matrix{Float64}, parent_traits::Vector{Float64}, 
						h2_vect::Matrix{Float64}, no_param::Int64, no_columns::Int64, col::Int64, j::Int64)

	pop_size = size(x_dist[x_dist[:, 1] .== col, 2:no_param+1])
	# QUANTITATIVE TRAITS
	pop_mean = mapslices(mean ∘ skipmissing, x_dist[x_dist[:, 1] .== col, 2:no_param+1], dims=1)

	if pop_size[1] == 1
		pop_var = Matrix{Float64}(fill(0, pop_size[1], pop_size[2]))
	else
		pop_var = mapslices(var ∘ skipmissing, x_dist[x_dist[:, 1] .== col, 2:no_param+1], dims=1)
	end
	pop_var = round.(pop_var, digits=6)
	pop_stdev = round.(sqrt.(pop_var), digits=6)

	exp_offspring_traits = h2_vect[j,col] .* reshape(parent_traits[1:no_param], size(pop_mean)[1], size(pop_mean)[2]) .+ (1-h2_vect[j,col]) .* pop_mean		
	sigma_offspring_traits = sqrt(1-(h2_vect[j,col])^2)*pop_stdev
	# MU and SIGMA for lognormal
	MU = log.(exp_offspring_traits.^2 ./ sqrt.(sigma_offspring_traits.^2 .+ exp_offspring_traits.^2)) 	
	SIGMA = sqrt.(log.(sigma_offspring_traits.^2 ./ exp_offspring_traits.^2 .+ 1))

	offspring_traits = fill(NaN, 1, length(exp_offspring_traits))
	for i = 1:length(exp_offspring_traits)
		if !isnan(MU[i]) && !isnan(SIGMA[i])
			if !iszero(SIGMA[i])
				offspring_traits[1,i] = rand(LogNormal(MU[i], SIGMA[i])) # pull traits when sigma != 0
			else
				offspring_traits[1,i] = exp_offspring_traits[i]
			end 
		end
	end
	
	offspring_genotypes = collect(transpose(parent_traits[no_param+1:no_columns-1]))
	##### what happens to genotype?

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