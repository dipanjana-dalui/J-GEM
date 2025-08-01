##############################################
#		    FUNCTION Birth-Death             #
##############################################
#= Birth and death functions are going to be 
unique to your model, and you will need to edit 
the following chunks of code accordingly.

For the gillespie algo to induce demographic 
stochasticity, you will separate out the growths
of each state as birth and deaths. It is more of 
a book-keeping exercise - you will have to account 
for all terms that add to a state in its birth func
(intrinsic birth, nutrient/prey aquisition etc), 
and all the terms that subtract from state in the 
death func (like competition and predation). 

Symmetric ODE systems can be generalized using 
linear algebraic matrix notation, but care must
be taken regarding dimensionality of the matrix
operations. If your birth or death function is 
outputting an array corresponding to multiple states,
then be extra careful while using the output downstream.

=#
# need: b_max, b_s, N, d_min, d_s

function Birth(b_max,b_s, R)
    b_new = max((b_max - b_s*R[1])*R[1],0)
end

function Death(d_min, d_s, R)
    d_new = (d_min + d_s*R[1])*R[1]
end

function event_terms(params_next::Matrix{Float64}, N::Vector{Int64})
    b_max = params_next[1,1] # max birth
    d_min = params_next[1,2] # min death
    b_s = params_next[1,3] # density dependence of birth
    d_s = params_next[1,4]

    birth_H =  Birth(b_max, b_s, N)
    death_H =  Death(d_min, d_s, N)
    return birth_H, death_H
end

