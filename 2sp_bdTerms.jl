# make birth death function 
# birth and death for H1 
function birth_prey(b_max::Float64, N::Vector{Int64})
    birth = max((b_max[1]*N[1]),0)
end

function death_prey(d_min::Float64, N::Vector{Int64}, scr::Float64)
    death =  d_min[1]*N[1]
    pred_death = scr[1]*N[1]*N[2]
    pred_death = pred_death[1]
    death = death + pred_death
end

#=function carrying_K(b_max, d_min, b_s, d_s)
    K = floor(vec((b_max - d_min)/(b_s + d_s))[1])
    return K
end=#

# birth and death of pathogens
function birth_pred(scr::Float64, N::Vector{Int64}, fec::Float64)
    birth = fec[1] .* scr[1] .* N[1] .* N[2]
end

function death_pred(m::Float64, N::Vector{Int64})
    death = m[1]*N[2]
end


function event_terms(params_next::Matrix{Float64}, N::Vector{Int64})
    b_max = params_next[1,1] # max birth
    d_min = params_next[1,2] # min death
    #b_s = params_next[1,3] # density dependence of birth
    #d_s = params_next[1,4]
    scr = params_next[2,3]
    fec = params_next[2,4]
    m = params_next[2,5]

    birth_H =  birth_prey(b_max, N)
    death_H =  death_prey(d_min, N, scr)
    birth_P =  birth_pred(scr, N, fec)
    death_P =  death_pred(m, N) 
    return birth_H, death_H, birth_P, death_P
end
