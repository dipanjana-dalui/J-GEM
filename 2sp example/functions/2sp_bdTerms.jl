# make birth death function 
# birth and death for H1 
function birth_host(b_max::Float64, b_s::Float64, N::Vector{Int64})
    birth = max((b_max[1] - b_s[1]*N[1])*N[1],0)
end

function death_host(d_min::Float64, d_s::Float64, N::Vector{Int64}, sp::Float64)
    comp_death =  (d_min[1] + d_s[1]*N[1])*N[1]
    infect_death = sp[1]*N[1]*N[2]
    infect_death = infect_death[1]
    death = comp_death + infect_death
end

#=function carrying_K(b_max, d_min, b_s, d_s)
    K = floor(vec((b_max - d_min)/(b_s + d_s))[1])
    return K
end=#

# birth and death of pathogens
function birth_path(sp::Float64, N::Vector{Int64}, fec::Float64)
    birth = fec[1] .* sp[1] .* N[1] .* N[2]
end

function death_path(m::Float64, N::Vector{Int64})
    death = m[1]*N[2]
end


function event_terms(params_next::Matrix{Float64}, N::Vector{Int64})
    b_max = params_next[1,1] # max birth
    d_min = params_next[1,2] # min death
    b_s = params_next[1,3] # density dependence of birth
    d_s = params_next[1,4]
    sp = params_next[2,5]
    fec = params_next[2,6]
    m = params_next[2,7]

    birth_H =  birth_host(b_max, b_s, N)
    death_H =  death_host(d_min, d_s, N, sp)
    birth_P =  birth_path(sp, N, fec)
    death_P =  death_path(m, N) 
    return birth_H, death_H, birth_P, death_P
end