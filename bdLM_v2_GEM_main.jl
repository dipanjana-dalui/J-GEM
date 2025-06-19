# birth-death logistic model GEM example

# all packages already loaded

# all AuxillaryFunc already loaded 


function GEM_sim(GEM_ver::Vector{Int64}, # Adjust type if known (e.g., Vector{String})
    t_max::Float64,
    min_time_step_to_store::Float64,
    no_species::Int64,
    num_rep::Int64,
    no_columns::Int64,
    no_params::Int64,
    R0::Vector{Float64}, # Adjust type if known (e.g., Vector{Int})
    which_par_quant::Matrix{Int64},
    state_geno_match::Matrix{Int64},
    state_par_match::Matrix{Int64},
    params::Vector{Float64},
    cv_vect::Vector{Float64},
    h2_vect::Vector{Float64},
    pop_stand_out_all::AbstractArray{Float64, 4} # Output array to be populated
    )   

    for j = 1:length(GEM_ver)
        println("\nj: $j")
    #    j = 1 # for debugging
        t = 0 # initiate
        stand_time = range(0, t_max, step = min_time_step_to_store)
        stand_time = collect(stand_time)
        num_time_steps = length(stand_time)

        ## memory preallocation should happen here for the entire simulations
        ## these are all the data structure we are storing
        pop_stand = zeros(no_species, num_time_steps, num_rep )
        #@show pop_stand
        x_stand = fill(NaN, no_columns-1,num_time_steps, no_species,num_rep) #trait
        x_var_stand = fill(NaN, no_columns-1,num_time_steps, no_species, num_rep) # trait variance 

        pop_data_out = fill(NaN, 3, num_time_steps, no_species)
        #@show pop_data_out
        pop_stand_out = fill(NaN, no_species, num_time_steps, num_rep)
        x_stand_out = fill(NaN, no_columns-1,num_time_steps, no_species,num_rep) #trait
        x_var_stand_out = fill(NaN, no_columns-1,num_time_steps, no_species, num_rep) # trait variance 

        for i = 1:num_rep
            println("\ni: $i")
    #       i=1
            t = 0
            #save a copy inside the loop to update
            R0 = copy(R0)
            R = copy(R0)
            println("\nR0: $R0")
            ## all structures made here are temp and exist only inside the
            ## replication loop
            #init_comm_mat =  MMatrix{R0, no_columns}(fill(NaN, comm_mat_rows, no_columns))
            init_comm_mat =  Array{Float64}(fill(NaN, Int(sum(R0)), no_columns))
            #@show init_comm_mat
            pop_slice = Array{Float64}(fill(0.0, no_species, num_time_steps))
            x_slice = fill(NaN, no_columns-1, num_time_steps, no_species)
            x_var_slice = fill(NaN, no_columns-1, num_time_steps, no_species)

            # store details of initial state 
            pop_slice[:,1] = R0 #first col/time step gets the initial pop 

            #draw initial population: 10x9 Matrix{Float64}
            x_dist_init = InitiatePop(R0, which_par_quant, state_geno_match, state_par_match,
                init_comm_mat, params, cv_vect) ## sum(R0) number of rows and no_col num of cols
            
            #@show x_dist_init
                
            for ii = 1:no_species #8×21×1 Array{Float64, 3}:[:, :, 1] 
                #ii = 1
                x_slice[:, 1, ii] = CalcAverageFreqs(ii, no_columns, no_params, x_dist_init)
                x_var_slice[1:no_params,1,ii] = mapslices(var ∘ skipmissing,x_dist_init[x_dist_init[:,1] .== ii,2:no_params+1],dims=1)

            end

            # count up each individual for all states after the first sampling
            for jj = 1:length(R)
                    x = init_comm_mat[:,1] #extract first col
                    #typeof(x)
                    R[jj] = count(.==(jj), x)
                    #R[jj] = count(x -> (x.==jj), init_comm_mat)
            end # this is suppose to give the count of individuals in each population
            #@show R

            x_dist = x_dist_init
            time_step_index = 2
            time_step = stand_time[time_step_index]


            while t < t_max  && 1 <=  sum(R)  
                
                # Find who is next
                FindWhoNext = WhoIsNext(x_dist,no_species,no_columns,no_params,R,
                                            state_par_match,state_geno_match)
                #@show FindWhoNext
                params_next = FindWhoNext[1]
                genotypes_next = FindWhoNext[2]
                whosnext = FindWhoNext[3]
                # R birth
                b_max = params_next[1] # max birth
                d_min = params_next[2] # min death
                b_s = params_next[3] # density dependence of birth
                d_s = params_next[4]
                #@show [b_max, d_min, b_s, d_s]
                
                # note[26.02.25] b_s and d_s should actually be calc using K (see DD derivation)
                # set b_s = b_max/K and d_s = d_min/K.
                birth_R =  max((b_max - b_s*R[1])*R[1],0)
                death_R =  (d_min + d_s*R[1])*R[1]
                #birth_R =  (b_max - (b_max/K)*R[1])*R[1]
                #death_R =  (d_min - (d_min/K)*R[1])*R[1]
                

                terms = [birth_R, death_R]
                #@show terms

                PickedEvent = PickEvent(terms, no_species)
                c_sum = PickedEvent[1]
                row = PickedEvent[2]
                col = PickedEvent[3]
                #@show PickedEvent
                ## row = 1 && col == 1 - sp1 birth
                ## row = 2 && col == 1 - sp1 death
                ## row = 1 && col == 2 - sp2 birth
                ## row = 2 && col == 2 - sp2 death 
                if row == 1 && col == 1
                    parent_traits = x_dist[Int(whosnext[col]), 2:no_columns] 
                            
                    new_trait = DrawNewTraits(x_dist,parent_traits,h2_vect,no_params,no_columns,col)
                    # new_trait = offspring_traits, offspring_genotypes
                    # ([4.676604793653037 1.0 0.0011999999999999995 0.0010000000000000002], [NaN NaN NaN NaN])
    
                    new_trait_row = hcat(col, new_trait)
                    #1  ([3.72285 1.0 0.0012 0.001], [NaN NaN NaN NaN])
                    new_trait_row = hcat(new_trait_row[1],new_trait_row[2][1],new_trait_row[2][2])

                    x_dist = vcat(x_dist, new_trait_row) 
                elseif row == 2 && col ==1 # death 
                    # delete the individual by making a new copy of the matrix
                    # without the row 
                    x_dist = x_dist[1:size(x_dist, 1) .!= Int(whosnext[col]), :]
                end

                ## UPDATE ABUNDANCES 
                for jj in 1:no_species
                    R[jj] = sum(x_dist[:,1].== jj)
                end

                if t > time_step
                    pop_slice[1:no_species,time_step_index] .= R # assign current values to sliced standard times
                    
                    for ii in 1:no_species
                        x_slice[:,time_step_index,ii] = CalcAverageFreqs(ii,no_columns,no_params,x_dist)
                        x_var_slice[1:no_params,time_step_index,ii] = mapslices(var∘ skipmissing, x_dist[x_dist[:,1].==ii,2:no_params+1],dims=1)
                    
                    end
                    time_step_index = time_step_index + 1 # advance to next standardized time
                    time_step = stand_time[time_step_index]
                    
                end

                #last thing to do before exiting the loop:
                # ADVANCE TIME 
                time_advance = exp(-1/c_sum[end])/(c_sum[end])
                if isnan(time_advance) == 0 && time_advance > 0 ##
                    t = t + time_advance
                else
                    break
                end
                #println("\nUpdated_R: $R")
                #println("\nt: $t")
            end #end of while loop running the core GEM

            # store the last value of the replicate
            pop_slice[1:no_species, time_step_index] = R
            for ii = 1:no_species
                x_slice[:,time_step_index,ii] = CalcAverageFreqs(ii,no_columns,no_params,x_dist);
                x_var_slice[1:no_params,time_step_index,ii] = mapslices(var∘ skipmissing, x_dist[x_dist[:,1].== ii,2:no_params+1],dims=1)
            end

            ## save the pop and param/geno values for the whole replicate in premade containers  c
            #i = 1
            pop_stand[:, :, i] = pop_slice
            x_stand[:, :, :, i] = x_slice 
            x_var_stand[:, :, :, i] = x_var_slice ##
            #@show pop_stand

        end ## end of for loop for single simulation replicate
        
        pop_stand_out_all[:, :, :, j] = pop_stand 
        x_stand_out_all[:, :, :, :, j] = x_stand
        x_var_stand_out_all[:, :, :, :, j] = x_var_stand
        #=
    This is where you will do any post simulation processing of the raw data
        =#
    end
    
    return nothing
end

out = vec(pop_stand_out_all);
out_mat = reshape(out, num_time_steps, num_rep, length(GEM_ver))
typeof(out_mat)
#=
out_dat_ver1 = DataFrame(out_mat[:,:,1], :auto)

dat = hcat(stand_time, out_dat_ver1, makeunique=true)

plot(dat.x1, dat.x1_1)
plot!(dat.x1, dat.x2)
plot!(dat.x1, dat.x3)
plot!(dat.x1, dat.x4)
plot!(dat.x1, dat.x5)
=#