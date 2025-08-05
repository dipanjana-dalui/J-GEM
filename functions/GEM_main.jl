# birth-death logistic model GEM example

# all packages already loaded

# all AuxillaryFunc already loaded 


function GEM_sim(GEM_ver::Vector{String},
    t_max::Float64,
    no_species::Int64,
    num_rep::Int64,
    no_columns::Int64,
    no_param::Int64,
    N0::Vector{Int64}, 
    which_par_quant::Matrix{Int64},
    state_geno_match::Matrix{Int64},
    state_par_match::Matrix{Int64},
    param_init::Vector{Float64},
    cv_vect::Matrix{Float64},
    h2_vect::Matrix{Float64},
    par_names::Vector{String},
    geno_names::Vector{String},
    pop_stand_out_all::Array{Float64},
    x_stand_out_all::Array{Float64},
    x_var_stand_out_all::Array{Float64}
    )   

    # @assert statements should go here



    for j = 1:length(GEM_ver)
        println("\nj: $j")
        t = 0 # initiate
        pop_stand = zeros(no_species, num_time_steps, num_rep ) # pop abundances
        x_stand = fill(NaN, no_columns-1,num_time_steps, no_species,num_rep) #trait
        x_var_stand = fill(NaN, no_columns-1,num_time_steps, no_species, num_rep) # trait variance 

        # start replication loop 
        for i = 1:num_rep
            println("\ni: $i")            
            t = 0
            #save a copy inside the loop to update
            N0 = copy(N0)
            N = copy(N0) 
            params = copy(param_init) 
            
            println("\nN0: $N0")
            # @show params        
            ## all structures made here are temp and exist only inside the
            ## replication loop
            
            init_comm_mat =  Array{Float64}(fill(NaN, Int(sum(N0)), no_columns))
            pop_slice = Array{Int64}(fill(0, no_species, num_time_steps))
            x_slice = fill(NaN, no_columns-1, num_time_steps, no_species)
            x_var_slice = fill(NaN, no_columns-1, num_time_steps, no_species)

            # store details of initial state 
            pop_slice[:,1] = N0 #first col/time step gets the initial pop 

            #draw initial population: 10x9 Matrix{Float64}
            x_dist_init = InitiatePop(N0, which_par_quant, state_geno_match, state_par_match, init_comm_mat, params, cv_vect, j) 
            #@show x_dist_init
                
            for ii = 1:no_species
                x_slice[:, 1, ii] = CalcMedian(ii, no_columns, no_param, x_dist_init)
                #x_slice[:, 1, ii] = CalcMean(ii, no_columns, no_param, x_dist_init)
                x_var_slice[1:no_param,1,ii] = CalcVar(ii, no_param, x_dist_init)
            end

            # count up each individual for all states after the first sampling
            for jj = 1:length(N)
                x = init_comm_mat[:,1] #extract first col
                N[jj] = count(.==(jj), x)
            end

            x_dist = x_dist_init
            time_step_index = 2
            time_step = stand_time[time_step_index]


            while t < t_max  && 1 <=  sum(N)  
                
                # Find who is next
                FindWhoNext = WhoIsNext(x_dist,no_species,no_columns,no_param,N0,state_par_match,state_geno_match)
                #@show FindWhoNext
                params_next = FindWhoNext[1]
                genotypes_next = FindWhoNext[2]
                whosnext = FindWhoNext[3]

                terms = collect(event_terms(params_next, N))
                
                PickedEvent = PickEvent(terms, no_species)
                c_sum = PickedEvent[1]
                row = PickedEvent[2]
                col = PickedEvent[3]
                #@show PickedEvent
                ## row = 1 && col == 1 - sp1 birth
                ## row = 2 && col == 1 - sp1 death
                ## row = 1 && col == 2 - sp2 birth
                ## row = 2 && col == 2 - sp2 death 
                if row == 1 #&& col == 1
                    parent_traits = x_dist[Int(whosnext[col]), 2:no_columns] 
                    #@show parent_traits        
                    new_trait = DrawNewTraits(x_dist,parent_traits,h2_vect,no_param,no_columns,col, j)
                    
                    new_trait_row = hcat(col, new_trait)
                    new_trait_row = hcat(new_trait_row[1],new_trait_row[2][1],new_trait_row[2][2])
                    x_dist = vcat(x_dist, new_trait_row) 

                elseif row == 2 #&& col ==1 # death 
                    # delete the individual by making a new copy of the matrix
                    # without the row 
                    x_dist = x_dist[1:size(x_dist, 1) .!= Int(whosnext[col]), :]
                end

                ## UPDATE ABUNDANCES 
                for jj in 1:no_species
                    N[jj] = sum(x_dist[:,1].== jj)
                end
                #println("updated \nN: $N")

                if t > time_step
                    pop_slice[1:no_species,time_step_index] .= N  # assign current values to sliced standard times
                    
                    for ii in 1:no_species
                        
                        x_slice[:,time_step_index,ii] = CalcMedian(ii,no_columns,no_param,x_dist)
                        ## needs editing
                        x_var_slice[1:no_param,time_step_index,ii] = CalcVar(ii, no_param, x_dist)
                    end
                    time_step_index = time_step_index + 1 # advance to next standardized time
                    time_step = stand_time[time_step_index]
                    
                end

                #last thing to do before exiting the loop:
                # ADVANCE TIME 
                time_advance = exp(-1/c_sum[end])/(c_sum[end])
                if !isnan(time_advance) && time_advance > 0 ##
                    t = t + time_advance
                else
                    break
                    println("Time advance error. Stopped at time:\nT $T")
                end
                #println("\nUpdated_R: $R")
                #println("\nt: $t")
            end #end of while loop running the core GEM

            # store the last value of the replicate
            pop_slice[1:no_species, time_step_index] = N
            for ii = 1:no_species
                x_slice[:,time_step_index,ii] = CalcMedian(ii,no_columns,no_param,x_dist);
                x_var_slice[1:no_param,time_step_index,ii] = CalcVar(ii, no_param,x_dist)
            end

        
            pop_stand[:, :, i] = pop_slice
            x_stand[:, :, :, i] = x_slice 
            x_var_stand[:, :, :, i] = x_var_slice ##
            #@show pop_stand

        end ## end of for loop for single simulation replicate
        
        pop_stand_out_all[:, :, :, j] = pop_stand 
        x_stand_out_all[:, :, :, :, j] = x_stand
        x_var_stand_out_all[:, :, :, :, j] = x_var_stand
        
    end # END of GEM_ver loop

## SAVING THE POPULATION TIME SERIES IN DATAFRAME
    pop_out = vec(pop_stand_out_all)
    pop_out_mat = reshape(pop_out, no_species, num_time_steps, num_rep, length(GEM_ver))

    pop_out_gem_v_store = Vector{DataFrame}(undef, length(GEM_ver)) 
    pop_out_spp_store = Vector{DataFrame}(undef, no_species)

    for k = 1:no_species
        for i = 1:length(GEM_ver)
            temp = DataFrame(hcat(stand_time, fill(i,num_time_steps),
                            fill(k,num_time_steps)),:auto)
            rename!(temp, :x1 => :time, :x2 => :GEM_ver, :x3 => :state_ID)
            pop_out_temp = DataFrame(pop_out_mat[k,:,:,i], :auto)
            pop_out_temp = hcat(temp, pop_out_temp, makeunique=true)
            pop_out_gem_v_store[i] = pop_out_temp
        end
        pop_out_spp_store[k] = vcat(pop_out_gem_v_store...)
    end
    pop_df = vcat(pop_out_spp_store...)

    ## SAVING THE PARAMETER AND GENOTYPE MEANS

    x_out = vec(x_stand_out_all)
    x_out_mat = reshape(x_out, no_columns-1,num_time_steps, no_species,num_rep, length(GEM_ver))
    
    ### these names should be imported from the setup file
    col_names = vcat(par_names, geno_names)

    x_out_spp_store = Vector{DataFrame}(undef, no_species)
    x_out_gem_v_store = Vector{DataFrame}(undef, length(GEM_ver)) 
    x_out_rep_store = Vector{DataFrame}(undef, num_rep)

    for k = 1:no_species
        #k = 2
        for i = 1:length(GEM_ver)
            #i = 1
            for j = 1:num_rep
                col_df = DataFrame(hcat(stand_time,
                        fill(j,num_time_steps),
                        fill(i,num_time_steps),
                        fill(k,num_time_steps)),:auto)
                rename!(col_df, :x1 => :time, :x2 => :rep, :x3 => :GEM_ver, :x4 => :state_ID)
                x_temp = x_out_mat[:,:,k,j,i]
                x_temp = DataFrame(transpose(x_temp), :auto)
                rename!(x_temp, col_names)
                x_out_dat = hcat(col_df, x_temp)
                x_out_rep_store[j] = x_out_dat
            end
            x_out_gem_v_store[i] = vcat(x_out_rep_store...)
        end
        x_out_spp_store[k] = vcat(x_out_gem_v_store...)
    end
    par_mean_df = vcat(x_out_spp_store...)

    #######################

    x_out_var = vec(x_var_stand_out_all)

    x_out_var_mat = reshape(x_out_var, no_columns-1,num_time_steps, no_species,num_rep, length(GEM_ver))
    
    x_out_var_spp_store = Vector{DataFrame}(undef, no_species)
    x_out_var_gem_v_store = Vector{DataFrame}(undef, length(GEM_ver)) 
    x_out_var_rep_store = Vector{DataFrame}(undef, num_rep)

    for k = 1:no_species
        for i = 1:length(GEM_ver)
            for j = 1:num_rep
                col_df = DataFrame(hcat(stand_time,
                        fill(j,num_time_steps),
                        fill(i,num_time_steps),
                        fill(k,num_time_steps)),:auto)
                rename!(col_df, :x1 => :time, :x2 => :rep, :x3 => :GEM_ver, :x4 => :state_ID)
                x_temp = x_out_var_mat[:,:,k,j,i]
                x_temp = DataFrame(transpose(x_temp), :auto)
                rename!(x_temp, col_names)
                x_out_var_dat = hcat(col_df, x_temp)
                x_out_var_rep_store[j] = x_out_var_dat
            end
            x_out_var_gem_v_store[i] = vcat(x_out_var_rep_store...)
        end
        x_out_var_spp_store[k] = vcat(x_out_var_gem_v_store...)
    end
    par_var_df = vcat(x_out_var_spp_store...)

    return     pop_df, par_mean_df, par_var_df;
end