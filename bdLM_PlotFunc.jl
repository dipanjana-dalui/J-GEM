##############################################
#		  DATA CLEANING & PLOT FUNC          #
##############################################

using  AlgebraOfGraphics
using Tidier




########
x_out = vec(x_stand_out_all)
#out_x_var = vec(x_var_stand_out_all)
typeof(x_out)

x_out_mat = reshape(x_out, no_columns-1,num_time_steps, no_species,num_rep, length(GEM_ver))
typeof(x_out_mat)
x_out_mat[:,:,1,:,1]

DataFrame(x_out_mat[:,:,1,2,1],:auto)

all_x_out = Vector{DataFrame}(undef, length(GEM_ver)*num_rep) 
x_out_rep_store = Vector{DataFrame}(undef, num_rep)
x_out_spp_store = Vector{DataFrame}(undef, spp)


for k = 1:no_species
    k=1
    #for i = 1:length(GEM_ver)
        i=1
    #    for j = 1:num_rep
        j = 2
        temp = DataFrame(hcat(stand_time,
                        fill(j,num_time_steps),
                        fill(GEM_ver[i],num_time_steps)),:auto)
        rename!(temp, :x1 => :time, :x2 => :rep, :x3 => :GEM_ver)
        x_out_temp = DataFrame(x_out_mat[:,:,k,j,i],:auto)
        transpose(x_out_temp)
        x_out_temp = vcat(temp, x_out_temp)
        x_out_rep_store[j] = x_out_temp

        end
    end

end

#=
out_dat_ver1 = DataFrame(out_mat[:,:,1], :auto)

dat = hcat(stand_time, out_dat_ver1, makeunique=true)

plot(dat.x1, dat.x1_1)
plot!(dat.x1, dat.x2)
plot!(dat.x1, dat.x3)
plot!(dat.x1, dat.x4)
plot!(dat.x1, dat.x5)
=#