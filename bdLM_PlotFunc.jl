##############################################
#		  DATA CLEANING & PLOT FUNC          #
##############################################

using  AlgebraOfGraphics
using Tidier


########




#=
out_dat_ver1 = DataFrame(out_mat[:,:,1], :auto)

dat = hcat(stand_time, out_dat_ver1, makeunique=true)

plot(dat.x1, dat.x1_1)
plot!(dat.x1, dat.x2)
plot!(dat.x1, dat.x3)
plot!(dat.x1, dat.x4)
plot!(dat.x1, dat.x5)
=#