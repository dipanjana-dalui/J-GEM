##############################################
#		        PLOT FUNCtions               #
##############################################
function Pop_Plot(pop_time::DataFrame, show::Bool)
    pop_stack = stack(pop_time, Not([:time, :GEM_ver, :state_ID]); 
                    variable_name = :replicate, value_name = :value)

    pop_plot = data(pop_stack) * mapping(
        :time,
        :value,
        #color = :replicate,
        group = :replicate,
        row = :GEM_ver => nonnumeric
        ) *
        visual(Lines)
    if show
        return draw(pop_plot)
    else
        return nothing
    end
end

#=
pop_time_series_df = GEM_run[1]
trait_mean_df = GEM_run[2]
trait_var_df = GEM_run[3]

using Tidier
@glimpse(trait_mean_df)

dat_t_mean_Gv2 = @filter(trait_mean_df, GEM_ver==2.0)
@filter(dat_t_mean_Gv2, state_ID == 1 )
dat_t_mean_Gv2.b_max
=#