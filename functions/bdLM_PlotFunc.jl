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



