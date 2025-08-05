##############################################
#		        PLOT FUNCtions               #
##############################################
function Pop_Plot(pop_time::DataFrame, i::Int64)
    pop_time = pop_time[pop_time.state_ID .== i, :]
    rep_col_name = names(pop_time)[4:size(pop_time)[2]]
    pop_stack = stack(pop_time, 
                    rep_col_name)
    pop_plot = data(pop_stack) * mapping(
        :time,
        :value,
        #color = :replicate,
        group = :variable,
        row = :GEM_ver => nonnumeric 
        ) *
        visual(Lines, color = :grey)
    draw(pop_plot, axis=(xlabel="Time", ylabel="Population abundances"))    
end



function Trait_Plot(trait_mean::DataFrame, trait_var::DataFrame, 
    spID::Int64, trait_to_plot::String )

    #=srcatch
    meandf = deci_threshold(trait_mean_df)
    vardf = deci_threshold(trait_var_df)
    spID = 1
    trait_to_plot = "b_max"
    =#

    meandf = deci_threshold(trait_mean)
    vardf = deci_threshold(trait_var)
    #pick the right GEM version number and state ID
    #meandf = meandf[meandf.GEM_ver .== GEMver, :]
    meandf = meandf[meandf.state_ID .== spID, :]
    #vardf = vardf[vardf.GEM_ver .== GEMver, :]
    vardf = vardf[vardf.state_ID .== spID, :]

    #pick the trait wanted
    meantemp = meandf[:, [:time, :rep, :GEM_ver]]
    meantemp2 = meandf[:,trait_to_plot]
    mean2plot = hcat(meantemp, meantemp2)

    #grab the right var
    vartemp2 = vardf[:,trait_to_plot]
    
    dftemp = hcat(mean2plot, vartemp2, makeunique=true)
    rename!(dftemp, :x1 => :mean, :x1_1 => :var)

    upper_bound = dftemp.mean .+ dftemp.var
    lower_bound = dftemp.mean .- dftemp.var
    df2plot = hcat(dftemp, upper_bound, lower_bound, makeunique = true)
    rename!(df2plot,:x1 => :upper_bound, :x1_1 => :lower_bound )

    mean_plot = data(df2plot) * mapping(
        :time,
        :mean,
        #color = :rep,
        group = :rep => nonnumeric,
        row = :GEM_ver => nonnumeric
    ) *
        visual(Lines, color = :grey)

    var_plot = data(df2plot) * mapping(
        :time,
        :upper_bound,
        :lower_bound,
        group = :rep => nonnumeric,
        row = :GEM_ver => nonnumeric
    )    * visual(Band, color = :lightgrey, alpha=0.5)

    tempplot = data(df2plot) * (mean_plot + var_plot)
    finalplot = draw(tempplot, axis=(xlabel="Time", ylabel="Trait value"))

    return finalplot
end

#phase plane plot
