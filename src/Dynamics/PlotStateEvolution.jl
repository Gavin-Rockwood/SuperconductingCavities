import QuantumToolbox as qt
import CairoMakie as cm
using YAXArrays

markers = ['+', '×', '∘', '⋆', '▿', '▵', '⊲', '⊳', '⏣', '⭔']


function PlotSingleModeEvolution(model, 
    tlist, 
    state_hist; 
    markers = markers,
    plot_every = 1,
    figure_kwargs = Dict{Any, Any}(),
    axis_kwargs = Dict{Any, Any}(),
    cmap_name = :seaborn_bright,
    scatterlines_kwargs = Dict{Any, Any}(),
    show_thresh = 1e-3, 
    emph_states = [], # This setting is for emphasizing states. If left empty all will have alpha = 1, otherwise, only the states is this list will have alpha = 1
    non_emph_alpha = 0.1) # the rest will have this alpha


    if !("size" in keys(figure_kwargs))
        figure_kwargs["size"] = (1000, 500)
    end
    if !("figure_padding" in keys(figure_kwargs))
        figure_kwargs["figure_padding"] = 10
    end
    if !("px_per_unit" in keys(figure_kwargs))
        figure_kwargs["px_per_unit"] = 4
    end

    if !("title" in keys(axis_kwargs))
        axis_kwargs["title"] = "Pulse Sequence Results"
    end
    if !("xlabel" in keys(axis_kwargs))
        axis_kwargs["xlabel"] = "Time"
    end
    if !("ylable" in keys(axis_kwargs))
        axis_kwargs["ylabel"] = "Probability"
    end

    if !("markersize" in keys(scatterlines_kwargs))
        scatterlines_kwargs["markersize"] = 20
    end
    if !("linewidth" in keys(scatterlines_kwargs))
        scatterlines_kwargs["linewidth"] = 0.2
    end

    # if !(:orientation in keys(legend_kwargs))
    #     legend_kwargs[:orientation] = :vertical
    # end
    # if !(:position in keys(legend_kwargs))
    #     legend_kwargs[:position] = :ct
    # end


    if typeof(state_hist) <: Dict
        EVs = state_hist
    else
        proj_dict = Utils.Get_Projection_Ops(model.dressed_states)
        EVs = Utils.Get_EVs(state_hist, proj_dict)
    end
    
    if length(emph_states) == 0
        emph_states = keys(model.dressed_states)
    end

    @info "Making Plot"
    f = cm.Figure(size = figure_kwargs["size"], figure_padding = figure_kwargs["figure_padding"], px_per_unit = figure_kwargs["px_per_unit"])

    ax = cm.Axis(f[1,1], xlabel = axis_kwargs["xlabel"], ylabel = axis_kwargs["ylabel"], title = axis_kwargs["title"])

    tlevels = model.Ĥ.dims[1]
    clevels = model.Ĥ.dims[2]
    cmap=(cm.cgrad(cmap_name, tlevels, categorical = true))

    #return(EVs)

    for t in 0:(tlevels-1)
        for c in 0:(clevels-1)
            state = (t, c)
            label = nothing
                if c == 1
                    label = string(t)
                end
            if (state in keys(EVs)) | (string(state) in keys(EVs))
                if string(state) in keys(EVs)
                    state = string(state)
                end

                x = tlist[1:plot_every:end]
                y = EVs[state][1:plot_every:end]

                if !(tlist[end] in x)
                    push!(x, tlist[end])
                    push!(y, EVs[state][end])
                end

                alpha = non_emph_alpha
                if maximum(abs.(y)) < show_thresh
                    alpha = 0.0
                end
                if (state in emph_states) | (state in string.(emph_states))
                    alpha = 1
                end

                cm.scatterlines!(ax, x, real.(y), marker = markers[t+1], color = (cmap[c+1], alpha), markersize = scatterlines_kwargs["markersize"], linewidth = scatterlines_kwargs["linewidth"])
            end
            cm.scatterlines!(ax, [0], [0], marker = markers[t+1], color = (cmap[c+1]), label = label, markersize = scatterlines_kwargs["markersize"], linewidth = scatterlines_kwargs["linewidth"], visible = false)
        end
    end

    #for i in 1:tlevels
    #    t = i-1
    #    label = string(t)
    #    cm.scatterlines!(ax, [0], [0], marker = markers[t+1], color = cmap[1], label = label, markersize = scatterlines_kwargs["markersize"], linewidth = scatterlines_kwargs["linewidth"])
    #end

    tick_loc = collect(0:tlevels-1).+0.5
    tick_name = string.(collect(0:tlevels-1))
    cm.Colorbar(f[1,2], colormap = cmap, limits = (0,tlevels), ticks = (tick_loc, tick_name), label = "Photon Number")

    label = "Transmon\nLevels"
    #if legend_kwargs[:orientation] == :vertical
    #    label = "Transmon\nLevels"
    #end
    f[1,3] = cm.Legend(f, ax, label, framevisible = false)#; legend_kwargs...)#orientation = :horizontal, position = :ct)
    cm.ylims!(ax, -0.1, 1.1)

    cm.display(f)
end