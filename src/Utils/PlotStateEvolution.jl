import QuantumToolbox as qt
import CairoMakie as cm
using YAXArrays

markers = ['+', '×', '∘', '⋆', '▿', '▵', '⊲', '⊳', '⏣', '⭔']

function Get_Projection_Ops(dict_of_wavefunctions)
    res = deepcopy(dict_of_wavefunctions)
    map!(x -> x*x', values(res))
    return res
end


function Get_EVs(list_of_states, dict_of_operators)
    res = Dict{Any, Any}()

    for key in keys(dict_of_operators)
        op = dict_of_operators[key]

        res[key] = [qt.expect(op, ψ) for ψ in list_of_states]
    end
    return res
end

function PlotSingleModeEvolution(model, tlist, state_hist; markers = markers, plot_every = 1, figure_kwargs = Dict{Any, Any}(), Axis_kwargs = Dict{Any, Any}(), cmap_name = :jet1, scatterlines_kwargs = Dict{Any, Any}(), show_thresh = 1e-3)

    if !("size" in keys(figure_kwargs))
        figure_kwargs["size"] = (1000, 500)
    end
    if !("figure_padding" in keys(figure_kwargs))
        figure_kwargs["figure_padding"] = 10
    end
    if !("px_per_unit" in keys(figure_kwargs))
        figure_kwargs["px_per_unit"] = 4
    end

    if !("title" in keys(Axis_kwargs))
        Axis_kwargs["title"] = "Pulse Sequence Results"
    end
    if !("xlabel" in keys(Axis_kwargs))
        Axis_kwargs["xlabel"] = "Time"
    end
    if !("ylable" in keys(Axis_kwargs))
        Axis_kwargs["ylabel"] = "Probability"
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
        proj_dict = Get_Projection_Ops(model.dressed_states)
        EVs = Get_EVs(state_hist, proj_dict)
    end
    

    @info "Making Plot"
    f = cm.Figure(size = figure_kwargs["size"], figure_padding = figure_kwargs["figure_padding"], px_per_unit = figure_kwargs["px_per_unit"])

    ax = cm.Axis(f[1,1], xlabel = Axis_kwargs["xlabel"], ylabel = Axis_kwargs["ylabel"], title = Axis_kwargs["title"])

    tlevels = model.Ĥ.dims[1]
    clevels = model.Ĥ.dims[2]
    cmap=(cm.cgrad(cmap_name, tlevels, categorical = true))

    #return(EVs)

    for t in 0:(tlevels-1)
        for c in 0:(clevels-1)
            state = (t, c)
            x = tlist[1:plot_every:end]
            y = EVs[state][1:plot_every:end]
            label = nothing
            if c == 1
                label = string(t)
            end
            alpha = 1
            if maximum(abs.(y)) < show_thresh
                alpha = 0.0
            end
            cm.scatterlines!(ax, x, real.(y), marker = markers[t+1], color = (cmap[c+1], alpha), markersize = scatterlines_kwargs["markersize"], linewidth = scatterlines_kwargs["linewidth"])
            cm.scatterlines!(ax, x[1], real.(y)[1], marker = markers[t+1], color = (cmap[c+1]), label = label, markersize = scatterlines_kwargs["markersize"], linewidth = scatterlines_kwargs["linewidth"], visible = false)
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

    f
end