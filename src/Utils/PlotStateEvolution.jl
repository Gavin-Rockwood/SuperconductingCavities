import QuantumToolbox as qt
import CairoMakie as cm
using YAXArrays

markers = ['+', '×', '∘', '⋆', '▿', '▵', '⊲', '⊳', '⏣', '⭔']


function ProbabilityPlotSingleModeEvolution(model, data :: Dataset; markers = markers, plot_every = 10, figure_kwargs = Dict{Any, Any}(), Axis_kwargs = Dict{Any, Any}(), cmap_name = :jet1, scatterlines_kwargs = Dict{Any, Any}())

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


    stepped = true
    if !(:Step_1 in keys(data.cubes))
        stepped = false
    end
    t0 = 0
    times = []
    ys = []
    data_keys = collect(keys(data.cubes))
    @info "Organizing Data"
    for i in 1:length(data.cubes)
        #@info "On Step $i"
        if stepped
            step_sym = Symbol("Step_$i")
        else
            step_sym = data_keys[i]
        end
        times = [times; t0 .+ data.cubes[step_sym].properties["Times"]]
        t0 = times[end]
        temp_data = collect(data.cubes[step_sym].data)
        if !(dims(data.cubes[step_sym], :Comp) == nothing)
            temp_data = temp_data[:, :, 1].^2+temp_data[:, :, 2].^2
        end
        push!(ys, temp_data)
    end
    ys = cat(ys...; dims = 2); 

    @info "Making Plot"
    f = cm.Figure(size = figure_kwargs["size"], figure_padding = figure_kwargs["figure_padding"], px_per_unit = figure_kwargs["px_per_unit"])

    ax = cm.Axis(f[1,1], xlabel = Axis_kwargs["xlabel"], ylabel = Axis_kwargs["ylabel"], title = Axis_kwargs["title"])

    tlevels = model.Ĥ.dims[1]
    cmap=(cm.cgrad(cmap_name, tlevels, categorical = true))

    states = data.cubes[data_keys[1]].axes[1]
    for i in 1:length(states)
        t = eval(Meta.parse(states[i]))[1]
        c = eval(Meta.parse(states[i]))[2]
        x = times[1:plot_every:end]
        y = ys[i, :][1:plot_every:end]

        cm.scatterlines!(ax, x, y, marker = markers[t+1], color = cmap[c+1], markersize = scatterlines_kwargs["markersize"], linewidth = scatterlines_kwargs["linewidth"])
    end

    for i in 1:tlevels
        t = i-1
        label = string(t)
        cm.scatterlines!(ax, [0], [0], marker = markers[t+1], color = cmap[1], label = label, markersize = scatterlines_kwargs["markersize"], linewidth = scatterlines_kwargs["linewidth"])
    end

    tick_loc = collect(0:tlevels-1).+0.5
    tick_name = string.(collect(0:tlevels-1))
    cm.Colorbar(f[1,2], colormap = cmap, limits = (0,tlevels), ticks = (tick_loc, tick_name), label = "Photon Number")

    cm.axislegend(ax, "Transmon Levels", orientation = :horizontal, position = :ct)

    f


end