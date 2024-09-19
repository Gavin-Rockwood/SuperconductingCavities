import QuantumToolbox as qt
using YAXArrays
import ProgressMeter as PM
using Dates


function LoadRunResults(file; cube_order = "Default", h_dims = "Default")
    ds = open_dataset(file)
    n_cubes = length(ds.cubes)
    
    
    to_return = Dict{Any, Any}()
    
    if cube_order == "Default"
        if "Order" in keys(ds.properties)
            cube_order = [Symbol(s) for s in split(ds.properties["Order"])]
        end
        if cube_order == "Steps"
            cube_order = []
            for i in 1:n_cubes 
                push!(cube_order, Symbol("Step_"*string(i)))
            end
        end
    end
    if "dims" in keys(ds.properties)
        h_dims = ds.properties["dims"]
    elseif h_dims == "Default"
        h_dims = nothing
    end

    state_list = []
    time_list = []
    t = 0

    totalstepnum = 0
    for i in 1:n_cubes
        totalstepnum += length(dims(ds.cubes[Symbol(cube_order[i])])[1])
    end
    prog = PM.Progress(totalstepnum)
    for i in 1:n_cubes
        cube = Symbol(cube_order[i])
        num_steps = length(dims(ds.cubes[cube])[1])
        Re = collect(ds.cubes[cube][P = At("Re")].data)
        Im = collect(ds.cubes[cube][P = At("Im")].data)

        for j in 1:num_steps
            vec1 = selectdim(Re, 1, j)
            vec2 = selectdim(Im, 1, j)
            vec = vec1+vec2*1im
            push!(state_list, qt.Qobj(vec, dims = h_dims))

            if "Times" in keys(ds.cubes[cube].properties)
                push!(time_list, ds.cubes[cube].properties["Times"][j]+t)
            end

            b = now()
            PM.next!(prog)
        end
        t += time_list[end]
    end

    return [state_list, time_list]

end