import QuantumOptics as qo
using YAXArrays
using Dates
using NetCDF

function RunSingleOperator(model::Transmon_Resonators, ψ:: qo.Ket, op_params; spps = 5, solver_kwargs = Dict{Any, Any}(), save_step = true, step_name = "DEFAULT", to_return = "All WFs", save_path = "Data/", run_name = "", save_as_seperate_file = false)
    step_name = replace(step_name, " " => "_")
    to_return_vals = ["All WFs", "Last WF", "Overlaps"]
    if !(to_return in to_return_vals)
        println("to_return (\"$to_return\") not one of the accepted values: ")
        for val in to_return_vals
            println("   $val")
        end
    end
    if run_name == ""
        run_name = "Single_Operator_Run_"*string(now())
    end

    if !("reltol" in keys(solver_kwargs))
        solver_kwargs["reltol"] = 1e-8
    end
    if !("abstol" in keys(solver_kwargs))
        solver_kwargs["abstol"] = 1e-8
    end
    if !("adaptive" in keys(solver_kwargs))
        solver_kwargs["adaptive"] = true
    end
    if !("tol" in keys(solver_kwargs))
        solver_kwargs["tol"] = 1e-5
    end
    if !("progress" in keys(solver_kwargs))
        solver_kwargs["progress"] = true
    end
    ν  = op_params["freq_d"]+op_params["shift"]


    Ĥt = f_for_schroedinger_dynamic(model, model.n̂ₜ, ν, op_params["epsilon"], envelope_name = op_params["Envelope"], envelope_kwargs = op_params["Envelope Args"])

    tspan = collect(LinRange(0, op_params["pulse_time"], Int(ceil(op_params["pulse_time"]*spps))+1))

    solver_kwargs_sym = Dict{Symbol, Any}()
    for key in keys(solver_kwargs)
        solver_kwargs_sym[Symbol(key)] = solver_kwargs[key]
    end
    @info "Running Time Evolution"
    sleep(0.01)
    res = qo.timeevolution.schroedinger_dynamic(tspan, ψ, Ĥt; solver_kwargs_sym...)
    @info "Time Evolution Complete"
    sleep(0.01)

    if (save_step) | (to_return == "Overlaps")
        @info "Getting Overlaps"
        sleep(0.01)
        states = string.(collect(keys(model.dressed_states)))
        steps_symbol = Symbol(step_name*"_Steps")
        axlist = (Dim{:State}(states), Dim{steps_symbol}(collect(0:length(tspan)-1)), Dim{:Comp}(["Re", "Im"]))
        
        properties = Dict{Any, Any}("Data"=>"Overlaps", "Time"=>string(now()), "Times" => tspan)
        data = zeros(size(axlist))
        overlaps = Dataset(; Symbol(step_name) => YAXArray(axlist, data, properties))

        coords = Iterators.product(axlist...)
        for coord in coords            
            state = coord[1]
            step = coord[2]
            ψt = res[2][step+1]
            ψ = model.dressed_states[eval(Meta.parse(state))]
            overlap = ψ'*ψt
            overlaps.cubes[Symbol(step_name)][At.([state, step, "Re"])...] = real(overlap)
            overlaps.cubes[Symbol(step_name)][At.([state, step, "Im"])...] = imag(overlap)
        end
    end

    if save_step
        file_name = save_path*run_name
        if save_as_seperate_file
            file_name = file_name*step_name
        end


        if !isdir(save_path)
            @info "Making Save Directory $save_path"
            mkdir(save_path)
        end
        savedataset(overlaps,path = file_name*".nc", append = true)
    end

    @info "Done with $step_name"

    if to_return == "All WFs"
        return res
    elseif to_return == "Last WF"
        return res[2][end]
    elseif to_return == "Overlaps"
        return overlaps
    end
end

function RunSingleOperator(model::Transmon_Resonators, ρ:: qo.Operator, op_params; c_ops = [],  spps = 5, solver_kwargs = Dict{Any, Any}(), save_step = true, step_name = "DEFAULT", to_return = "All DMs", save_path = "Data/", run_name = "", save_as_seperate_file = false)
    step_name = replace(step_name, " " => "_")
    to_return_vals = ["All DMs", "Last DM", "Probabilities"]
    if !(to_return in to_return_vals)
        println("to_return (\"$to_return\") not one of the accepted values: ")
        for val in to_return_vals
            println("   $val")
        end
        return nothing
    end
    if run_name == ""
        run_name = "Single_Operator_Run_"*string(now())
    end

    if !("reltol" in keys(solver_kwargs))
        solver_kwargs["reltol"] = 1e-8
    end
    if !("abstol" in keys(solver_kwargs))
        solver_kwargs["abstol"] = 1e-8
    end
    if !("adaptive" in keys(solver_kwargs))
        solver_kwargs["adaptive"] = true
    end
    if !("tol" in keys(solver_kwargs))
        solver_kwargs["tol"] = 1e-5
    end
    if !("progress" in keys(solver_kwargs))
        solver_kwargs["progress"] = true
    end
    ν  = op_params["freq_d"]+op_params["shift"]


    Ĥt = f_for_master_dynamic(model, model.n̂ₜ, ν, op_params["epsilon"], c_ops = c_ops, envelope_name = op_params["Envelope"], envelope_kwargs = op_params["Envelope Args"])

    tspan = collect(LinRange(0, op_params["pulse_time"], Int(ceil(op_params["pulse_time"]*spps))+1))

    solver_kwargs_sym = Dict{Symbol, Any}()
    for key in keys(solver_kwargs)
        solver_kwargs_sym[Symbol(key)] = solver_kwargs[key]
    end
    @info "Running Time Evolution"
    sleep(0.01)
    res = qo.timeevolution.master_dynamic(tspan, ρ, Ĥt; solver_kwargs_sym...)
    @info "Time Evolution Complete"
    sleep(0.01)

    if (save_step) | (to_return == "Probabilities")
        @info "Getting Probabilities"
        sleep(0.01)
        states = string.(collect(keys(model.dressed_states)))
        steps_symbol = Symbol(step_name*"_Steps")
        axlist = (Dim{:State}(states), Dim{steps_symbol}(collect(0:length(tspan)-1)))
        
        properties = Dict{Any, Any}("Data"=>"Probabilities", "Time"=>string(now()), "Times" => tspan)
        data = zeros(size(axlist))
        probabilities = Dataset(; Symbol(step_name) => YAXArray(axlist, data, properties))

        coords = Iterators.product(axlist...)
        for coord in coords            
            state = coord[1]
            step = coord[2]
            ρt = res[2][step+1]
            ψ = model.dressed_states[eval(Meta.parse(state))]
            probability = abs(ψ'*ρt*ψ)
            probabilities.cubes[Symbol(step_name)][At.([state, step])...] = probability
        end
    end

    if save_step
        file_name = save_path*run_name
        if save_as_seperate_file
            file_name = file_name*step_name
        end


        if !isdir(save_path)
            @info "Making Save Directory $save_path"
            mkdir(save_path)
        end
        savedataset(probabilities,path = file_name*".nc", append = true)
    end

    @info "Done with $step_name"
    if to_return == "All DMs"
        return res
    elseif to_return == "Last DM"
        return res[2][end]
    elseif to_return == "Probabilities"
        return probabilities
    end
end

function RunPulseSequence(model::Transmon_Resonators, ψ::qo.Ket, op_sequence; spps = 5, solver_kwargs = Dict{Any, Any}(), run_name = "", save_path = "Data/")
    if run_name == ""
        run_name = "Operator_Sequence_"*string(now())
    end

    for i in 1:length(op_sequence)
        op = op_sequence[i]
        RunSingleOperator(model, ψ, model.Stuff["op_drive_params"][op]; spps = spps, solver_kwargs = solver_kwargs, run_name = run_name, save_path = save_path, step_name = "Step_"*string(i), to_return = "All WFs")
    end
end


function RunPulseSequence(model::Transmon_Resonators, ρ::qo.Operator, op_sequence; spps = 5, solver_kwargs = Dict{Any, Any}(), run_name = "", save_path = "Data/")
    if run_name == ""
        run_name = "Operator_Sequence_"*string(now())
    end

    for i in 1:length(op_sequence)
        op = op_sequence[i]
        RunSingleOperator(model, ρ, model.Stuff["op_drive_params"][op]; spps = spps, solver_kwargs = solver_kwargs, run_name = run_name, save_path = save_path, step_name = "Step_"*string(i), to_return = "All DMs")
    end
end



