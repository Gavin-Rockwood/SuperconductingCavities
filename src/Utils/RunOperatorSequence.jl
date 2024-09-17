import QuantumToolbox as qt
using YAXArrays
using Dates
using NetCDF
import OrdinaryDiffEq as ODE
import ProgressMeter as PM

function RunSingleOperator(model::Transmon_Resonators, 
    ψ::qt.QuantumObject{<:AbstractVector{T1},qt.KetQuantumObject}, 
    op_params; 
    spps = 5, 
    solver_kwargs = Dict{Any, Any}(), 
    save_step = false, 
    step_name = "DEFAULT", 
    op_name = "DEFAULT", 
    to_return = "All WFs", 
    save_path = "Data/", 
    run_name = "", 
    save_as_seperate_file = false, 
    tspan = [], 
    overlaps_progbar = true
    ) where T1<:Number
    #-------------------------------------------------------------

    step_name = replace(step_name, " " => "_")
    to_return_vals = ["All WFs", "Last WF", "Overlaps", "Nothing"]
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

    ν  = op_params["freq_d"]+op_params["shift"]
    ε = op_params["epsilon"]

    drive_coef = Get_Drive_Coef(ν, ε, envelope = Get_Envelope(op_params["Envelope"], op_params["Envelope Args"]))

    Ĥ_D = Get_Ĥ_D(model.n̂ₜ, drive_coef)


    if length(tspan) == 0
        tspan = collect(LinRange(0, op_params["pulse_time"], Int(ceil(op_params["pulse_time"]*spps))+1))
    end

    solver_kwargs_sym = Dict{Symbol, Any}()
    for key in keys(solver_kwargs)
        solver_kwargs_sym[Symbol(key)] = solver_kwargs[key]
    end
    @info "Running Time Evolution"
    sleep(0.01)
    res = qt.sesolve(2*π*model.Ĥ, ψ, tspan, H_t = Ĥ_D; solver_kwargs_sym...)
    @debug println(tostr(res))
    
    @info "Time Evolution Complete"
    sleep(0.01)

    if (save_step) | (to_return == "Overlaps")
        @info "Getting Overlaps"
        sleep(0.01)
        states = string.(collect(keys(model.dressed_states)))
        steps_symbol = Symbol(step_name*"_Steps")
        axlist = (Dim{:State}(states), Dim{steps_symbol}(collect(0:length(res.times)-1)), Dim{:Comp}(["Re", "Im"]))
        
        properties = Dict{Any, Any}("Data"=>"Overlaps", "Time"=>string(now()), "Times" => res.times, "Operator" => op_name)
        data = zeros(size(axlist))

        ds_properties = Dict{Any, Any}("Basis Shape" => collect(model.hilbertspace.Ĥ.dims), "Data"=>"Wave Function Amplitudes", "Solver_Args" => string(solver_kwargs))
        overlaps = Dataset(; Symbol(step_name) => YAXArray(axlist, data, properties), properties = ds_properties)

        coords = Iterators.product(axlist...)
        @debug print(tostr(overlaps.cubes[Symbol(step_name)]))
        p = PM.Progress(length(coords), enabled = overlaps_progbar)
        for coord in coords            
            state = coord[1]
            step = coord[2]
            ψt = res.states[step+1]
            ψ = model.dressed_states[eval(Meta.parse(state))]
            overlap = ψ'*ψt
            overlaps.cubes[Symbol(step_name)][At.([state, step, "Re"])...] = real(overlap)
            overlaps.cubes[Symbol(step_name)][At.([state, step, "Im"])...] = imag(overlap)
            PM.next!(p)
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
        return res.states[end]
    elseif to_return == "Overlaps"
        return overlaps
    end
end

function RunSingleOperator(model::Transmon_Resonators, 
    ρ::qt.QuantumObject{<:AbstractArray{T1},qt.OperatorQuantumObject}, 
    op_params; 
    c_ops = nothing,  
    spps = 5, 
    solver_kwargs = Dict{Any, Any}(), 
    save_step = false, 
    step_name = "DEFAULT", 
    op_name = "DEFAULT",
    to_return = "All DMs", 
    save_path = "Data/", 
    run_name = "", 
    save_as_seperate_file = false, 
    tspan = []
    ) where T1<:Number
    #-------------------------------------------------------------
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

    ν  = op_params["freq_d"]+op_params["shift"]
    ε = op_params["epsilon"]
    drive_coef = Get_Drive_Coef(ν, ε, envelope = Get_Envelope(op_params["Envelope"], op_params["Envelope Args"]))

    Lₜ = Get_Lₜ(model.n̂ₜ, drive_coef)

    if length(tspan) == 0
        tspan = collect(LinRange(0, op_params["pulse_time"], Int(ceil(op_params["pulse_time"]*spps))+1))
    end
    
    solver_kwargs_sym = Dict{Symbol, Any}()
    for key in keys(solver_kwargs)
        solver_kwargs_sym[Symbol(key)] = solver_kwargs[key]
    end
    @info "Running Time Evolution"
    sleep(0.01)
    res = qt.mesolve(2*π*model.Ĥ, ρ, tspan, c_ops, H_t = Lₜ; solver_kwargs_sym...)
    @info "Time Evolution Complete"
    sleep(0.01)

    if (save_step) | (to_return == "Probabilities")
        @info "Getting Probabilities"
        sleep(0.01)
        states = string.(collect(keys(model.dressed_states)))
        steps_symbol = Symbol(step_name*"_Steps")
        axlist = (Dim{:State}(states), Dim{steps_symbol}(collect(0:length(res.times)-1)))
        
        properties = Dict{Any, Any}("Data"=>"Probabilities", "Time"=>string(now()), "Times" => res.times, "Operator" => op_name)
        data = zeros(size(axlist))
        
        ds_properties = Dict{Any, Any}("Basis Shape" => collect(model.hilbertspace.Ĥ.dims), "Data"=>"Probabilities From Density Matrix", "Solver_Args"=>string(solver_kwargs))
        probabilities = Dataset(; Symbol(step_name) => YAXArray(axlist, data, properties), properties = ds_properties)

        coords = Iterators.product(axlist...)
        for coord in coords            
            state = coord[1]
            step = coord[2]
            ρt = res.states[step+1]
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
        return res.states[end]
    elseif to_return == "Probabilities"
        return probabilities
    end
end

function RunPulseSequence(model::Transmon_Resonators, 
    ψ::qt.QuantumObject{<:AbstractVector{T1}, qt.KetQuantumObject, 2}, 
    op_sequence; 
    spps = 5, 
    solver_kwargs = Dict{Any, Any}(), 
    run_name = "", 
    save_path = "Data/"
    ) where T1<:Number
    #-------------------------------------------------------------
    
    if run_name == ""
        run_name = "Operator_Sequence_"*string(now())
    end

    for i in 1:length(op_sequence)
        op = op_sequence[i]
        @info "Running operator $op"
        ψ = RunSingleOperator(model, ψ, model.Stuff["op_drive_params"][op]; spps = spps, solver_kwargs = solver_kwargs, run_name = run_name, save_path = save_path, step_name = "Step_"*string(i), to_return = "Last WF", save_step = true, op_name = op)
    end
end


function RunPulseSequence(model::Transmon_Resonators, 
    ρ::qt.QuantumObject{<:AbstractArray{T1},qt.OperatorQuantumObject}, 
    op_sequence; 
    spps = 5, 
    solver_kwargs = Dict{Any, Any}(), 
    run_name = "", 
    save_path = "Data/", 
    c_ops = []
    ) where T1<:Number
    #-------------------------------------------------------------
    if run_name == ""
        run_name = "Operator_Sequence_"*string(now())
    end

    for i in 1:length(op_sequence)
        op = op_sequence[i]
        @info "Running operator $op"
        ρ = RunSingleOperator(model, ρ, model.Stuff["op_drive_params"][op]; spps = spps, solver_kwargs = solver_kwargs, run_name = run_name, save_path = save_path, step_name = "Step_"*string(i), to_return = "Last DM", c_ops = c_ops, save_step = true, op_name = op)
    end
end



