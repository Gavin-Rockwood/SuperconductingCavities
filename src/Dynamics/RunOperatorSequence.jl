import QuantumToolbox as qt
using YAXArrays
using Dates
using NetCDF
import OrdinaryDiffEq as ODE
import ProgressMeter as PM
using OrdinaryDiffEqVerner: Vern9

function RunSingleOperator(Ĥ::qt.QuantumObject, Ô_D::qt.QuantumObject, 
    ψ::qt.QuantumObject{<:AbstractVector{T1},qt.KetQuantumObject}, 
    op_params; 
    spns = 5, 
    solver_kwargs = Dict{Any, Any}(), 
    save_step = false, 
    step_name = "DEFAULT", 
    op_name = "DEFAULT", 
    to_return = "All", 
    save_path = "Data/", 
    run_name = "", 
    save_as_seperate_file = false, 
    tspan = [],
    strob_skip = 1,
    other_ds_properties = Dict{Any, Any}(),
    use_logging = true,
    progress_bar = true,
    return_drive_coef = false # this is a flag to return the drive coef instead of running the time evolution. Used for debugging!
    ) where T1<:Number
    #-------------------------------------------------------------

    step_name = replace(step_name, " " => "_")
    to_return_vals = ["All", "Last", "Nothing"]
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

    if !("alg" in keys(solver_kwargs))
        solver_kwargs["alg"] = Vern9()
    end
    if !("abstol" in keys(solver_kwargs))
        solver_kwargs["abstol"] = 1e-6
    end 
    if !("reltol" in keys(solver_kwargs))
        solver_kwargs["reltol"] = 1e-6
    end 

    if !(progress_bar in keys(solver_kwargs))
        solver_kwargs["progress_bar"] = progress_bar
    end 
    
    chirp = false
    if "chirp_params" in keys(op_params)
        if !(op_params["chirp_params"] == nothing)
            chirp = true
        end
    end

    
    ν0  = op_params["freq_d"]
    if chirp
        @info "Using Chirp"
        α = 1
        if "chirp_factor" in keys(op_params)
            α = op_params["chirp_factor"]
        end
        ν = chirper(ν0, op_params["chirp_params"])# +(1-α)*op_params["shift"]+op_params["shift"]
    else
        ν = ν0.+op_params["shift"]
    end

    ε = op_params["epsilon"]

    digitize = false
    if "digitize" in keys(op_params)
        digitize = op_params["digitize"]
    end
    step_length = 2.3
    if "step_length" in keys(op_params)
        step_length = op_params["step_length"]
    end

    envelope = Envelopes.Get_Envelope(op_params["Envelope"], op_params["Envelope Args"], digitize = digitize, step_length = step_length)

    drive_coef = Get_Drive_Coef(ν, ε; envelope = envelope, drive_time = op_params["pulse_time"])
    
    if "filter_params" in keys(op_params)
        filter_params = Dict(Symbol(key)=>val for (key, val) in op_params["filter_params"]) # turns the string keys into symbols. 
        drive_coef = Get_Low_Pass_Filtered_Drive_Coef(drive_coef, op_params["pulse_time"]; filter_params...)
    end
    
    if return_drive_coef
        return drive_coef
    end
    Ĥ_D = Get_Ĥ_D(Ô_D, drive_coef)


    if (length(tspan) == 0) & (spns != "Stroboscopic")
        tspan = collect(LinRange(0, op_params["pulse_time"], Int(ceil(op_params["pulse_time"]*spns))+1))
    end

    if spns == "Stroboscopic"
        tspan = Get_Stroboscopic_Times(op_params)[1:strob_skip:end]
    end

    solver_kwargs_sym = Dict{Symbol, Any}()
    for key in keys(solver_kwargs)
        solver_kwargs_sym[Symbol(key)] = solver_kwargs[key]
    end
    if (use_logging) @info "Running Time Evolution" end
    sleep(0.01)
    res = qt.sesolve(2*π*Ĥ, ψ, tspan, H_t = Ĥ_D; solver_kwargs_sym...)
    @debug println(tostr(res))
    
    if (use_logging)  @info "Time Evolution Complete" end
    sleep(0.01)

    if (save_step) | (to_return == "DS")
        @info "Saving Steps"
        sleep(0.01)
        
        properties = Dict{Any, Any}("Data"=>"Wave Functions", "Time"=>string(now()), "Times" => res.times, "Operator" => op_name, "Operator_Parameters" => string(op_params))
        
        ds_properties = Dict{Any, Any}("dims" => collect(Ĥ.dims), "Data"=>"Wave Functions", "Solver_Args" => string(solver_kwargs))
        ds_properties = merge(ds_properties, other_ds_properties)

        ds = Utils.Qobj_List_To_DS(res.states; cube_name = Symbol(step_name), cube_properties = properties, ds_properties = ds_properties, step_name = Symbol(step_name*"_steps"))
        
        if save_step
            file_name = save_path*run_name
            if save_as_seperate_file
                file_name = file_name*step_name
            end

            if !isdir(save_path)
                @info "Making Save Directory $save_path"
                mkdir(save_path)
            end
            savedataset(ds,path = file_name*".nc", append = true)
        end
    end

    if (use_logging) @info "Done with $step_name" end

    if to_return == "All"
        return res
    elseif to_return == "Last"
        return res.states[end]
    end
end

"""
    RunSingleOperator(Ĥ::qt.QuantumObject, Ô_D::qt.QuantumObject, 
    ρ::qt.QuantumObject{<:AbstractArray{T1},qt.OperatorQuantumObject}, 
    op_params; 
    c_ops = nothing,  
    spns = 5, 
    solver_kwargs = Dict{Any, Any}(), 
    save_step = false, 
    step_name = "DEFAULT", 
    op_name = "DEFAULT",
    to_return = "All", 
    save_path = "Data/", 
    run_name = "", 
    save_as_seperate_file = false, 
    tspan = [],
    strob_skip = 1,
    other_ds_properties = Dict{Any, Any}()
    ) where T1<:Number

TBW
"""
function RunSingleOperator(Ĥ::qt.QuantumObject, Ô_D::qt.QuantumObject, 
    ρ::qt.QuantumObject{<:AbstractArray{T1},qt.OperatorQuantumObject}, 
    op_params; 
    c_ops = nothing,  
    spns = 5, 
    solver_kwargs = Dict{Any, Any}(), 
    save_step = false, 
    step_name = "DEFAULT", 
    op_name = "DEFAULT",
    to_return = "All", 
    save_path = "Data/", 
    run_name = "", 
    save_as_seperate_file = false, 
    tspan = [],
    strob_skip = 1,
    other_ds_properties = Dict{Any, Any}()
    ) where T1<:Number
    #-------------------------------------------------------------
    step_name = replace(step_name, " " => "_")
    to_return_vals = ["All", "Last", "Nothing", "DS"]
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

    if !("alg" in keys(solver_kwargs))
        solver_kwargs["alg"] = Vern9()
    end
    if !("abstol" in keys(solver_kwargs))
        solver_kwargs["abstol"] = 1e-6
    end 
    if !("reltol" in keys(solver_kwargs))
        solver_kwargs["reltol"] = 1e-6
    end 

    chirp = false
    if "chirp_params" in keys(op_params)
        if !(op_params["chirp_params"] == nothing)
            chirp = true
        end
    end

    ν  = op_params["freq_d"].+op_params["shift"]
    ε = op_params["epsilon"]

    if chirp
        @info "Using Chirp"
        α = 1
        if "chirp_factor" in keys(op_params)
            α = op_params["chirp_factor"]
        end
        ν(ε) = op_params["freq_d"] + α*sum(op_params["chirp_params"][n]*ε^n for n in 1:length(op_params["chirp_params"]))+(1-α)*op_params["shift"]
    end

    digitize = false
    if "digitize" in keys(op_params)
        digitize = op_params["digitize"]
    end
    step_length = 0
    if "step_length" in keys(op_params)
        step_length = op_params["step_length"]
    end

    envelope = Envelopes.Get_Envelope(op_params["Envelope"], op_params["Envelope Args"], digitize = digitize, step_length = step_length)
    drive_coef = Get_Drive_Coef(ν, ε, envelope = envelope)

    if "filter_params" in keys(op_params)
        drive_coef = Get_Low_Pass_Filtered_Drive_Coef(drive_coef, op_params["pulse_time"]; op_params["filter_params"]...)
    end

    Lₜ = Get_Lₜ(Ô_D, drive_coef)

    if (length(tspan) == 0) & (spns != "Stroboscopic")
        tspan = collect(LinRange(0, op_params["pulse_time"], Int(ceil(op_params["pulse_time"]*spns))+1))
    end
    if spns == "Stroboscopic"
        tspan = Get_Stroboscopic_Times(op_params)[1:strob_skip:end]
    end
    
    solver_kwargs_sym = Dict{Symbol, Any}()
    for key in keys(solver_kwargs)
        solver_kwargs_sym[Symbol(key)] = solver_kwargs[key]
    end
    @info "Running Time Evolution"
    sleep(0.01)
    res = qt.mesolve(2*π*Ĥ, ρ, tspan, c_ops, H_t = Lₜ; solver_kwargs_sym...)
    @info "Time Evolution Complete"
    sleep(0.01)

    if (save_step) | (to_return == "DS")
        @info "Saving Steps"
        sleep(0.01)
        
        properties = Dict{Any, Any}("Data"=>"Density Matrices", "Time"=>string(now()), "Times" => res.times, "Operator" => op_name, "Operator_Parameters" => string(op_params))
        
        ds_properties = Dict{Any, Any}("dims" => collect(Ĥ.dims), "Data"=>"Density Matrices", "Solver_Args" => string(solver_kwargs))
        ds_properties = merge(ds_properties, other_ds_properties)
        
        ds = Utils.Qobj_List_To_DS(res.states; cube_name = Symbol(step_name), cube_properties = properties, ds_properties = ds_properties, step_name = Symbol(step_name*"_steps"))

        if save_step
            file_name = save_path*run_name
            if save_as_seperate_file
                file_name = file_name*step_name
            end

            if !isdir(save_path)
                @info "Making Save Directory $save_path"
                mkdir(save_path)
            end
            savedataset(ds,path = file_name*".nc", append = true)
        end
    end
    
    @info "Done with $step_name"
    if to_return == "All"
        return res
    elseif to_return == "Last"
        return res.states[end]
    end
end

function RunPulseSequence(Ĥ::qt.QuantumObject, Ô_D::qt.QuantumObject, 
    ψ::qt.QuantumObject{<:AbstractVector{T1}, qt.KetQuantumObject}, 
    op_sequence,
    op_params_dict; 
    spns = 5, 
    solver_kwargs = Dict{Any, Any}(), 
    run_name = "", 
    save_path = "Data/",
    Return = true,
    clean_up = true,
    strob_skip = 1,
    ) where T1<:Number
    #-------------------------------------------------------------
    df = DateFormat("e-u-d-yy.H.M")
    t = now()
    the_time = string(Dates.format(t, df))
    if run_name == ""
        run_name = "Operator_Sequence_"*the_time
    end
    println("The Name for this run is: "*run_name)
    println("It is being saved at: "*save_path)
    cube_order = ""
    for i in 1:length(op_sequence)
        step_name = "Step_"*string(i)
        cube_order = cube_order*" "*step_name
    end
    other_ds_properties = Dict{Any, Any}("Order" => cube_order)
    for i in 1:length(op_sequence)
        op = op_sequence[i]
        @info "Running operator $op"
        step_name = "Step_"*string(i)

        op_params = Dict{Any, Any}()
        if op[1:5] == "wait_"
            wait_time = parse(Float64, op[6:end])
            op_params = Dict("pulse_time" => wait_time, "freq_d" => 0, "shift" => 0, "epsilon" => 0, "Envelope" => "Square", "Envelope Args" => Dict{Any, Any}())
        else
            op_params = op_params_dict[op]
        end
        ψ = RunSingleOperator(Ĥ,Ô_D, ψ, op_params; spns = spns, solver_kwargs = solver_kwargs, run_name = run_name, save_path = save_path, step_name = step_name, to_return = "Last", save_step = true, op_name = op, other_ds_properties = other_ds_properties, strob_skip = strob_skip)
    end
    @info "Done With Running Sequence"
    if Return
        @info "Loading Data"
        dat = Utils.LoadRunResults(save_path*run_name*".nc")
        if clean_up
            rm(save_path*run_name*".nc")
        end
        return dat
    end
end

function RunPulseSequence(Ĥ::qt.QuantumObject, Ô_D::qt.QuantumObject, 
    ρ::qt.QuantumObject{<:AbstractArray{T1},qt.OperatorQuantumObject}, 
    op_sequence,
    op_params_dict;
    spns = 5, 
    solver_kwargs = Dict{Any, Any}(), 
    run_name = "", 
    save_path = "Data/", 
    c_ops = [],
    other_ds_properties = Dict{Any, Any}(),
    strob_skip = 1
    ) where T1<:Number
    #-------------------------------------------------------------
    df = DateFormat("e-u-d-yy.H.M")
    t = now()
    the_time = string(Dates.format(t, df))
    if run_name == ""
        run_name = "Operator_Sequence_"*the_time
    end
    println("The Name for this run is: "*run_name)
    println("It is being saved at: "*save_path)
    cube_order = ""
    for i in 1:length(op_sequence)
        step_name = "Step_"*string(i)
        cube_order = cube_order*" "*step_name
    end
    
    other_ds_properties["Order"] = cube_order

    num_steps = length(op_sequence)

    for i in 1:length(op_sequence)
        op = op_sequence[i]
        op_params = Dict{Any, Any}()
        if op[1:5] == "wait_"
            wait_time = parse(Float64, op[6:end])
            op_params = Dict("pulse_time" => wait_time, "freq_d" => 0, "shift" => 0, "epsilon" => 0, "Envelope" => "Square", "Envelope Args" => Dict{Any, Any}())
        else
            op_params = op_params_dict[op]
        end
        @info "Step $i/$num_steps: Running operator $op"
        step_name = "Step_"*string(i)
        ρ = RunSingleOperator(Ĥ,Ô_D, ρ, op_params; spns = spns, solver_kwargs = solver_kwargs, run_name = run_name, save_path = save_path, step_name = step_name, to_return = "Last", c_ops = c_ops, save_step = true, op_name = op, other_ds_properties = other_ds_properties, strob_skip = strob_skip)
    end
    @info "Done With Running Sequence"
end



