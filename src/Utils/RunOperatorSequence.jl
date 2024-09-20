import QuantumToolbox as qt
using YAXArrays
using Dates
using NetCDF
import OrdinaryDiffEq as ODE
import ProgressMeter as PM
using OrdinaryDiffEqVerner: Vern9

function RunSingleOperator(model::Transmon_Resonators, 
    ψ::qt.QuantumObject{<:AbstractVector{T1},qt.KetQuantumObject}, 
    op_params; 
    spps = 5, 
    solver_kwargs = Dict{Any, Any}(), 
    save_step = false, 
    step_name = "DEFAULT", 
    op_name = "DEFAULT", 
    to_return = "All", 
    save_path = "Data/", 
    run_name = "", 
    save_as_seperate_file = false, 
    tspan = [],
    other_ds_properties = Dict{Any, Any}()
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

    if (save_step) | (to_return == "DS")
        @info "Saving Steps"
        sleep(0.01)
        
        properties = Dict{Any, Any}("Data"=>"Wave Functions", "Time"=>string(now()), "Times" => res.times, "Operator" => op_name)
        
        ds_properties = Dict{Any, Any}("dims" => collect(model.hilbertspace.Ĥ.dims), "Data"=>"Wave Functions", "Solver_Args" => string(solver_kwargs))
        for key in keys(other_ds_properties)
            ds_properties[key] = other_ds_properties[key]
        end
        
        ds = Qobj_List_To_DS(res.states; cube_name = Symbol(step_name), cube_properties = properties, ds_properties = ds_properties, step_name = Symbol(step_name*"_steps"))
        
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

function RunSingleOperator(model::Transmon_Resonators, 
    ρ::qt.QuantumObject{<:AbstractArray{T1},qt.OperatorQuantumObject}, 
    op_params; 
    c_ops = nothing,  
    spps = 5, 
    solver_kwargs = Dict{Any, Any}(), 
    save_step = false, 
    step_name = "DEFAULT", 
    op_name = "DEFAULT",
    to_return = "All", 
    save_path = "Data/", 
    run_name = "", 
    save_as_seperate_file = false, 
    tspan = [],
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

    if (save_step) | (to_return == "DS")
        @info "Saving Steps"
        sleep(0.01)
        
        properties = Dict{Any, Any}("Data"=>"Density Matrices", "Time"=>string(now()), "Times" => res.times, "Operator" => op_name)
        
        ds_properties = Dict{Any, Any}("dims" => collect(model.hilbertspace.Ĥ.dims), "Data"=>"Density Matrices", "Solver_Args" => string(solver_kwargs))
        
        ds = Qobj_List_To_DS(res.states; cube_name = Symbol(step_name), cube_properties = properties, ds_properties = ds_properties, step_name = Symbol(step_name*"_steps"))
        
        for key in keys(other_ds_properties)
            ds_properties[key] = other_ds_properties[key]
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

function RunPulseSequence(model::Transmon_Resonators, 
    ψ::qt.QuantumObject{<:AbstractVector{T1}, qt.KetQuantumObject, 2}, 
    op_sequence; 
    spps = 5, 
    solver_kwargs = Dict{Any, Any}(), 
    run_name = "", 
    save_path = "Data/",
    Return = false,
    clean_up = false
    ) where T1<:Number
    #-------------------------------------------------------------
    if run_name == ""
        run_name = "Operator_Sequence_"*string(now())
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
        ψ = RunSingleOperator(model, ψ, model.Stuff["op_drive_params"][op]; spps = spps, solver_kwargs = solver_kwargs, run_name = run_name, save_path = save_path, step_name = step_name, to_return = "Last", save_step = true, op_name = op, other_ds_properties = other_ds_properties)
    end

    if Return
        dat = LoadRunResults(save_path*run_name*".nc")
        if clean_up
            rm(save_path*run_name*".nc")
        end
        return dat
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
        ρ = RunSingleOperator(model, ρ, model.Stuff["op_drive_params"][op]; spps = spps, solver_kwargs = solver_kwargs, run_name = run_name, save_path = save_path, step_name = step_name, to_return = "Last", c_ops = c_ops, save_step = true, op_name = op, other_ds_properties = other_ds_properties)
    end
end



