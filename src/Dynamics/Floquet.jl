import QuantumToolbox as qt
import OrdinaryDiffEq as ODE
import LsqFit as LF
using ProgressMeter
using DimensionalData
import Peaks

#export Get_Floquet_eigsys, Floquet_0_Sweep

function Get_Floquet_t0_Eigsys(hilbertspace::Hilbertspaces.Hilbertspace, Ĥ_D, T; t0 = 0)
    
    U = Propagator(hilbertspace, Ĥ_D, T+t0, ti = t0)

    λs, λ⃗s = qt.eigenstates(U)
    λs = -angle.(λs)/T#imag(log.(λs))
    return λs, λ⃗s
end

function Floquet_t0_Sweep(hilbertspace::Hilbertspaces.Hilbertspace,
    drive_op,
    list_of_params;
    use_logging=true,
    states_to_track::T1 =  Dict{Any, Any}()
    )where T1<:Dict
    STEPS = length(list_of_params)

    F_Modes = []
    F_Energies = []

    if (use_logging) @info "Beginning Floquet Sweep" end
    P = Progress(STEPS)
    for i in 1:STEPS

        # Check if we have already done this parameter set. This way I can just reuse the values instead of recalculating. 
        # Note that this relies on dictionary equality, so it could get a bit confused if the dictionary was defined in different orders 
        # but have the same parameters
        checking_if_done = findall(x->x == list_of_params[i], list_of_params[1:i-1])
        if length(checking_if_done) >0
            idx = checking_if_done[1]
            λs = F_Energies[idx]
            λ⃗s = F_Modes[idx]
            push!(F_Energies, λs)
            push!(F_Modes, λ⃗s)
            next!(P)
        else
            if (use_logging) @debug "On Param Set Number $i" end
            ν = list_of_params[i]["ν"]
            ε = list_of_params[i]["ε"]

            if "t0" in keys(list_of_params[i])
                t0 = list_of_params[i]["t0"]
            else    
                t0 = 0
            end
            drive_coef = Get_Drive_Coef(ν, ε)
            Ĥ_D = Get_Ĥ_D(drive_op, drive_coef)
            λs, λ⃗s = Get_Floquet_t0_Eigsys(hilbertspace, Ĥ_D, 1/ν, t0 = t0)
            if ("t" in keys(list_of_params[i])) & (length(states_to_track) == 0)
                for j in 1:length(λ⃗s)
                    op_params = Dict{Any, Any}()
                    op_params["epsilon"] = ε
                    op_params["freq_d"] = ν
                    op_params["shift"] = 0
                    T = (abs(1/op_params["freq_d"]))
                    op_params["pulse_time"] = list_of_params[i]["t"]#%T
                    n = floor(list_of_params[i]["t"]/(abs(1/op_params["freq_d"])))*T
                    
                    if "Envelope" in keys(list_of_params[i])
                        op_params["Envelope"] = list_of_params[i]["Envelope"]
                    else
                        op_params["Envelope"] = "Square"
                    end
                    if "Envelope Args" in keys(list_of_params[i])
                        op_params["Envelope Args"] = list_of_params[i]["Envelope Args"]
                    else
                        op_params["Envelope Args"] = Dict{Any, Any}()
                    end
                    λ⃗s[j] = RunSingleOperator(hilbertspace.Ĥ, drive_op, λ⃗s[j], op_params; to_return = "Last", progress_bar = false, use_logging = false)#*ℯ^(-1im*λs[j]*n*T)
                end
            end

            push!(F_Energies, λs)
            push!(F_Modes, λ⃗s)
            next!(P)
        end
    end

    res = Dict{Any, Any}()
    res["F_Modes"] = F_Modes
    res["F_Energies"] = F_Energies

    if (use_logging) @info "Done With Floquet Sweep" end
    
    if length(states_to_track) > 0
        if (use_logging) @info "Tracking State" end
        other_sorts = Dict("Quasienergies" => F_Energies)
        tracking_res = Utils.State_Tracker(F_Modes, states_to_track, other_sorts = other_sorts, use_logging = use_logging)
        res["Tracking"] = tracking_res


        if (use_logging) @info "Running the necessary time evolutions" end
        P = Progress(STEPS)
        for step in 1:STEPS   
            if ("t" in keys(list_of_params[step]))
                op_params = Dict{Any, Any}()
                op_params["epsilon"] = list_of_params[step]["ε"]
                op_params["freq_d"] = list_of_params[step]["ν"]
                op_params["shift"] = 0
                T = (abs(1/op_params["freq_d"]))
                n = floor(list_of_params[step]["t"]/(abs(1/op_params["freq_d"])))*T
                op_params["pulse_time"] = list_of_params[step]["t"]%T
                op_params["Envelope"] = "Square"
                op_params["Envelope Args"] = Dict{Any, Any}()
                for j in 1:length(states_to_track) 
                    state_key = string(collect(keys(states_to_track))[j])
                    ψ = tracking_res[State = At(state_key), Step = At(step)]["ψ"]
                    
                    qe = tracking_res[State = At(state_key), Step = At(step)]["Quasienergies"]
                    ψ = RunSingleOperator(hilbertspace.Ĥ, drive_op, ψ, op_params; to_return = "Last", progress_bar = false, use_logging = false)
                    tracking_res[State = At(state_key), Step = At(step)]["ψ"] = ℯ^(-1im*qe*n*T)*ψ
                end
            end
            next!(P)
        end

        return tracking_res
    end


    return res
end


function Get_Floquet_t0_Table(hilbertspace::Hilbertspaces.Hilbertspace, Ĥ_D, T, states_to_track; use_logging = false)
    λs, λ⃗s = Get_Floquet_t0_Eigsys(hilbertspace, Ĥ_D, T)

    other_sorts = Dict("Quasienergies" => [λs])
    tracking_res = Utils.State_Tracker([λ⃗s], states_to_track, other_sorts = other_sorts, use_logging = use_logging)

    return tracking_res
end


function Get_Stroboscopic_Times(
    ϕ::T1,
    pulse_time;
    sample_frequency = 10000,
    sample_every = 1,
    start_at = 2,
    return_all = false,
    recieving_drive_coef = false
    ) where T1<:Union{Function, ODE.ODESolution}

    times_to_sample = collect(LinRange(0, pulse_time, ceil(Int, sample_frequency*pulse_time)))

    
    f_for_peak_dc(t) = -abs.(ϕ(t))
    f_for_peak(t) = -abs.(sin.(ϕ(t)))
    

    if recieving_drive_coef
        for_peak = f_for_peak_dc.(times_to_sample)
    else
        for_peak = f_for_peak.(times_to_sample)
    end
    indices, heighs = Peaks.findmaxima(for_peak)

    times = vcat(times_to_sample[indices][start_at:2*sample_every:end], [pulse_time])

    if return_all
        return [times, indices, heighs]
    else
        return times
    end
end

function Get_Stroboscopic_Times(
    pulse::T1;kwargs...
    ) where T1<:Dict

    ε0 = pulse["epsilon"]
    ν0 = pulse["freq_d"]
    νε = ν0+pulse["shift"]
    
    digitize = false
    step_length = 2.3
    if "digitize" in keys(pulse)
        digitize = pulse["digitize"]
        if "step_length" in keys(pulse)
            step_length = pulse["step_length"]
        end
    end
    envelope = Envelopes.Get_Envelope(pulse["Envelope"], pulse["Envelope Args"], digitize = digitize, step_length = step_length)
    if "chirp_params" in keys(pulse)
        νε=chirper(ν0, pulse["chirp_params"])
    end

    filtered = false
    drive_coef = Get_Drive_Coef(νε, ε0, envelope = envelope, drive_time = pulse["pulse_time"]) 
    if "filter_params" in keys(pulse)
        filtered = true
        filter_params = Dict(Symbol(key)=>val for (key, val) in pulse["filter_params"]) # turns the string keys into symbols. 
        drive_coef = Get_Low_Pass_Filtered_Drive_Coef(drive_coef, pulse["pulse_time"]; filter_params...)
    end

    if filtered
        recieving_drive_coef = true
        ϕ(t) = drive_coef(t)
    else
        ϕ = Get_Drive_Coef(νε, ε0, envelope = envelope, drive_time = pulse["pulse_time"], return_ϕ = true)
        recieving_drive_coef = false
    end
    Get_Stroboscopic_Times(ϕ, pulse["pulse_time"]; recieving_drive_coef = recieving_drive_coef, kwargs...)
end

function Get_Pulse_Floquet_Sweep(hilbertspace::Hilbertspaces.Hilbertspace,
    Ĥ_D,
    pulse_params;
    stroboscopic_times = [],
    sample_frequency = 100,
    states_to_track = Dict{Any, Any}()
    )
    drive_time = pulse_params["pulse_time"]
    ε0 = pulse_params["epsilon"]

    digitize = false
    step_length = 2.3
    if "digitize" in keys(pulse_params)
        digitize = pulse_params["digitize"]
        if "step_length" in keys(pulse_params)
            step_length = pulse_params["step_length"]
        end
    end

    envelope = Envelopes.Get_Envelope(pulse_params["Envelope"], pulse_params["Envelope Args"], digitize = digitize, step_length = step_length)

    εt(t) = ε0*envelope(t)
    ν0 = pulse_params["freq_d"]
    νε(ε) = ν0+pulse_params["shift"]

    drive_coef = Get_Drive_Coef(νε, ε0, envelope = envelope, drive_time = drive_time) 

    filtered = false
    if "filter_params" in keys(pulse_params)
        filtered = true
        filter_params = Dict(Symbol(key)=>val for (key, val) in pulse_params["filter_params"]) # turns the string keys into symbols. 
        println("Filtering: $filtered")
        drive_coef = Get_Low_Pass_Filtered_Drive_Coef(drive_coef, pulse_params["pulse_time"]; filter_params...)
    end

    if "chirp_params" in keys(pulse_params)
        νε=chirper(ν0, pulse_params["chirp_params"])
    end

    if length(stroboscopic_times) == 0
        ϕ = Get_Drive_Coef(νε, ε0, envelope=envelope, return_ϕ = true, drive_time = drive_time)
        stroboscopic_times = Get_Stroboscopic_Times(ϕ, pulse_params["pulse_time"]; sample_frequency = sample_frequency)
    end

    list_of_params = []
    rt = pulse_params["Envelope Args"]["ramp_time"]
    pt = pulse_params["pulse_time"]
    for i in 1:length(stroboscopic_times)
        t = stroboscopic_times[i]
        if filtered
            ti = t
            if i>1
                ti = stroboscopic_times[i-1]
            end
            tf = t
            if i<length(stroboscopic_times)
                tf = stroboscopic_times[i+1]
            end
            times_to_sample = collect(ti:1/sample_frequency:tf)

            signal = drive_coef.(times_to_sample)
            
            to_fit(t, p) = 2*π.*p[1].*εt.(t).*sin.(2*π.*p[2].*t.+p[3])
            #to_fit(t,p) = 2*π.*εt.(t).*p[1].*sin.(2*π.*p[2].*t.+p[3])
            fit = LF.curve_fit(to_fit, times_to_sample, signal, [maximum(signal), νε(εt(2*rt)), 0.0])
            #println("Fit Params: "*Utils.tostr(fit.param))
            push!(list_of_params, Dict{Any,Any}("ν" => abs(fit.param[2]), "ε" => εt(t)*fit.param[1]))#/(2*π))) # This is the wrong frequency to use. Maybe need a phase shift???? 

        else
            push!(list_of_params, Dict{Any,Any}("ν" => νε(εt(t)), "ε" => εt(t)))
        end
    end

    if length(states_to_track) == 0
        states_to_track = hilbertspace.dressed_states
    end
    floq_sweep_res = Floquet_t0_Sweep(hilbertspace, Ĥ_D, list_of_params; states_to_track = states_to_track)

    return floq_sweep_res
end


function Pulse_Floquet_Projections(
    pulse_res,
    floq_sweep::T2
    ) where T2<:DimMatrix


    states = collect(floq_sweep.dims[findall(x -> x == :State, name(floq_sweep.dims))[1]])

    projections = Dict{Any, Any}()
    for state in states
        projections[state] = []
        for i in 1:length(pulse_res)
            Φ = floq_sweep[State = At(state), Step = At(i)]["ψ"]
            push!(projections[state], abs(Φ'*pulse_res[i])^2)
        end
    end
    return projections
end


"""
    This reformats the floquet sweep results from an YAXARRAY into a set of nested dictionaries that can be saved as a .json file.
"""
function Reformat_Sweep_Results_To_Save(floq_sweep)
    dims = floq_sweep.dims;
    tracked_items = collect(keys(floq_sweep[1]));
    tracked_items[findall(x->x=="ψ", tracked_items)[1]] = "Quasimodes"

    dim_names = Dimensions.label.(dims)
    floq_sweep_new = Dict{Any, Any}()

    for state in dims[1].val
        floq_sweep_new[state] = Dict{Any, Any}()
        for step in dims[2].val
            floq_sweep_new[state][step] = Dict{Any, Any}()
            for item in tracked_items
                item_name = item
                if item_name == "Quasimodes"
                    item_name = "ψ"
                end
                dat = floq_sweep[State = At(state), Step = At(step)][item_name]
                if typeof(dat) <: qt.QuantumObject
                    dat = dat.data
                end
                floq_sweep_new[state][step][item] = dat
            end
        end
    end
    return floq_sweep_new
end