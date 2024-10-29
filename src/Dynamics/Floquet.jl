import QuantumToolbox as qt
import OrdinaryDiffEq as ODE
using ProgressMeter

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


function Get_Pulse_Floquet_Basis(hilbertspace::Hilbertspaces.Hilbertspace, Ĥ_D, op_params)
    ε0 = op_params["epsilon"]
    εt = ε0*Envelopes.Envelope_Dict[op_params["Envelope"]](t; op_params["Envelope Args"]...)
end