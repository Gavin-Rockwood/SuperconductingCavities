import QuantumToolbox as qt
import OrdinaryDiffEq as ODE

export Get_Floquet_t0_eigsys, Floquet_0_Sweep

function Get_Floquet_t0_eigsys(hilbertspace::Hilbertspaces.Hilbertspace, Ĥ_D, T)
    
    U = Propagator(hilbertspace, Ĥ_D, T)

    λs, λ⃗s = qt.eigenstates(U)
    λs = -angle.(λs)/T#imag(log.(λs))
    return λs, λ⃗s
end

function Floquet_0_Sweep(hilbertspace::Hilbertspaces.Hilbertspace, drive_op, list_of_params; use_logging=true)
    STEPS = length(list_of_params)

    F_Modes = []
    F_Energies = []

    if (use_logging) @info "Beginning Floquet Sweep" end
    for i in 1:STEPS
        if (use_logging) @debug "On Param Set Number $i" end
        ν = list_of_params[i]["ν"]
        ε = list_of_params[i]["ε"]
        drive_coef = Get_Drive_Coef(ν, ε)
        Ĥ_D = Get_Ĥ_D(drive_op, drive_coef)
        λs, λ⃗s = Get_Floquet_t0_eigsys(hilbertspace, Ĥ_D, 1/ν)

        push!(F_Energies, λs)
        push!(F_Modes, λ⃗s)
    end

    res = Dict{Any, Any}()
    res["F_Modes"] = F_Modes
    res["F_Energies"] = F_Energies

    if (use_logging) @info "Done With Floquet Sweep" end
    return res
end