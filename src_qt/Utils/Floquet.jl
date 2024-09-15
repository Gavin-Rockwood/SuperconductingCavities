import QuantumOptics as qo
import OrdinaryDiffEq as ODE

export Get_Floquet_t0_eigsys, Floquet_0_Sweep

function Get_Floquet_t0_eigsys(f, T)
    
    FloqOp = qo.dense(qo.identityoperator(qo.basis(f(0,0))))
    
    tspan = collect(LinRange(0, T, 2))

    res = qo.timeevolution.schroedinger_dynamic(tspan, FloqOp, f, alg = ODE.Vern9())

    λs, λ⃗s = qo.eigenstates(res[2][end], warning = false)
    λs = -angle.(λs)/T#imag(log.(λs))
    return λs, λ⃗s
end

function Floquet_0_Sweep(model, drive_op, list_of_params)
    STEPS = length(list_of_params)

    F_Modes = []
    F_Energies = []

    @info "Beginning Floquet Sweep"
    for i in 1:STEPS
        @debug "On Param Set Number $i"
        f = f_for_schroedinger_dynamic(model, drive_op, list_of_params[i]["ν"], list_of_params[i]["ε"])
        λs, λ⃗s = Get_Floquet_t0_eigsys(f, 1/list_of_params[i]["ν"])

        push!(F_Energies, λs)
        push!(F_Modes, λ⃗s)
    end

    res = Dict{Any, Any}()
    res["F_Modes"] = F_Modes
    res["F_Energies"] = F_Energies

    @info "Done With Floquet Sweep"
    return res
end