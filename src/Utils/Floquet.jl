import QuantumOptics as qo

export Get_Floquet_t0_eigsys, Floquet_0_Sweep

function Get_Floquet_t0_eigsys(Ĥₜ, T; N_Steps = 100)
    δt = T/N_Steps

    FloqOp = qo.dense(qo.identityoperator(qo.basis(Ĥₜ(0))))

    t = 0
    for i in 1:N_Steps
        t += δt
        FloqOp = qo.exp(qo.dense(-1im*Ĥₜ(t)*δt))*FloqOp
    end

    λs, λ⃗s = qo.eigenstates(qo.dense(FloqOp), warning = false)
    λs = imag(log.(λs))
    return λs, λ⃗s
end

function Floquet_0_Sweep(model, drive_op, list_of_params; Floq_N_Steps = 100)
    STEPS = length(list_of_params)

    F_Modes = []
    F_Energies = []

    @info "Beginning Floquet Sweep"
    for i in 1:STEPS
        @debug "On Param Set Number $i"
        Ĥₜ = Get_Drive_Hamiltonian(model, drive_op, list_of_params[i]["ν"], list_of_params[i]["ε"])
        λs, λ⃗s = Get_Floquet_t0_eigsys(Ĥₜ, 1/list_of_params[i]["ν"], N_Steps = Floq_N_Steps)

        push!(F_Energies, λs)
        push!(F_Modes, λ⃗s)
    end

    res = Dict{Any, Any}()
    res["F_Modes"] = F_Modes
    res["F_Energies"] = F_Energies

    @info "Done With Floquet Sweep"
    return res
end