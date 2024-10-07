import QuantumToolbox as qt
import ProgressMeter as PM
using OrdinaryDiffEqVerner: Vern9

export Get_Drive_Coef, Get_Evelope, Get_Ĥ_D, Get_Lₜ, Propagator

function Get_Drive_Coef(ν::T1,
    ε;
    envelope = Envelopes.Square_Envelope
    )where T1<:Number

    function drive_coef(t, params...)
        return 2*π*ε*envelope(t)*sin(2π*ν*t)
    end
    return drive_coef
end

function Get_Drive_Coef(ν::Function,
    ε;
    envelope = Envelopes.Square_Envelope)

    function drive_coef(t, params...)
        εt = ε*envelope(t)
        return 2*π*εt*sin(2π*ν(εt)*t)
    end
    return drive_coef
end


function Get_Ĥ_D(op::qt.QuantumObject, drive_coef::Union{Nothing, Function}; TDOS = true, params = nothing, init_time = 0.0)
    drive_coef = (drive_coef isa Nothing) ? (t)->1 : drive_coef

    if TDOS
        return qt.TimeDependentOperatorSum([drive_coef], [op], params = params, init_time = init_time)
    else
        return t -> drive_coef(t)*op
    end
end

function Get_Lₜ(op, drive_coef; params = nothing, init_time = 0.0)
    return qt.TimeDependentOperatorSum([drive_coef], [qt.liouvillian(op)], params = params, init_time = init_time)
end


function Propagator(hilbertspace, Ĥₜ, tf; progress_meter = false, ti = 0)
    U = 0*Utils.eye_like(hilbertspace.Ĥ)

    p = PM.Progress(length(hilbertspace.dressed_states), enabled = progress_meter)
    for state in keys(hilbertspace.dressed_states)
        ψi = hilbertspace.dressed_states[state]
        se_res = qt.sesolve(2*π*hilbertspace.Ĥ, ψi, [ti, tf], H_t = Ĥₜ, progress_bar = false, alg = Vern9())
        ψf = se_res.states[end]
        U += ψf*ψi'
        PM.next!(p)
    end
    return U
end