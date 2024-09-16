import QuantumToolbox as qt
include("ExtraStuff/Envelopes.jl")

export Get_Drive_Hamiltonian, Get_Drive_Hamiltonian_With_Envelope, f_for_schroedinger_dynamic, f_for_schroedinger_dynamic, f_for_master_dynamic

function Get_Drive_Coef(ν, ε, envelope = Square_Envelope)
    return t->2*pi*(ε*envelope(t)*sin(2π*ν*t))
end

function Get_Envelope(envelope_name, envelope_kwargs)
    envelope_kwargs_sym = Dict{Symbol, Any}()
    for key in keys(envelope_kwargs)
        envelope_kwargs_sym[Symbol(key)] = envelope_kwargs[key]
    end
    function envelope(t)
        return Envelope_Dict[envelope_name](t; envelope_kwargs_sym...)
    end
    return envelope
end

function Get_Ĥₜ(op::qt.QuantumObject, drive_coef::Union{Nothing, Function}; TDOS = true, params = nothing, init_time = 0.0)
    drive_coef = (drive_coef isa Nothing) ? (t)->1 : drive_coef

    if TDOS
        return qt.TimeDependentOperatorSum([cf], [op], params = params, init_time = init_time)
    else
        return t -> drive_coef(t)*op
    end
end

function Get_Lₜ(op, drive_coef; params = nothing, init_time = 0.0)
    return qt.TimeDependentOperatorSum([drive_coef], [qt.liouvillian(op)], params = params, init_time = init_time)
end
