include("ExtraStuff/Envelopes.jl")

export Get_Drive_Hamiltonian, Get_Drive_Hamiltonian_With_Envelope, f_for_schroedinger_dynamic, f_for_schroedinger_dynamic, f_for_master_dynamic

function Get_Drive_Hamiltonian(model, op, ν, ε)
    return Get_Drive_Hamiltonian_With_Envelope(model, op, ν, ε, Square_Envelope)
end

function Get_Drive_Hamiltonian_With_Envelope(model, op, ν, ε, envelope)
    return t->2*pi*(ε*envelope(t)*sin(2π*ν*t)*op)
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
