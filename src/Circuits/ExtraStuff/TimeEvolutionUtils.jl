
export Get_Drive_Hamiltonian, Get_Drive_Hamiltonian_With_Envelope, f_for_schroedinger_dynamic, f_for_schroedinger_dynamic, f_for_master_dynamic

function Get_Drive_Hamiltonian(model, op, ν, ε)
    return Get_Drive_Hamiltonia_With_Envelope(model, op, ν, ε, Square_Envelope)
end
function Get_Drive_Hamiltonian_With_Envelope(model, op, ν, ε, envelope)
    return t->ε*envelope(t)*sin(2π*ν*t)*op+model.hilbertspace.Ĥ
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

function f_for_schroedinger_dynamic(model, op, ν, ε; envelope_name = "Square Envelope", envelope_kwargs = Dict{Any, Any}())
    
    envelope = Get_Envelope(envelope_name, envelope_kwargs)

    op = qo.sparse(op)
    Ĥ = qo.sparse(model.Ĥ)

    return (t, ψ) -> 2*π*(ε*envelope(t)*sin(2*π*ν*t)*op+Ĥ)
end

function f_for_master_dynamic(model, op, ν, ε; c_ops = [], envelope_name = "Square Envelope", envelope_kwargs = Dict{Any, Any}())
    
    envelope = Get_Envelope(envelope_name, envelope_kwargs)

    op = qo.sparse(op)
    Ĥ = qo.sparse(model.Ĥ)

    c_ops = collect(values(c_ops))

    return (t, ψ) -> [2*π*(ε*envelope(t)*sin(2*π*ν*t)*op+Ĥ), c_ops, qo.dagger.(c_ops)]
end