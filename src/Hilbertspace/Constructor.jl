@kwdef struct Hilbertspace
    Components :: Dict
    Interactions :: Vector
    𝕀̂_Dict :: Dict
    𝕀̂ :: qt.QuantumObject
    Ĥ :: qt.QuantumObject

    dressed_states :: Dict
    dressed_energies :: Dict
end

function init(Components, Interactions; order = [])
    if length(order) == length(Components)
        key_list = order
    else
        key_list = collect(keys(Components))
    end
    
    𝕀̂_Dict = Dict{Any, Any}()
    for key in key_list
        𝕀̂_Dict[key] = qt.eye(Components[key].dim)
    end
    Ĥ_comp_vec = []
    for key in key_list
        push!(Ĥ_comp_vec, Components[key].Ĥ)
    end

    Ĥ_non_int_list = []
    for key in key_list
        op_dict = Dict(Components[key].name => Components[key].Ĥ)
        push!(Ĥ_non_int_list, Utils.IdentityWrapper(𝕀̂_Dict, op_dict, order = order))
    end
    
    Ĥ_non_int = sum(Ĥ_non_int_list)
    Ĥ_int = 0*Ĥ_non_int
    for key in keys(Interactions)
        term = Utils.IdentityWrapper(𝕀̂_Dict, Interactions[key]["ops"], order = order)*Interactions[key]["g"]
        Ĥ_int += term
    end

    Ĥ = Ĥ_non_int+Ĥ_int
    λ_dressed, ψ_dressed = qt.eigenstates(Ĥ)

    dims = []
    for key in key_list
        push!(dims, Components[key].dim)
    end

    for_iter = []
    for i in 1:length(dims)
        push!(for_iter, collect(1:dims[i]))
    end
    states_to_iter = Iterators.product(for_iter...)

    dressed_states = Dict{Any, Any}()
    dressed_energies = Dict{Any, Any}()

    for state in states_to_iter
        overlaps = zeros(length(ψ_dressed))
        bare_ψ_list = []
        for i in 1:length(key_list)
            key = key_list[i]
            push!(bare_ψ_list, qt.Qobj(Components[key].eigsys.vectors[:, state[i]]))
        end
        ψ_bare = qt.tensor(bare_ψ_list...)

        for i in 1:length(overlaps)
            overlaps[i] = norm((ψ_dressed[i]'*ψ_bare))^2
        end
        max_idx = argmax(overlaps)
        dressed_states[state.-1] = ψ_dressed[max_idx]
        dressed_energies[state.-1] = Real(λ_dressed[max_idx])
    end

    𝕀̂_vec = []
    for i in key_list
        push!(𝕀̂_vec, 𝕀̂_Dict[i])
    end

    if length(𝕀̂_vec) == 1
        𝕀̂ = 𝕀̂_vec[1]
    else
        𝕀̂ = qt.tensor(𝕀̂_vec...)
    end
    return Hilbertspace(Components=Components, Interactions = Interactions, 𝕀̂_Dict =𝕀̂_Dict, Ĥ = Ĥ, dressed_states = dressed_states, dressed_energies = dressed_energies, 𝕀̂ = 𝕀̂)

end