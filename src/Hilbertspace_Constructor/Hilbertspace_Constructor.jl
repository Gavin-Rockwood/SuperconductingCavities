import QuantumOptics as qo
using LinearAlgebra
#using ProtoStructs

#path_to_utils = join(split(@__DIR__, "/")[1:end-1], "/")*"/Utils"
#include(path_to_utils*"/HilbertSpaceUtils.jl")
include("Elements/Transmon.jl")
include("Elements/Resonator.jl")
include("Hilbertspace_Struct.jl")

export Hilbertspace_Constructor


function Hilbertspace_Constructor(Components, Interactions; order = [])
    if length(order) == length(Components)
        key_list = order
    else
        key_list = collect(keys(Components))
    end
    
    𝕀̂_Dict = Dict{Any, Any}()
    for key in key_list
        𝕀̂_Dict[key] = qo.identityoperator(Components[key].𝔹)
    end
    Ĥ_comp_vec = []
    for key in key_list
        push!(Ĥ_comp_vec, Components[key].Ĥ)
    end

    Ĥ_non_int_list = []
    for key in key_list
        op_dict = Dict(Components[key].name => Components[key].Ĥ)
        push!(Ĥ_non_int_list, IdentityWrapper(𝕀̂_Dict, op_dict, order = order))
    end
    
    Ĥ_non_int = sum(Ĥ_non_int_list)
    Ĥ_int = 0*Ĥ_non_int
    for key in keys(Interactions)
        term = IdentityWrapper(𝕀̂_Dict, Interactions[key]["ops"], order = order)*Interactions[key]["g"]
        Ĥ_int += term
    end

    Ĥ = Ĥ_non_int+Ĥ_int

    dressed_eigsys = qo.eigenstates(qo.dense(Ĥ))

    dims = []
    for key in key_list
        push!(dims, size(Components[key].Ĥ)[1])
    end

    for_iter = []
    for i in 1:length(dims)
        push!(for_iter, collect(1:dims[i]))
    end
    states_to_iter = Iterators.product(for_iter...)

    dressed_states = Dict{Any, Any}()
    dressed_energies = Dict{Any, Any}()

    for state in states_to_iter
        overlaps = zeros(length(dressed_eigsys[1]))
        bare_ψ_list = []
        for i in 1:length(key_list)
            key = key_list[i]
            push!(bare_ψ_list, Components[key].eigsys[2][state[i]])
        end
        ψ_bare = qo.tensor(bare_ψ_list...)

        for i in 1:length(overlaps)
            overlaps[i] = norm((dressed_eigsys[2][i]'*ψ_bare))^2
        end
        max_idx = argmax(overlaps)
        dressed_states[state.-1] = dressed_eigsys[2][max_idx]
        dressed_energies[state.-1] = dressed_eigsys[1][max_idx]
    end


    return HilbertSpace(Components=Components, Interactions = Interactions, 𝕀̂_Dict =𝕀̂_Dict, Ĥ = Ĥ, dressed_states = dressed_states, dressed_energies = dressed_energies)

end