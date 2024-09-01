#path_to_CC = join(split(@__DIR__, "/")[1:end-2], "/")*"/CircuitConstruction"
#include(path_to_CC*"/Hilbertspace_Struct.jl")

export IdentityWrapper


function IdentityWrapper(𝕀̂_Dict::Dict, Operator_Dict; order = [])
    if length(order) == length(𝕀̂_Dict)
        key_list = order
    else
        key_list = collect(keys(𝕀̂_Dict))
    end

    op_Dict = deepcopy(𝕀̂_Dict)

    for key in keys(Operator_Dict)
        op_Dict[key] = Operator_Dict[key]
    end

    op_vec = []
    for i in 1:length(key_list)
        push!(op_vec, op_Dict[key_list[i]])
    end

    return qo.tensor(op_vec...)
end

function IdentityWrapper(hilbertspace::HilbertSpace, Operator_Dict; order = [])
    if length(order) == length(hilbertspace.𝕀̂_Dict)
        key_list = order
    else
        key_list = collect(keys(hilbertspace.𝕀̂_Dict))
    end

    op_Dict = deepcopy(hilbertspace.𝕀̂_Dict)

    for key in keys(Operator_Dict)
        op_Dict[key] = Operator_Dict[key]
    end
    op_vec = []
    for i in 1:length(key_list)
        push!(op_vec, op_Dict[key_list[i]])
    end

    return qo.tensor(op_vec...)
end