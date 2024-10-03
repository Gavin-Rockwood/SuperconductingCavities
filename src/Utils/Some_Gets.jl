import QuantumToolbox as qt

function Get_Projection_Ops(dict_of_wavefunctions)
    res = deepcopy(dict_of_wavefunctions)
    map!(x -> x*x', values(res))
    return res
end


function Get_EVs(list_of_states, dict_of_operators)
    res = Dict{Any, Any}()

    for key in keys(dict_of_operators)
        op = dict_of_operators[key]

        res[key] = [qt.expect(op, ψ) for ψ in list_of_states]
    end
    return res
end