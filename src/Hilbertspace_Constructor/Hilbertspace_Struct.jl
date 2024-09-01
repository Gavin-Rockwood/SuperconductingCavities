import QuantumOptics as qo
using ProtoStructs

export HilbertSpace

@kwdef struct HilbertSpace
    Components :: Dict
    Interactions :: Vector
    𝕀̂_Dict :: Dict
    Ĥ :: qo.Operator

    dressed_states :: Dict
    dressed_energies :: Dict
end
