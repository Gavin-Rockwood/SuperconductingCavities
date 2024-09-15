import QuantumToolbox as qt
using ProtoStructs

export HilbertSpace

@kwdef struct HilbertSpace
    Components :: Dict
    Interactions :: Vector
    𝕀̂_Dict :: Dict
    𝕀̂ :: qt.QuantumObject
    Ĥ :: qt.QuantumObject

    dressed_states :: Dict
    dressed_energies :: Dict
end
