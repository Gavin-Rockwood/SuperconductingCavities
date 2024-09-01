import QuantumOptics as qo
using LinearAlgebra
#using ProtoStructs

export Resonator, Init_Resonator 

@kwdef struct Resonator
    name :: String
    E :: Float64
    N :: Int

    𝔹 :: qo.FockBasis{Int64}
    Ĥ :: qo.Operator
    â :: qo.Operator
    N̂ :: qo.Operator

    eigsys :: Tuple
end

function Init_Resonator(E, N, name)
    𝔹 = qo.FockBasis(N-1)

    â = qo.destroy(𝔹)
    N̂ = â'*â
    Ĥ = E*N̂

    eigsys = qo.eigenstates(qo.dense(Ĥ))

    return Resonator(E = E, N = N, name = name, 𝔹 = 𝔹, Ĥ = Ĥ, N̂ = N̂, eigsys = eigsys, â = qo.destroy(𝔹))


end
