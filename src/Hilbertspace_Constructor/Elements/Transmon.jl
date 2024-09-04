import QuantumOptics as qo
using LinearAlgebra
#using ProtoStructs

export Transmon, Init_Transmon

@kwdef struct Transmon
    name :: String
    Eᶜ :: Float64
    Eʲ :: Float64
    ng :: Real

    N_cut :: Int # U(1) Charge Number cutoff
    N :: Int # Number of Truncated Levels
    
    𝔹_cut :: qo.NLevelBasis{Int64} # Cut U(1) Basis
    𝔹 :: qo.NLevelBasis{Int64} # Truncated Basis

    Ĥ_cut :: qo.Operator
    Ĥ :: qo.Operator

    n̂_cut :: qo.Operator # Cut U(1) charge operator
    n̂ :: qo.Operator # Truncated n operator

    eigsys_cut ::  Tuple
    eigsys :: Tuple
end

function Init_Transmon(Eᶜ, Eʲ, N_cut, N, name;  ng = 0)
    cut_dim = 2*N_cut+1
    𝔹_cut = qo.NLevelBasis(cut_dim)
    𝕀̂_cut = qo.identityoperator(𝔹_cut)
    
    jump_cut = 0*𝕀̂_cut
    for i in 1:(cut_dim-1)
        jump_cut +=  qo.transition(𝔹_cut, i,  i+1)
    end

    n̂_cut = 0*𝕀̂_cut
    for i in 1:cut_dim
        n̂_cut += qo.transition(𝔹_cut,  i, i)*(i-N_cut-1)
    end

    Ĥ_cut = 4*Eᶜ*(ng*𝕀̂_cut - n̂_cut)^2 - 0.5*Eʲ*(jump_cut + jump_cut')

    eigsys_cut = qo.eigenstates(qo.dense(Ĥ_cut))

    Π = zeros(ComplexF64, N, cut_dim)
    for i in 1:N
        Π[i, :] = eigsys_cut[2][i].data
    end

    H⃗_cut = Ĥ_cut.data
    H⃗ = Π*H⃗_cut*Π'
    n⃗_cut = n̂_cut.data
    n⃗ = Π*n⃗_cut*Π'

    𝔹 = qo.NLevelBasis(N)
    
    Ĥ = qo.dense(qo.identityoperator(𝔹))
    Ĥ.data = H⃗
    
    herm_check = norm((Ĥ - Ĥ').data)
    if herm_check > 1e-9
        print("Herm_check for Ĥ Failed with value $herm_check")
    end

    Ĥ = 0.5*(Ĥ+Ĥ')
    
    n̂ = 0*Ĥ
    for i in 1:N
        n̂ += (i-1)*qo.transition(𝔹, i, i)
    end

    n̂.data = n⃗
    
    herm_check = norm((n̂ - n̂').data)
    if herm_check > 1e-9
        print("Herm_check for n̂ Failed with value $herm_check")
    end

    n̂ = 0.5*(n̂+n̂')


    eigsys = qo.eigenstates(qo.dense(Ĥ))
    
    return Transmon(name = name, Eᶜ = Eᶜ, Eʲ = Eʲ, ng = ng, N_cut = N_cut, N = N, 𝔹_cut = 𝔹_cut, 𝔹 = 𝔹, Ĥ_cut = Ĥ_cut, Ĥ = Ĥ, n̂_cut = n̂_cut, n̂ = n̂, eigsys_cut = eigsys_cut, eigsys = eigsys)
end
