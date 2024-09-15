import QuantumToolbox as qt
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
    dim :: Int
    
    Ĥ_cut :: qt.QuantumObject
    Ĥ :: qt.QuantumObject

    n̂_cut :: qt.QuantumObject # Cut U(1) charge operator
    n̂ :: qt.QuantumObject # Truncated n operator

    eigsys_cut ::  qt.EigsolveResult
    eigsys :: qt.EigsolveResult
end

function Init_Transmon(Eᶜ, Eʲ, N_cut, N, name;  ng = 0)
    cut_dim = 2*N_cut+1
    𝕀̂_cut = qt.eye(cut_dim)
    
    jump_cut = qt.tunneling(cut_dim, 1)

    n̂_cut = qt.num(cut_dim) - N_cut

    Ĥ_cut = 4*Eᶜ*(ng*𝕀̂_cut - n̂_cut)^2 - 0.5*Eʲ*(jump_cut + jump_cut')

    eigsys_cut = qt.eigenstates(Ĥ_cut)

    Π = zeros(ComplexF64, cut_dim, N)
    for i in 1:N
        Π[:, i] = eigsys_cut.vectors[:, i]
    end

    H⃗_cut = Ĥ_cut.data
    H⃗ = Π'*H⃗_cut*Π
    n⃗_cut = n̂_cut.data
    n⃗ = Π'*n⃗_cut*Π

    Ĥ = qt.eye(N)
    Ĥ = qt.Qobj(H⃗)
    
    herm_check = norm(Ĥ - Ĥ')
    if herm_check > 1e-9
        print("Herm_check for Ĥ Failed with value $herm_check")
    end

    Ĥ = 0.5*(Ĥ+Ĥ')
    
    n̂ = qt.Qobj(n⃗)
    
    herm_check = norm(n̂ - n̂')
    if herm_check > 1e-9
        print("Herm_check for n̂ Failed with value $herm_check")
    end

    n̂ = 0.5*(n̂+n̂')


    eigsys = qt.eigenstates(Ĥ)
    
    return Transmon(name = name, Eᶜ = Eᶜ, Eʲ = Eʲ, ng = ng, N_cut = N_cut, N = N, dim = N, Ĥ_cut = Ĥ_cut, Ĥ = Ĥ, n̂_cut = n̂_cut, n̂ = n̂, eigsys_cut = eigsys_cut, eigsys = eigsys)
end
