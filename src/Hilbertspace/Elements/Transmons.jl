module Transmons

    import QuantumToolbox as qt
    using LinearAlgebra
    #using ProtoStructs

    export Transmon, init

    @kwdef struct Transmon
        name :: String
        Eᶜ :: Float64
        Eʲ :: Float64
        ng :: Real

        full_N :: Int # U(1) Charge Number cutoff
        N :: Int # Number of Truncated Levels
        dim :: Int
        
        full_Ĥ :: qt.QuantumObject
        Ĥ :: qt.QuantumObject

        full_n̂ :: qt.QuantumObject # Cut U(1) charge operator
        n̂ :: qt.QuantumObject # Truncated n operator

        full_eigsys ::  qt.EigsolveResult
        eigsys :: qt.EigsolveResult
    end

    function init(Eᶜ, Eʲ, N_full, N, name;  ng = 0)
        dim_full = 2*N_full+1
        𝕀̂_full = qt.eye(dim_full)
        
        jump_full = qt.tunneling(dim_full, 1)

        n̂_full = qt.num(dim_full) - N_full

        Ĥ_full = 4*Eᶜ*(ng*𝕀̂_full - n̂_full)^2 - 0.5*Eʲ*(jump_full)

        eigsys_full = qt.eigenstates(Ĥ_full)

        Π = zeros(ComplexF64, dim_full, N)
        for i in 1:N
            Π[:, i] = eigsys_full.vectors[:, i]
        end

        H⃗_full = Ĥ_full.data
        H⃗ = Π'*H⃗_full*Π
        n⃗_full = n̂_full.data
        n⃗ = Π'*n⃗_full*Π

        Ĥ = qt.Qobj(H⃗)
        
        herm_check = norm(Ĥ - Ĥ')
        if herm_check > 1e-9
            println("Herm_check for Ĥ Failed with value $herm_check")
        end

        Ĥ = 0.5*(Ĥ+Ĥ')
        
        n̂ = qt.Qobj(n⃗)
        
        herm_check = norm(n̂ - n̂')
        if herm_check > 1e-9
            println("Herm_check for n̂ Failed with value $herm_check")
        end

        n̂ = 0.5*(n̂+n̂')


        eigsys = qt.eigenstates(Ĥ)
        
        return Transmon(name = name, Eᶜ = Eᶜ, Eʲ = Eʲ, ng = ng, N_full = N_full, N = N, dim = N, Ĥ_full = Ĥ_full, Ĥ = Ĥ, n̂_full = n̂_full, n̂ = n̂, eigsys_full = eigsys_full, eigsys = eigsys)
    end

end