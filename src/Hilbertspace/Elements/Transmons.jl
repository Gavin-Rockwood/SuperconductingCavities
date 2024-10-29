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

        N_full :: Int # U(1) Charge Number cutoff
        N :: Int # Number of Truncated Levels
        dim :: Int
        
        Ĥ_full :: qt.QuantumObject
        Ĥ :: qt.QuantumObject

        n̂_full :: qt.QuantumObject # Cut U(1) charge operator
        n̂ :: qt.QuantumObject # Truncated n operator

        eigsys_full ::  qt.EigsolveResult
        eigsys :: qt.EigsolveResult

        loss_ops :: Dict
    end

    function init(Eᶜ, Eʲ, N_full, N, name;  ng = 0, κᶜ = 1/(56*1000), κᵈ = 1.2348024316109425e-5)
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

        D̂ = 0*Ĥ
        for i in 1:(N-1) # this skips the 0 state becasue the coefficient is 0
            ψi = qt.fock(N, i)
            D̂ += sqrt(2*κᵈ)*sqrt(i)*ψi*ψi'
        end

        Ĉ = 0*Ĥ
        for i in 0:(N-2)
            ip1 = i+1
            ψi = qt.fock(N, i)
            ψip1 = qt.fock(N, ip1)
            Ĉ += sqrt(κᶜ)*sqrt(ip1)*ψi*ψip1'
        end
        
        loss_ops = Dict("Collapse" => Ĉ, "Dephasing" => D̂)

        return Transmon(name = name, Eᶜ = Eᶜ, Eʲ = Eʲ, ng = ng, N_full = N_full, N = N, dim = N, Ĥ_full = Ĥ_full, Ĥ = Ĥ, n̂_full = n̂_full, n̂ = n̂, eigsys_full = eigsys_full, eigsys = eigsys, loss_ops = loss_ops)
    end

end