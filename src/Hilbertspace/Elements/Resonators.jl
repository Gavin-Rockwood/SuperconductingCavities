module Resonators
    import QuantumToolbox as qt
    using LinearAlgebra


    export Resonator, init
    @kwdef struct Resonator
        name :: String
        E :: Float64
        N :: Int
        dim :: Int

        Ĥ :: qt.QuantumObject
        â :: qt.QuantumObject
        N̂ :: qt.QuantumObject

        eigsys :: qt.EigsolveResult

        loss_ops :: Dict
    end

    function init(E, N, name; κᶜ = 1/(1000*1000), κᵈ =  0)

        â = qt.destroy(N)
        N̂ = â'*â
        Ĥ = E*N̂

        eigsys = qt.eigenstates(Ĥ)

        # Cavity Collapse
        Ĉ = 0*Ĥ
        for i in 0:(N-2)
            ip1 = i+1
            ψi = qt.fock(N, i)
            ψip1 = qt.fock(N, ip1)
            Ĉ += sqrt(κᶜ)*sqrt(ip1)*ψi*ψip1'
        end

        D̂ = 0*Ĥ
        for i in 1:(N-1) # this skips the 0 state becasue the coefficient is 0
            ψi = qt.fock(N, i)
            D̂ += sqrt(2*κᵈ)*sqrt(i)*ψi*ψi'
        end

        loss_ops = Dict("Collapse" => Ĉ, "Dephasing" => D̂)

        return Resonator(E = E, N = N, dim = N, name = name, Ĥ = Ĥ, N̂ = N̂, eigsys = eigsys, â = â, loss_ops = loss_ops)
    end

end