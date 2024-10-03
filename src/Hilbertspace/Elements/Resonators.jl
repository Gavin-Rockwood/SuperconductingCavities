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
    end

    function init(E, N, name)

        â = qt.destroy(N)
        N̂ = â'*â
        Ĥ = E*N̂

        eigsys = qt.eigenstates(Ĥ)

        return Resonator(E = E, N = N, dim = N, name = name, Ĥ = Ĥ, N̂ = N̂, eigsys = eigsys, â = â)
    end

end