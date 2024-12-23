module SNAILs

    import QuantumToolbox as qt
    using LinearAlgebra
    using Symbolics
    import NonlinearSolve as NLS
    #using ProtoStructs
    using ..Utils
    export Transmon, Init_Transmon

    @kwdef struct SNAIL
        name :: String
        Eᶜ :: Float64
        Eʲ :: Float64
        Eˡ :: Float64
        α :: Float64 #
        Φᵉ :: Float64 # Flux

        dim_full :: Int # U(1) Charge Number cutoff
        N :: Int # Number of Truncated Levels
        dim :: Int
        
        Ĥ_full :: qt.QuantumObject
        Ĥ :: qt.QuantumObject

        n̂_full :: qt.QuantumObject # dim for truncated U(1) charge operator
        n̂ :: qt.QuantumObject # Truncated n operator

        φ̂_full :: qt.QuantumObject
        φ̂ :: qt.QuantumObject

        eigsys_full ::  qt.EigsolveResult
        eigsys :: qt.EigsolveResult
    end


    function c(n)
        @variables φ α N Φᵉ

        U = -α*cos(φ) - N * cos((φ - Φᵉ)/N)

        D = Differential(φ)^n

        return expand_derivatives(D(U))
    end

    function get_c_coeffs_bare(N_val, α_val, Φᵉ_val)
        @variables φ N α Φᵉ

        c_syms = []
        for n in 1:6
            c_sym = c(n)
            push!(c_syms, c_sym)
        end

        f_to_min(φ_val, p) = [Symbolics.value(substitute(c_syms[1], Dict(φ=>φ_val[1], N=>N_val, α => α_val, Φᵉ => Φᵉ_val)))]

        prob = NLS.NonlinearProblem(f_to_min, [0.0], [])
        φ_min = NLS.solve(prob)[1]
        @debug "φ_min: $φ_min"


        cs = [0.0]
        
        for n in 2:6
            push!(cs, Symbolics.value(substitute(c_syms[n], Dict(φ=>φ_min, N=>N_val, α=>α_val, Φᵉ=>Φᵉ_val))))
        end

        return cs
    end

    function get_c_coeffs_dressed(N, α, Φᵉ, Eʲ, Eˡ)
        cs = get_c_coeffs_bare(N, α, 2*π*Φᵉ)
        @debug "bare cs: $cs"
        p = Eˡ/(Eˡ+cs[2]*Eʲ)

        c2_dr = p*cs[2]
        c3_dr = cs[3]
        c4_dr = cs[4]-3*cs[3]^2/cs[2]*(1-p)/p
        c5_dr = cs[5]-10*cs[3]*cs[4]/cs[2]*(1-p)/p+15*cs[3]^3/cs[2]^2*(1-p)^2/p^2
        c6_dr = cs[6] - (10*cs[4]^2+15*cs[5]*cs[3])/(cs[2]*p)*(1-p) + (105*cs[4]*cs[3]^2)/(cs[2]*p)^2*(1-p)^2-(105*cs[3]^4)/(cs[2]*p)^3*(1-p)^3

        return [0, c2_dr, c3_dr, c4_dr, c5_dr, c6_dr]
    end

    function init(Eᶜ, Eʲ, Eˡ, α, Φᵉ,  dim_full, N; name = "SNAIL", N_Junc = 3)
        cs = get_c_coeffs_dressed(N_Junc, α, Φᵉ, Eʲ, Eˡ)
        @debug "dressed cs: $cs"
        
        ν = sqrt(8*cs[2]*Eᶜ*Eʲ)
        φ_zpf = (2*Eᶜ/Eʲ/cs[2])^(1/4)
        @debug "φ_zpf: $φ_zpf"
        n_zpf = 1/2/φ_zpf
        @debug "n_zpf: $n_zpf"

        â_full = qt.destroy(dim_full)
        n̂_full = 1im*(â_full'-â_full)
        φ̂_full = φ_zpf * (â_full'+â_full)

        herm_check = norm(φ̂_full-φ̂_full')
     
        if herm_check > 1e-9
            println("Herm_check for φ̂_full: $herm_check")
        end
        Ĥ_full = ν*â_full'*â_full

        den = 2
        for n in 3:6
            den *= n
            Ĥ_full+=Eʲ*cs[n]/den*(φ̂_full^n)
        end

        eigsys_full = qt.eigenstates(Ĥ_full)
        
        Π = zeros(ComplexF64, dim_full, N)

        for i in 1:N
            Π[:, i] = eigsys_full.vectors[:,i]
        end

        H⃗_full = Ĥ_full.data
        H⃗ = Π'*H⃗_full*Π
        n⃗_full = n̂_full.data
        n⃗ = Π'*n⃗_full*Π
        φ⃗_full = φ̂_full.data
        φ⃗ = Π'*φ⃗_full*Π

        Ĥ = qt.Qobj(H⃗)
        herm_check = norm(Ĥ - Ĥ')
        if herm_check > 1e-9
            println("Herm_check for Ĥ Failed with value $herm_check")
        end
        Ĥ =  0.5*(Ĥ+Ĥ')

        n̂ =  qt.Qobj(n⃗)
        herm_check = norm(n̂-n̂')
        if herm_check > 1e-9
            println("Herm_check for n̂ Failed with value $herm_check")
        end
        n̂ = 0.5*(n̂+n̂')



        φ̂ =  qt.Qobj(φ⃗)
        herm_check = norm(φ̂ - φ̂')
        if herm_check > 1e-9
            println("Herm_check for n̂ Failed with value $herm_check")
        end
        φ̂ = 0.5*(φ̂+φ̂')

        eigsys = qt.eigenstates(Ĥ)

        return SNAIL(name = name,  Eᶜ = Eᶜ, Eʲ = Eʲ, Eˡ = Eˡ, α = α, Φᵉ = Φᵉ, dim_full = dim_full, N = N, dim = N,  Ĥ_full = Ĥ_full, Ĥ = Ĥ, n̂_full = n̂_full, n̂ = n̂, φ̂_full = φ̂_full, φ̂ = φ̂, eigsys_full = eigsys_full, eigsys = eigsys)
        

    end

end