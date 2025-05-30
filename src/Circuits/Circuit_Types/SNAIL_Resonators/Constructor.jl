@kwdef struct SNAILResonators <: Model
    params :: Dict
    order :: Vector
    hilbertspace :: HS.Hilbertspace
    
    Ĥ :: qt.QuantumObject
    n̂ₛ :: qt.QuantumObject
    n̂ᵣs :: Vector

    CandD_Ops :: Dict

    Stuff :: Dict
    
    dressed_states :: Dict
    dressed_energies :: Dict
end


function init(Eᶜ, Eʲ, Eˡ, α, Φᵉ, Eᵒˢᶜs, gs, Nₛ, Nᵣs; Nₛ_full = 120, Cavity_Names = [], κᵗᶜ = 1/(50*1000), κᵗᵈ = 1/(192.5*1000), κᶜᶜ = 1/(1000*1000), Model_Name = "Default", Save_Path = "", N_Junc = 3)
    if typeof(Nᵣs) <: Number
        Nᵣs = [Nᵣs]
    end

    if !(length(Eᵒˢᶜs) == length(Nᵣs) & length(Nᵣs) == length(gs))
        @error "The lengths of g, ω and Nᵣs not the same. Please Try again"
        return Nothing
    end 
    
    
    Components = Dict{Any, Any}()
    Interactions = []
    
    
    if length(Cavity_Names) != length(Nᵣs)
        Cavity_Names = []
        for i in 1:length(Nᵣs)
            name = "abcdefghijklmnop"[i]
            push!(Cavity_Names, "Mode "*name)
        end
    end
    
    order = [["SNAIL"];Cavity_Names]
    Components["SNAIL"] = HS.Elements.SNAILs.init(Eᶜ, Eʲ, Eˡ, α, Φᵉ, Nₛ_full, Nₛ; name = "SNAIL", N_Junc = N_Junc)

    for i in 1:length(Nᵣs)
        name = Cavity_Names[i]
        Components[name] = HS.Elements.Resonators.init(Eᵒˢᶜs[i], Nᵣs[i], name)
    end

    for i in 1:length(gs)
        name = Cavity_Names[i]
        interaction = Dict("ops"=>Dict(Components["SNAIL"].name => Components["SNAIL"].n̂, Components[name].name => Components[name].â+Components[name].â'), "g"=>gs[i])
        push!(Interactions, interaction)
    end
    
    hilbertspace = HS.init(Components, Interactions, order = order)
    
    params = Dict{Any, Any}("E_C"=>Eᶜ, "E_J"=>Eʲ, "E_L"=>Eˡ, "alpha" => α, "Flux" => Φᵉ, "E_oscs"=>Eᵒˢᶜs, "gs"=>gs, "Ns"=>Nₛ, "Nrs"=>Nᵣs, "Ns_full"=>Nₛ_full, "kappa_tc" => κᵗᶜ, "kappa_td" => κᵗᵈ, "kappa_cc" => κᶜᶜ, "Cavity_Names" => Cavity_Names, "Model_Name" => Model_Name, "Save_Path"=>Save_Path, "ModelType" => "SNAILResonators")

    n̂ₛ = HS.IdentityWrapper(hilbertspace, Dict("SNAIL"=>Components["SNAIL"].n̂), order = order)

    n̂ᵣs = []
    for i in 1:length(Nᵣs)
        name = Cavity_Names[i]
        push!(n̂ᵣs, HS.IdentityWrapper(hilbertspace, Dict(name=>Components[name].N̂), order = order))
    end
    
    CandD_Ops = Dict{Any, Any}()

    # Bare Collapse Operators
    #------------------------------------------------------------------------------------------------------------------------------------------------
    SNAIL_Dephasing = 0*Components["SNAIL"].Ĥ
    for i in 0:(Nₛ-1)
        SNAIL_Dephasing+=sqrt(2*κᵗᵈ)*(i)*qt.projection(Components["SNAIL"].dim, i, i)
    end
    CandD_Ops["Bare SNAIL Dephasing"] = Utils.IdentityWrapper(hilbertspace, Dict("SNAIL"=>SNAIL_Dephasing), order = order)

    SNAIL_Collapse = 0*Components["SNAIL"].Ĥ
    for i in 0:Nₛ-2
        ip1 = i+1
        SNAIL_Collapse+=sqrt(κᵗᶜ)*(ip1)*qt.projection(Components["SNAIL"].dim, i, ip1)
    end
    CandD_Ops["Bare SNAIL Collapse"] = Utils.IdentityWrapper(hilbertspace, Dict("SNAIL"=>SNAIL_Collapse), order = order)

    for mode in 1:length(Nᵣs)
        name = Cavity_Names[mode]
        Cavity_Collapse = 0*Components[name].Ĥ
        for i in 1:Nᵣs[mode]-1
            ip1 = i+1
            ψi = qt.fock(Components[name].dim, i-1)
            ψip1 = qt.fock(Components[name].dim, i-1+1)
            Cavity_Collapse += sqrt(κᶜᶜ)*(ip1)*ψi*ψip1'
        end
        CandD_Ops["Bare "*name*" Collapse"] = Utils.IdentityWrapper(hilbertspace, Dict(name=>Cavity_Collapse), order = order)
    end

    # Dressed Collapse Operators
    #------------------------------------------------------------------------------------------------------------------------------------------------
    # SNAIL_Dephasing = 0
    # list_for_iter = []
    # for N in Nᵣs
    #     push!(list_for_iter, 0:(N-1))
    # end
    # iter =  Iterators.product(list_for_iter...)
    # for i in 0:(Nₛ-1)
    #     for j in iter
    #         ψ =  hilbertspace.dressed_states[(i, j...)]
    #         SNAIL_Dephasing+=sqrt(2*κᵗᵈ)*(i)*ψ*ψ'
    #     end
    # end
    # CandD_Ops["Dressed SNAIL Dephasing"] = SNAIL_Dephasing

    # SNAIL_Collapse = 0
    # for i in 0:Nₛ-2
    #     ip1 = i+1
    #     for j in iter
    #         ψ =  hilbertspace.dressed_states[(i, j...)]
    #         ψp1 = hilbertspace.dressed_states[(i+1, j...)]
    #         SNAIL_Collapse+=sqrt(κᵗᶜ)*(i+1)*ψ*ψp1'
    #     end
    # end
    # CandD_Ops["Dressed SNAIL Collapse"] = SNAIL_Collapse


    # for N in Nᵣs
    #     push!(list_for_iter, 0:(N-1))
    # end
    # iter =  Iterators.product(list_for_iter...)    
    # for mode in 1:length(Nᵣs)
    #     name = Cavity_Names[mode]
    #     dims_to_itert = [Nₛ, deleteat!(copy(Nᵣs), mode)...]
    #     list_for_iter = []
    #     for i in length(dims_to_itert)
    #         push!(list_for_iter, 0:(dims_to_itert[i]-1))
    #     end
    #     iter = Iterators.product(list_for_iter...)

    #     Cavity_Collapse = 0
    #     for i in 0:(Nᵣs[mode]-2)
    #         for j in iter
    #             k = [j...]
    #             kp1 = [j...]
    #             ψ = hilbertspace.dressed_states[tuple(insert!(k, mode+1, i)...)]
    #             ψp1 = hilbertspace.dressed_states[tuple(insert!(kp1, mode+1, i+1)...)]
    #             Cavity_Collapse += sqrt(κᶜᶜ)*(i+1)*ψ*ψp1'
    #         end
    #     end
    #     CandD_Ops["Dressed "*name*" Collapse"] = Cavity_Collapse
    # end
    
    return SNAILResonators(params = params, hilbertspace=hilbertspace, n̂ₛ=n̂ₛ, Stuff = Dict{Any, Any}(), dressed_states = hilbertspace.dressed_states, dressed_energies = hilbertspace.dressed_energies, order = order, n̂ᵣs = n̂ᵣs, CandD_Ops = CandD_Ops, Ĥ = hilbertspace.Ĥ)

end