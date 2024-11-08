@kwdef struct TransmonResonators <: Model
    params :: Dict
    order :: Vector
    hilbertspace :: HS.Hilbertspace
    
    Ĥ :: qt.QuantumObject
    n̂ₜ :: qt.QuantumObject
    n̂ᵣs :: Vector

    CandD_Ops :: Dict

    Stuff :: Dict
    
    dressed_states :: Dict
    dressed_energies :: Dict
end


function init(Eᶜ, Eʲ, Eᵒˢᶜs, gs, Nₜ, Nᵣs; Nₜ_cut=60, ng = 0, Cavity_Names = [], κᵗᶜ = 1/(56*1000), κᵗᵈ = 1.2348024316109425e-5, κᶜᶜ = [1/(1000*1000)], κᶜᵈ = [0], Model_Name = "Default", Save_Path = "")
    if typeof(Nᵣs) <: Number
        Nᵣs = [Nᵣs]
    end
    if typeof(Eᵒˢᶜs) <: Number
        Eᵒˢᶜs = [Eᵒˢᶜs]
    end
    if typeof(gs) <: Number
        gs = [gs]
    end
    if typeof(κᶜᶜ) <: Number
        κᶜᶜ = [κᶜᶜ]
    end
    if typeof(κᶜᵈ) <: Number
        κᶜᵈ = [κᶜᵈ]
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
    order = [["Transmon"];Cavity_Names]
    
    Components["Transmon"] = HS.Elements.Transmons.init(Eᶜ, Eʲ, Nₜ_cut, Nₜ, "Transmon", κᵈ = κᵗᵈ, κᶜ = κᵗᶜ, ng = ng)

    for i in 1:length(Nᵣs)
        name = Cavity_Names[i]
        Components[name] = HS.Elements.Resonators.init(Eᵒˢᶜs[i], Nᵣs[i], name, κᵈ = κᶜᵈ[i], κᶜ = κᶜᶜ[i])
    end

    for i in 1:length(gs)
        name = Cavity_Names[i]
        interaction = Dict("ops"=>Dict(Components["Transmon"].name => Components["Transmon"].n̂, Components[name].name => Components[name].â+Components[name].â'), "g"=>gs[i])
        push!(Interactions, interaction)
    end
    
    hilbertspace = HS.init(Components, Interactions, order = order)
    
    params = Dict{Any, Any}("E_C"=>Eᶜ, "E_J"=>Eʲ, "E_oscs"=>Eᵒˢᶜs, "gs"=>gs, "Nt"=>Nₜ, "Nrs"=>Nᵣs, "Nt_cut"=>Nₜ_cut, "ng"=>ng, "kappa_tc" => κᵗᶜ, "kappa_td" => κᵗᵈ, "kappa_cc" => κᶜᶜ, "Cavity_Names" => Cavity_Names, "Model_Name" => Model_Name, "Save_Path"=>Save_Path, "ModelType" => "TransmonResonators")

    n̂ₜ = HS.IdentityWrapper(hilbertspace, Dict("Transmon"=>Components["Transmon"].n̂), order = order)

    n̂ᵣs = []
    for i in 1:length(Nᵣs)
        name = Cavity_Names[i]
        push!(n̂ᵣs, HS.IdentityWrapper(hilbertspace, Dict(name=>Components[name].N̂), order = order))
    end
    
    # CandD_Ops = Dict{Any, Any}()

    # # Bare Collapse Operators
    # #------------------------------------------------------------------------------------------------------------------------------------------------
    # # Transmon_Dephasing = 0*Components["Transmon"].Ĥ
    # # for i in 0:(Nₜ-1)
    # #     Transmon_Dephasing+=sqrt(2*κᵗᵈ)*sqrt(i)*qt.projection(Components["Transmon"].dim, i, i)
    # # end
    # CandD_Ops["Transmon Dephasing"] = Utils.IdentityWrapper(hilbertspace, Dict("Transmon"=>Components["Transmon"].loss_ops["Dephasing"]), order = order)

    # # Transmon_Collapse = 0*Components["Transmon"].Ĥ
    # # for i in 0:Nₜ-2
    # #     ip1 = i+1
    # #     Transmon_Collapse+=sqrt(κᵗᶜ)*sqrt(ip1)*qt.projection(Components["Transmon"].dim, i, ip1)
    # # end
    # CandD_Ops["Transmon Collapse"] = Utils.IdentityWrapper(hilbertspace, Dict("Transmon"=>Components["Transmon"].loss_ops["Collapse"]), order = order)

    # for mode in 1:length(Nᵣs)
    #     name = Cavity_Names[mode]
    #     CandD_Ops[name*" Collapse"] = Utils.IdentityWrapper(hilbertspace, Dict(name=>Components[name].loss_ops["Collapse"]), order = order)
    # end

    loss_ops = Dict{Any, Any}()
    for comp in keys(Components)
        loss_ops[comp*" Collapse"] = Utils.IdentityWrapper(hilbertspace, Dict(comp=>Components[comp].loss_ops["Collapse"]), order = order)
        loss_ops[comp*" Dephasing"] = Utils.IdentityWrapper(hilbertspace, Dict(comp=>Components[comp].loss_ops["Dephasing"]), order = order)
    end


    return TransmonResonators(params = params, hilbertspace=hilbertspace, n̂ₜ=n̂ₜ, Stuff = Dict{Any, Any}(), dressed_states = hilbertspace.dressed_states, dressed_energies = hilbertspace.dressed_energies, order = order, n̂ᵣs = n̂ᵣs, CandD_Ops = loss_ops, Ĥ = hilbertspace.Ĥ)

end