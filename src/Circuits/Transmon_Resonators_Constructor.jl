#path_to_utils= join(split(@__DIR__, "/")[1:end-1], "/")*"/Utils"
#include(path_to_utils*"/HilbertSpaceUtils.jl")

include("Struct_Defs/Transmon_Resonators.jl")

export Transmon_Resonators_Constructor

function Transmon_Resonators_Constructor(Eᶜ, Eʲ, Eᵒˢᶜs, gs, Nₜ, Nᵣs; Nₜ_cut=60, ng = 0, Cavity_Names = [], κᵗᶜ = 1/(50*1000), κᵗᵈ = 1/(192.5*1000), κᶜᶜ = 1/(1000*1000))
    if !(length(Eᵒˢᶜs) == length(Nᵣs) & length(Nᵣs) == length(gs))
        @error "The lengths of g, ω and Nᵣs not the same. Please Try again"
        return Nothing
    end 
    
    order = [["Transmon"];Cavity_Names]

    Components = Dict{Any, Any}()
    Interactions = []



    if length(Cavity_Names) != length(Nᵣs)
        Cavity_Names = []
        for i in 1:length(Nᵣs)
            name = "abcdefghijklmnop"[i]
            push!(Cavity_Names, "Mode "*name)
        end
    end

    Components["Transmon"] = Init_Transmon(Eᶜ, Eʲ, Nₜ_cut, Nₜ, "Transmon")

    for i in 1:length(Nᵣs)
        name = Cavity_Names[i]
        Components[name] = Init_Resonator(Eᵒˢᶜs[i], Nᵣs[i], name)
    end

    for i in 1:length(gs)
        name = Cavity_Names[i]
        interaction = Dict("ops"=>Dict(Components["Transmon"].name => Components["Transmon"].n̂, Components[name].name => Components[name].â+Components[name].â'), "g"=>gs[i])
        push!(Interactions, interaction)
    end
    
    hilbertspace = Hilbertspace_Constructor(Components, Interactions, order = order)
    
    params = Dict{Any, Any}("Eᶜ"=>Eᶜ, "Eʲ"=>Eʲ, "Eᵒˢᶜs"=>Eᵒˢᶜs, "gs"=>gs, "Nₜ"=>Nₜ, "Nᵣs"=>Nᵣs, "Nₜ_cut"=>Nₜ_cut, "ng"=>ng)

    n̂ₜ = IdentityWrapper(hilbertspace, Dict("Transmon"=>Components["Transmon"].n̂), order = order)

    n̂ᵣs = []
    for i in 1:length(Nᵣs)
        name = Cavity_Names[i]
        push!(n̂ᵣs, IdentityWrapper(hilbertspace, Dict(name=>Components[name].N̂), order = order))
    end
    
    Collapse_Ops = Dict{Any, Any}()

    Transmon_Dephasing = 0*Components["Transmon"].Ĥ
    for i in 1:Nₜ
        Transmon_Dephasing+=sqrt(κᵗᵈ)*(i-1)*qo.transition(Components["Transmon"].𝔹, i, i)
    end
    Collapse_Ops["Transmon Dephasing"] = IdentityWrapper(hilbertspace, Dict("Transmon"=>Transmon_Dephasing), order = order)

    Transmon_Collapse = 0*Components["Transmon"].Ĥ
    for i in 1:Nₜ-1
        ip1 = i+1
        Transmon_Collapse+=sqrt(κᵗᶜ)*(ip1 - 1)*qo.transition(Components["Transmon"].𝔹, i, ip1)
    end
    Collapse_Ops["Transmon Collapse"] = IdentityWrapper(hilbertspace, Dict("Transmon"=>Transmon_Collapse), order = order)

    for mode in 1:length(Nᵣs)
        name = Cavity_Names[mode]
        Cavity_Collapse = 0*Components[name].Ĥ
        for i in 1:Nᵣs[mode]-1
            ip1 = i+1
            ψi = qo.fockstate(Components[name].𝔹, i-1)
            ψip1 = qo.fockstate(Components[name].𝔹, i-1+1)
            Cavity_Collapse += sqrt(κᶜᶜ)*(ip1-1)*ψi*ψip1'
        end
        Collapse_Ops[name*" Collapse"] = IdentityWrapper(hilbertspace, Dict(name=>Cavity_Collapse), order = order)
    end
    


    return Transmon_Resonators(params = params, hilbertspace=hilbertspace, n̂ₜ=n̂ₜ, Stuff = Dict{Any, Any}(), dressed_states = hilbertspace.dressed_states, dressed_energies = hilbertspace.dressed_energies, order = order, n̂ᵣs = n̂ᵣs, Collapse_Ops = Collapse_Ops, Ĥ = hilbertspace.Ĥ)
 
end

function Transmon_Resonators_Loader(file)
    saved_dict = JSON.parsefile(file)

    Eᶜ = saved_dict["Main_Config"]["E_C"]
    Eʲ = saved_dict["Main_Config"]["E_J"]
    Eᵒˢᶜs = saved_dict["Main_Config"]["E_osc"]
    gs = saved_dict["Main_Config"]["g"]
    ng = saved_dict["Main_Config"]["transmon_ng"]
    
    Nₜ = saved_dict["Main_Config"]["transmon_truncated_dim"]
    Nᵣs = saved_dict["Main_Config"]["resonator_truncated_dim"]
    Nₜ_cut = saved_dict["Main_Config"]["transmon_ncut"]

    model_name = saved_dict["Main_Config"]["model_name"]
    save_path = saved_dict["Main_Config"]["save_path"]

    Cavity_Names = saved_dict["Cavity_Names"]

    κᵗᶜ = saved_dict["Extra_Stuff"]["CandD_Coefs"]["T_C"]
    κᵗᵈ = saved_dict["Extra_Stuff"]["CandD_Coefs"]["T_D"]
    κᶜᶜ = saved_dict["Extra_Stuff"]["CandD_Coefs"]["C_C"]

    model = Transmon_Resonators_Constructor(Eᶜ, Eʲ, Eᵒˢᶜs, gs, Nₜ, Nᵣs, Nₜ_cut=Nₜ_cut, ng=ng, Cavity_Names=Cavity_Names, κᵗᶜ=κᵗᶜ, κᵗᵈ=κᵗᵈ, κᶜᶜ = κᶜᶜ)

    model.Stuff["op_drive_params"] = saved_dict["op_drive_params_dict"]
    model.Stuff["Drive Sequences"] = saved_dict["Drive Sequences"]
    model.Stuff["model_name"] = model_name
    model.Stuff["save_path"] = save_path

    return model
end