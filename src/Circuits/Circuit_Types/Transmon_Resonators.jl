import QuantumToolbox as qt
import JSON3
export Transmon_Resonators, Transmon_Resonators_Constructor, Transmon_Resonators_Loader


@kwdef struct Transmon_Resonators <: Model
    params :: Dict
    order :: Vector
    hilbertspace :: HilbertSpace
    
    Ĥ :: qt.QuantumObject
    n̂ₜ :: qt.QuantumObject
    n̂ᵣs :: Vector

    CandD_Ops :: Dict

    Stuff :: Dict
    
    dressed_states :: Dict
    dressed_energies :: Dict
end


function Transmon_Resonators_Constructor(Eᶜ, Eʲ, Eᵒˢᶜs, gs, Nₜ, Nᵣs; Nₜ_cut=60, ng = 0, Cavity_Names = [], κᵗᶜ = 1/(50*1000), κᵗᵈ = 1/(192.5*1000), κᶜᶜ = 1/(1000*1000), Model_Name = "Default", Save_Path = "")
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
    
    params = Dict{Any, Any}("E_C"=>Eᶜ, "E_J"=>Eʲ, "E_oscs"=>Eᵒˢᶜs, "gs"=>gs, "Nt"=>Nₜ, "Nrs"=>Nᵣs, "Nt_cut"=>Nₜ_cut, "ng"=>ng, "kappa_tc" => κᵗᶜ, "kappa_td" => κᵗᵈ, "kappa_cc" => κᶜᶜ, "Cavity_Names" => Cavity_Names, "Model_Name" => Model_Name, "Save_Path"=>Save_Path)

    n̂ₜ = IdentityWrapper(hilbertspace, Dict("Transmon"=>Components["Transmon"].n̂), order = order)

    n̂ᵣs = []
    for i in 1:length(Nᵣs)
        name = Cavity_Names[i]
        push!(n̂ᵣs, IdentityWrapper(hilbertspace, Dict(name=>Components[name].N̂), order = order))
    end
    
    CandD_Ops = Dict{Any, Any}()

    Transmon_Dephasing = 0*Components["Transmon"].Ĥ
    for i in 0:(Nₜ-1)
        Transmon_Dephasing+=sqrt(2*κᵗᵈ)*(i)*qt.projection(Components["Transmon"].dim, i, i)
    end
    CandD_Ops["Transmon Dephasing"] = IdentityWrapper(hilbertspace, Dict("Transmon"=>Transmon_Dephasing), order = order)

    Transmon_Collapse = 0*Components["Transmon"].Ĥ
    for i in 0:Nₜ-2
        ip1 = i+1
        Transmon_Collapse+=sqrt(κᵗᶜ)*(ip1)*qt.projection(Components["Transmon"].dim, i, ip1)
    end
    CandD_Ops["Transmon Collapse"] = IdentityWrapper(hilbertspace, Dict("Transmon"=>Transmon_Collapse), order = order)

    for mode in 1:length(Nᵣs)
        name = Cavity_Names[mode]
        Cavity_Collapse = 0*Components[name].Ĥ
        for i in 1:Nᵣs[mode]-1
            ip1 = i+1
            ψi = qt.fock(Components[name].dim, i-1)
            ψip1 = qt.fock(Components[name].dim, i-1+1)
            Cavity_Collapse += sqrt(κᶜᶜ)*(ip1)*ψi*ψip1'
        end
        CandD_Ops[name*" Collapse"] = IdentityWrapper(hilbertspace, Dict(name=>Cavity_Collapse), order = order)
    end
    
    return Transmon_Resonators(params = params, hilbertspace=hilbertspace, n̂ₜ=n̂ₜ, Stuff = Dict{Any, Any}(), dressed_states = hilbertspace.dressed_states, dressed_energies = hilbertspace.dressed_energies, order = order, n̂ᵣs = n̂ᵣs, CandD_Ops = CandD_Ops, Ĥ = hilbertspace.Ĥ)
 
end

function Transmon_Resonators_Loader_Python(file)
    #saved_dict = JSON.parsefile(file)
    saved_dict = JSON3.read(file)
    Eᶜ = saved_dict["Main_Config"]["E_C"]
    Eʲ = saved_dict["Main_Config"]["E_J"]
    Eᵒˢᶜs = saved_dict["Main_Config"]["E_osc"]
    gs = saved_dict["Main_Config"]["g"]
    ng = saved_dict["Main_Config"]["transmon_ng"]
    
    Nₜ = saved_dict["Main_Config"]["transmon_truncated_dim"]
    Nᵣs = saved_dict["Main_Config"]["resonator_truncated_dim"]
    Nₜ_cut = saved_dict["Main_Config"]["transmon_ncut"]

    Model_Name = saved_dict["Main_Config"]["model_name"]
    Save_Path = saved_dict["Main_Config"]["save_path"]

    Cavity_Names = saved_dict["Cavity_Names"]

    κᵗᶜ = saved_dict["Extra_Stuff"]["CandD_Coefs"]["T_C"]
    κᵗᵈ = saved_dict["Extra_Stuff"]["CandD_Coefs"]["T_D"]
    κᶜᶜ = saved_dict["Extra_Stuff"]["CandD_Coefs"]["C_C"]

    model = Transmon_Resonators_Constructor(Eᶜ, Eʲ, Eᵒˢᶜs, gs, Nₜ, Nᵣs, Nₜ_cut=Nₜ_cut, ng=ng, Cavity_Names=Cavity_Names, κᵗᶜ=κᵗᶜ, κᵗᵈ=κᵗᵈ, κᶜᶜ = κᶜᶜ, Model_Name = Model_Name, Save_Path =Save_Path)

    if "op_drive_params_dict" in keys(saved_dict)
        model.Stuff["op_drive_params"] = saved_dict["op_drive_params_dict"] 
    else
        model.Stuff["op_drive_params"] = Dict{Any, Any}()
    end
    if "Drive_Sequences" in keys(saved_dict)
        model.Stuff["Drive_Sequences"] = saved_dict["Drive_Sequences"]
    else
        model.Stuff["Drive_Sequences"] = Dict{Any, Any}()
    end

    return model
end


function Save_Model(model::Transmon_Resonators; model_name = nothing, save_path = nothing)
    save_dict = Dict{Any, Any}()
    save_dict["Main_Config"] = model.params
    save_dict["Stuff"] = Dict{Any, Any}()
    
    if "op_drive_params" in keys(model.Stuff)
        save_dict["Stuff"]["op_drive_params"] = model.Stuff["op_drive_params"]
    end
    if "Drive_Sequences" in keys(model.Stuff)
        save_dict["Stuff"]["Drive_Sequences"] = model.Stuff["Drive_Sequences"]
    end

    save_path = model.params["Save_Path"]
    filename = save_path*model.params["Model_Name"]

    
    if !isdir(save_path)
        @info "Making Save Directory $save_path"
        mkdir(save_path)
    end

    #### Making Backup
    try 
        cp(filename*".json", filename*"_BACKUP"*".json")
        rm(filename*".json")
    catch
    end

    open(filename*".json", "w") do f
        JSON3.pretty(f, save_dict)
    end
    try 
        rm(filename*"_BACKUP"*".json")
    catch
    end


end

function Transmon_Resonators_Loader(file)
    #saved_dict = JSON.parsefile(file)
    saved_dict = JSON3.read(file)
    
    Eᶜ = saved_dict["Main_Config"]["E_C"]
    Eʲ = saved_dict["Main_Config"]["E_J"]
    Eᵒˢᶜs = saved_dict["Main_Config"]["E_oscs"]
    gs = saved_dict["Main_Config"]["gs"]
    ng = saved_dict["Main_Config"]["ng"]
    
    Nₜ = saved_dict["Main_Config"]["Nt"]
    Nᵣs = saved_dict["Main_Config"]["Nrs"]
    Nₜ_cut = saved_dict["Main_Config"]["Nt_cut"]

    κᵗᶜ = saved_dict["Main_Config"]["kappa_tc"]
    κᵗᵈ = saved_dict["Main_Config"]["kappa_td"]
    κᶜᶜ = saved_dict["Main_Config"]["kappa_cc"]

    Model_Name = saved_dict["Main_Config"]["Model_Name"]
    Save_Path = saved_dict["Main_Config"]["Save_Path"]

    Cavity_Names = saved_dict["Main_Config"]["Cavity_Names"]

    

    model = Transmon_Resonators_Constructor(Eᶜ, Eʲ, Eᵒˢᶜs, gs, Nₜ, Nᵣs, Nₜ_cut=Nₜ_cut, ng=ng, Cavity_Names=Cavity_Names, κᵗᶜ=κᵗᶜ, κᵗᵈ=κᵗᵈ, κᶜᶜ = κᶜᶜ, Model_Name = Model_Name, Save_Path = Save_Path)
    
    if "op_drive_params" in keys(saved_dict["Stuff"])
        model.Stuff["op_drive_params"] = saved_dict["Stuff"]["op_drive_params"]
    else 
        model.Stuff["op_drive_params"] = Dict{Any, Any}()

    end
    if "Drive_Sequences" in keys(saved_dict["Stuff"])
        model.Stuff["Drive_Sequences"] = saved_dict["Stuff"]["Drive_Sequences"]
    else
        model.Stuff["Drive_Sequences"] = Dict{Any, Any}()
    end

    return model
end
