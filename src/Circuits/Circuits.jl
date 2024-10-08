module Circuits
    using ..SuperconductingCavities: Model
    using ..Utils
    import ..Hilbertspaces
    import ..Dynamics
    include("Circuit_Types/Transmon_Resonators/Transmon_Resonators.jl")
    include("Circuit_Types/SNAIL_Resonators/SNAIL_Resonators.jl")
    include("CircuitUtils.jl")

    Circuit_Constructors = Dict{Any, Any}()
    Circuit_Constructors["TransmonResonators"] = Transmon_Resonators.init
    Circuit_Constructors["SNAILResonators"] = SNAIL_Resonators.init
end


