export Transmon_Resonators


@kwdef struct Transmon_Resonators <: Model
    params :: Dict
    order :: Vector
    hilbertspace :: HilbertSpace
    
    Ĥ :: qo.Operator
    n̂ₜ :: qo.Operator
    n̂ᵣs :: Vector

    Collapse_Ops :: Dict

    Stuff :: Dict
    
    dressed_states :: Dict
    dressed_energies :: Dict
end