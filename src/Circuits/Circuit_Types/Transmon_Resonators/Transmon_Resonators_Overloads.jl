function Dynamics.FindStarkShift(model::TransmonResonators, 
    state1, 
    state2, 
    args...; kwargs...
    )

    ν = model.dressed_energies[state2]-model.dressed_energies[state1]
    ψ1 = model.dressed_states[state1]
    ψ2 = model.dressed_states[state2]

    state_names = [string(state1), string(state2)]
    
    Dynamics.FindStarkShift(model.hilbertspace, model.n̂ₜ, ψ1, ψ2, ν, args...;state_names = state_names, kwargs...)
end

function Dynamics.OptimizePulse(model::TransmonResonators, args...; kwargs...)
    Dynamics.OptimizePulse(model.Ĥ, model.n̂ₜ, args...; kwargs... )
end


function Dynamics.RunSingleOperator(model::TransmonResonators, args...; kwargs...)
    Dynamics.RunSingleOperator(model.Ĥ, model.n̂ₜ, args...; kwargs...)
end

function Dynamics.RunPulseSequence(model::TransmonResonators,
    state::Union{qt.QuantumObject{<:AbstractVector{T1}, qt.KetQuantumObject, 2}, qt.QuantumObject{<:AbstractArray{T1},qt.OperatorQuantumObject}},
    op_sequence,
    args...;
    kwargs...
    ) where T1<:Number

    Dynamics.RunPulseSequence(model.Ĥ, model.n̂ₜ, state, op_sequence, model.Stuff["op_drive_params"], args...; kwargs...)
end


function Utils.save_model(model::TransmonResonators; kwargs...)
    save_model(model; kwargs...)
end

function Dynamics.Get_Floquet_t0_Eigsys(model::TransmonResonators, args...; kwargs...)
    Dynamics.Get_Floquet_t0_Eigsys(model.hilbertspace, args...; kwargs...)
end

function Dynamics.Floquet_t0_Sweep(model::TransmonResonators, args...; kwargs...)
    Dynamics.Floquet_t0_Sweep(model.hilbertspace, model.n̂ₜ, args...; kwargs...)
end

function Dynamics.Get_Floquet_t0_Table(model::TransmonResonators, args...; kwargs...)
    Dynamics.Get_Floquet_t0_Table(model.hilbertspace, args...; kwargs...)
end

function Dynamics.Propagator(model::TransmonResonators, args...; kwargs...)
    Dynamics.Propagator(model.hilbertspace, args...; kwargs...)
end