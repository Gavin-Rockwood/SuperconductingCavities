import QuantumOptics as qo
using Logging

using LinearAlgebra
using SparseArrays
import CairoMakie as cm

using ProtoStructs

import QuantumOptics.⊗
import QuantumOptics.*

import CSV
using JSON
import Tables

using YAXArrays

using MiniLoggers
using Dates
using Revise

import SuperconductingCavities as SC

# MiniLogger(minlevel = MiniLoggers.Info) |> global_logger
InfoLogger = MiniLogger(minlevel = MiniLoggers.Info)
ProgressLogger = MiniLogger(minlevel = LogLevel(-1))
DebugLogger = MiniLogger(minlevel = MiniLoggers.Debug)

global_logger(ProgressLogger)

function tostr(obj)
    io = IOBuffer()
    show(io, "text/plain", obj)
    String(take!(io))
end


Mode3 = SC.Transmon_Resonators_Loader("ModelSaves/Mode_3/Mode_3.json");

solver_kwargs = Dict{Any, Any}("reltol" => 1e-8, "abstol" => 1e-8, "tol"=>1e-8)
ψ = Mode3.dressed_states[(0,0)]
ρ = ψ*ψ'
start_time = now()

with_loss = false
c_ops = []
if with_loss
    c_ops = collect(values(Mode3.CandD_Ops))
end

@info "With Loss: $with_loss"

SC.RunPulseSequence(Mode3, ρ, Mode3.Stuff["Drive Sequences"]["Binomial_Code"], c_ops = c_ops, run_name = "Run_Loss_"*string(with_loss)*"_"*string(now()), solver_kwargs = solver_kwargs)
end_time = now()

@info "Total Run Time: "*string(end_time - start_time)