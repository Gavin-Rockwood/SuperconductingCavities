import QuantumToolbox as qt
using Logging

using MiniLoggers
using Dates


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

solver_kwargs = Dict{Any, Any}("reltol" => 1e-8, "abstol" => 1e-8)#, "tol"=>1e-8)
ψ = (Mode3.dressed_states[(0,0)]+Mode3.dressed_states[(1,0)])/sqrt(2)
ρ = ψ*ψ'
start_time = now()

with_loss = true
c_ops = []
if with_loss
    c_ops = collect(values(Mode3.CandD_Ops))
end

pulse_args = deepcopy(Mode3.Stuff["op_drive_params"]["q_g_0"])
pulse_args["pulse_time"] = 30*pulse_args["pulse_time"]
pulse_args["epsilon"] = 0.0

@info "Running Ramsey"

SC.RunSingleOperator(Mode3, ρ, pulse_args, c_ops = c_ops, run_name = "Ramsey_00_plus_10_"*string(now()), solver_kwargs = solver_kwargs, spps = 1, save_step = true)
end_time = now()

@info "Total Run Time: "*string(Dates.canonicalize(end_time - start_time))