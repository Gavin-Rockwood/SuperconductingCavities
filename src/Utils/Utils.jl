import QuantumToolbox as qt

include("HilbertSpaceUtils.jl")
include("IO.jl")
include("RunOperatorSequence.jl")
include("Floquet.jl")
include("StateTracking.jl")
include("PulseFinder.jl")
include("TimeEvolutionUtils.jl")
include("PlotStateEvolution.jl")
#include("RunExperiment.jl")


function tostr(obj)
    io = IOBuffer()
    show(io, "text/plain", obj)
    String(take!(io))
end


function eye_like(op::qt.QuantumObject)
    return qt.eye(size(op)[1], dims = op.dims)
end