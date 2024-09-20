import QuantumToolbox as qt
using YAXArrays
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

function Base.collect(X::qt.QuantumObject)
    return X.data
end

function Qobj_List_To_DS(QOL; cube_name = :Default, cube_properties = Dict{Any, Any}(), ds_properties = Dict{Any, Any}(), step_name = :Step)
    step1 = collect.(QOL)
    
    
    if length(size(step1[1])) == 1
        step2 = collect.(reim.(step1));
        step3 = permutedims(cat(reduce.(hcat, step2)..., dims = 3), (3,1,2));
        dim = size(step3)
        axlist = (Dim{step_name}(1:dim[1]),  Dim{:i}(1:dim[2]), Dim{:P}(["Re", "Im"]))
    elseif length(size(step1[1])) > 1
        step2 = collect.(reim(step1))
        step3 = permutedims(cat(reduce.((x,y) -> cat(x, y, dims = 3), step2)..., dims = 4), (3, 1,2,4));
        
        dim = size(step3)
        
        axlist = (Dim{step_name}(1:dim[1]),  Dim{:i}(1:dim[2]), Dim{:j}(1:dim[3]), Dim{:P}(["Re", "Im"]))
    end
    
    ds = Dataset(; cube_name => YAXArray(axlist, step3, cube_properties), properties = ds_properties)
end