module Envelopes
    Envelope_Dict = Dict{Any, Any}()
    Envelope_Dict_Cal = Dict{Any, Any}()

    include("Square.jl")
    include("Guassian.jl")
    include("Guassian_Ramp.jl")
    include("Sine_Squared.jl")
    include("Sine_Squared_Ramp.jl")

    function Get_Envelope(envelope_name, envelope_kwargs)
        envelope_kwargs_sym = Dict{Symbol, Any}()
        for key in keys(envelope_kwargs)
            envelope_kwargs_sym[Symbol(key)] = envelope_kwargs[key]
        end
        function envelope(t)
            return Envelope_Dict[envelope_name](t; envelope_kwargs_sym...)
        end
        return envelope
    end

end