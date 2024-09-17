
function Square_Envelope(t)
    return 1.0
end

function Guassian_Envelope(t; sigma=1, mu = 0)
    σ = sigma
    μ = mu
    return exp(-(t-μ)^2/(2*σ^2))
end

function Linear_Envelope(t; t0=0, y0=0, m=0)
    return m*(t-t0)+y0
end

function Sine_Squared_Envelope(t; ramp_time = 0, offset = 0, phi = 0)
    return sin((π/2)*(t-offset)/ramp_time + phi)^2
end

function Guassian_Ramp_Envelope(t; pulse_time = 0, ramp_time = 10, sigma_factor = 2)
    σ = ramp_time/sigma_factor
    flat_top_time = pulse_time - 2*ramp_time

    if t<= ramp_time
        return Guassian_Envelope(t; sigma = 1, mu = ramp_time) 
    elseif (t>ramp_time) & (t<flat_top_time+ramp_time)
        return 1.0
    elseif (t>=flat_top_time+ramp_time)
        return Guassian_Envelope(t; sigma = 0, mu = flat_top_time+ramp_time)
    end
end

function Sine_Squared_Ramp_Envelope(t; pulse_time = 0, ramp_time = 10)
    flat_top_time = pulse_time-2*ramp_time
    if t<=ramp_time
        return Sine_Squared_Envelope(t; ramp_time = ramp_time)
    elseif (t>ramp_time) & (t<flat_top_time+ramp_time)
        return 1.0
    elseif t>=flat_top_time+ramp_time
        return Sine_Squared_Envelope(t; ramp_time = ramp_time, offset = flat_top_time+ramp_time, phi = π/2)
    end
end



Envelope_Dict = Dict{Any, Any}()
Envelope_Dict["Guassian_Ramp"] = Guassian_Ramp_Envelope
Envelope_Dict["Sine_Squared_Ramp"] = Sine_Squared_Ramp_Envelope
Envelope_Dict["Square"] = Square_Envelope
Envelope_Dict["Guassian"] = Guassian_Envelope
Envelope_Dict["Sine_Squared"] = Sine_Squared_Envelope
