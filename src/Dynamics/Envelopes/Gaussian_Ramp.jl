function Gaussian_Ramp_Envelope(t; pulse_time = 0, ramp_time = 10, sigma_factor = 4)
    σ = ramp_time/sigma_factor
    flat_top_time = pulse_time - 2*ramp_time

    if t<= ramp_time
        return Guassian_Envelope(t; sigma = σ, mu = ramp_time)
    elseif (t>ramp_time) & (t<flat_top_time+ramp_time)
        return 1.0+0*t # the 0*t is for autodif :D
    elseif (t>=flat_top_time+ramp_time)
        return Guassian_Envelope(t; sigma = σ, mu = flat_top_time+ramp_time)
    end
end
Envelope_Dict["Gaussian_Ramp"] = Gaussian_Ramp_Envelope

function Gaussian_Ramp_Envelope_Cal(x...)
    t = x[1]
    Envelope_Args = x[2]

    Envelope_Args["pulse_time"] = t

    return Envelope_Args
end
Envelope_Dict_Cal["Gaussian_Ramp"] = Gaussian_Ramp_Envelope_Cal