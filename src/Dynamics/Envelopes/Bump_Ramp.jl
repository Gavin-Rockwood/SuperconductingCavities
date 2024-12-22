function Bump_Ramp_Envelope(t; pulse_time = 10, ramp_time = 1, k = 2)
    if t<ramp_time
        return Bump_Envelope(t, pulse_time = 2*ramp_time, k = k)
    elseif (t<=(pulse_time-ramp_time)) && (t>=ramp_time)
        return 1
    elseif t>(pulse_time-ramp_time)
        return Bump_Envelope(t, pulse_time = 2*ramp_time, k = k, center = pulse_time-ramp_time)
    end
end
Envelope_Dict["Bump_Ramp"] = Bump_Ramp_Envelope

function Bump_Ramp_Envelope_Cal(x...)
    t = x[1]
    Envelope_Args = x[2]

    Envelope_Args["pulse_time"] = t
    if !("k" in keys(Envelope_Args))
        Envelope_Args["k"] = 2
    end
    
    return Envelope_Args
end

Envelope_Dict_Cal["Bump_Ramp"] = Bump_Ramp_Envelope_Cal