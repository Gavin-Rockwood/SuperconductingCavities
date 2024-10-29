function Sine_Squared_Ramp_Envelope(t; pulse_time = 0, ramp_time = 10)
    flat_top_time = pulse_time-2*ramp_time
    if t<=ramp_time
        return Sine_Squared_Envelope(t; ramp_time = ramp_time)
    elseif (t>ramp_time) & (t<flat_top_time+ramp_time)
        return 1.0+0*t # the 0*t is for autodif :D
    elseif t>=flat_top_time+ramp_time
        return Sine_Squared_Envelope(t; ramp_time = ramp_time, offset = flat_top_time+ramp_time, phi = π/2)
    end
end
Envelope_Dict["Sine_Squared_Ramp"] = Sine_Squared_Ramp_Envelope


function Sine_Squared_Ramp_Envelope_Cal(x...)
    t = x[1]
    Envelope_Args = x[2]
    Envelope_Args["pulse_time"] = t
    return Envelope_Args
end
Envelope_Dict_Cal["Sine_Squared_Ramp"] = Sine_Squared_Ramp_Envelope_Cal