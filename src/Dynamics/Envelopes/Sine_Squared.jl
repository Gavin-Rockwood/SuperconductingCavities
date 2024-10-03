function Sine_Squared_Envelope(t; ramp_time = 0, offset = 0, phi = 0)
    return sin((π/2)*(t-offset)/ramp_time + phi)^2
end
Envelope_Dict["Sine_Squared"] = Sine_Squared_Envelope

function Sine_Squared_Envelope_Cal(x...)
    t = x[1]
    Envelope_Args = x[2]

    Envelope_Args["offset"] = t/2
    Envelope_Args["ramp_time"] = t/2

    return Envelope_Args
end
Envelope_Dict_Cal["Sine_Squared"] = Sine_Squared_Envelope_Cal