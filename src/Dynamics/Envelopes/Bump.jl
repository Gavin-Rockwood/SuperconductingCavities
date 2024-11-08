function Bump_Envelope(t; pulse_time = 10, k = 2, center = "None")
    if center == "None"
        center = pulse_time/2
    end
    x = (t-center)/(pulse_time/2)
    if x<=-1
        return 0
    elseif x>=1
        return 0
    elseif x == 0
        return 1
    else
        return exp(k*x^2/(x^2-1))
    end
    
end
Envelope_Dict["Bump"] = Bump_Envelope

function Bump_Envelope_Cal(x...)
    t = x[1]
    Envelope_Args = x[2]

    Envelope_Args["pulse_time"] = t
    if !("k" in keys(Envelope_Args))
        Envelope_Args["k"] = 2
    end

    return Envelope_Args
end
Envelope_Dict_Cal["Bump"] = Bump_Envelope_Cal