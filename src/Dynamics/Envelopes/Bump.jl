function Bump_Envelope(t; pulse_time = 10, k = 1, a = 5)
    x = (t-pulse_time/2)/(pulse_time/2)
    if x<=-1
        return 0
    elseif x>=1
        return 0
    elseif x == 0
        return 1
    else
        1/(1+exp(k*(a-(a+1)*abs(x))/(x^2-abs(x))))
    end
    
end
Envelope_Dict["Bump"] = Bump_Envelope

function Bump_Envelope_Cal(x...)
    t = x[1]
    Envelope_Args = x[2]

    Envelope_Args["pulse_time"] = t
    if !("k" in keys(Envelope_Args))
        Envelope_Args["k"] = 1
    end

    if !("a" in keys(Envelope_Args))
        Envelope_Args["a"] = 5
    end

    return Envelope_Args
end
Envelope_Dict_Cal["Bump"] = Bump_Envelope_Cal