function Guassian_Envelope(t; sigma=1, mu = 0)
    σ = sigma
    μ = mu
    return exp(-(t-μ)^2/(2*σ^2))
end
Envelope_Dict["Guassian"] = Guassian_Envelope

function Guassian_Envelope_Cal(x...)
    t = x[1]
    Envelope_Args = x[2]

    if "sigma_factor" in keys(x[3])
        sigma_factor = x[3]["sigma_factor"]
    else
        sigma_factor = 4
    end

    Envelope_Args["sigma"] = t/sigma_factor
    Envelope_Args["mu"] = t/2

    return Envelope_Args
end
Envelope_Dict_Cal["Guassian"] = Guassian_Envelope_Cal