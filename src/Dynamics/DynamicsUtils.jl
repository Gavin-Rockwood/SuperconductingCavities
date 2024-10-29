import QuantumToolbox as qt
import ProgressMeter as PM
using DifferentialEquations
using OrdinaryDiffEqVerner: Vern9
import FFTW
import Interpolations as Interp

export Get_Drive_Coef, Get_Evelope, Get_Ĥ_D, Get_Lₜ, Propagator

"""
Gets the coefficient for when ν is a constant
"""
function Get_Drive_Coef(ν::T1,
    ε;
    envelope = Envelopes.Square_Envelope,
    drive_time = 0,
    return_ℂ = false
    )where T1<:Number
    
    function drive_coef(t, params...; return_ℂ = return_ℂ, kwargs...)
        if return_ℂ
            return 2*π*ε*envelope(t)*exp(-2π*ν*t*1im)
        else 
            return 2*π*ε*envelope(t)*sin(2π*ν*t)
        end
    end

    return drive_coef
end

"""
Gets the coefficient for when ν is a function of ε
"""
function Get_Drive_Coef(ν::T1,
    ε;
    envelope = Envelopes.Square_Envelope,
    drive_time = 0,
    t0 = 0,
    ϕ = nothing, # Incase you want to use a precomputed phase
    return_ϕ = false,
    return_ℂ = false
    )where T1<:Function

    εt(t) = ε*envelope(t)

    if !(typeof(ϕ) <: Function)
        ω(u, p, t) = 2*π*ν(εt(t))

        u0 = 0
        tspan = [t0, t0+drive_time]
        prob = ODEProblem(ω, u0, tspan)

        ϕ = solve(prob, Vern9(), abstol = 1e-9, reltol = 1e-9)
        if return_ϕ
            return ϕ
        end
    end

    function drive_coef(t, params...; return_ℂ = return_ℂ, kwargs...)
        if return_ℂ
            return 2*π*εt(t)*exp(-1im*ϕ(t))
        else
            return 2*π*εt(t)*sin(ϕ(t))
        end
    end
    return drive_coef
end

"""
This applies a low pass butterworth filter to the drive coefficient. The default is a cut off of 4 GHz with 20 poles. The sampling rate for the 
FFT defaults to 100GHz. This function using the imaginary part of the drive coefficient as I want the sine part! This setting is controlled
by the re_or_im keyword argument.
"""
function Get_Low_Pass_Filtered_Drive_Coef(drive_coef, drive_time; freq_cutoff = 4, poles = 20, fs = 1e2, re_or_im = imag, return_T = false)
    time = collect(0:1/fs:(drive_time+1/fs))
    signal = drive_coef.(time, return_ℂ = true)

    fft = FFTW.fft(signal)
    
    νs = collect(1:1:length(fft))*fs/length(fft)

    T(ν) = 1/sqrt(1+(ν/freq_cutoff)^(2*poles))

    if return_T
        return T
    end

    Y = fft.*T.(νs)

    ifft = FFTW.ifft(Y)
    interp_res = Interp.linear_interpolation(time, re_or_im.(ifft))
    return (t, params...) -> interp_res(t)
end


"""
This gets the drive hamiltonian for a given operator and drive coefficient
"""
function Get_Ĥ_D(op::qt.QuantumObject, drive_coef::Union{Nothing, Function}; TDOS = true, params = nothing, init_time = 0.0)
    drive_coef = (drive_coef isa Nothing) ? (t)->1 : drive_coef
    if TDOS
        return qt.TimeDependentOperatorSum([drive_coef], [op], params = params, init_time = init_time)
    else
        return t -> drive_coef(t)*op
    end
end

function Get_Lₜ(op, drive_coef; params = nothing, init_time = 0.0)
    return qt.TimeDependentOperatorSum([drive_coef], [qt.liouvillian(op)], params = params, init_time = init_time)
end

"""
Calculates the propagator for a given Hamiltonian Ĥₜ over a time [ti,tf]
"""
function Propagator(hilbertspace, Ĥₜ, tf; progress_meter = false, ti = 0)
    U = 0*Utils.eye_like(hilbertspace.Ĥ)

    p = PM.Progress(length(hilbertspace.dressed_states), enabled = progress_meter)
    for state in keys(hilbertspace.dressed_states)
        ψi = hilbertspace.dressed_states[state]
        se_res = qt.sesolve(2*π*hilbertspace.Ĥ, ψi, [ti, tf], H_t = Ĥₜ, progress_bar = false, alg = Vern9())
        ψf = se_res.states[end]
        U += ψf*ψi'
        PM.next!(p)
    end
    return U
end



"""
This function returns a function that can be used to chirp a frequency ν0 by the parameters chirp_params. 
"""
function chirper(ν0, chirp_params)
    return ε -> ν0+sum(chirp_params[i]*ε^i for i in 1:length(chirp_params))
end




