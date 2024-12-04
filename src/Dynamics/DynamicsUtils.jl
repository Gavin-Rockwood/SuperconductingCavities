import QuantumToolbox as qt
import ProgressMeter as PM
using DifferentialEquations
using OrdinaryDiffEqVerner: Vern9
import FFTW
import Interpolations as Interp

export Get_Drive_Coef, Get_Evelope, Get_Ĥ_D, Get_Lₜ, Propagator

"""
    Get_Drive_Coef(ν::T1, ε; envelope=Envelopes.Square_Envelope, drive_time=0, return_ϕ=false, return_ℂ=false) where T1 <: Number

Generates a function to compute the drive coefficient for a given frequency and envelope.

# Arguments
- `ν::T1`: The frequency of the drive (in Hz).
- `ε`: The amplitude or strength of the drive.
- `envelope=Envelopes.Square_Envelope`: A callable defining the envelope shape of the drive. Default is a square envelope.
- `drive_time=0`: An optional parameter that could represent the duration of the drive, though it's currently unused in the function.
- `return_ϕ=false`: If true, the function returns the phase function ϕ(t) instead of the drive coefficient.
- `return_ℂ=false`: If true, the function returns the complex form of the drive coefficient.

# Returns
A function that computes the drive coefficient based on the provided parameters. If `return_ϕ` is true,
it returns a function for ϕ(t). Otherwise, it returns a function for the drive coefficient.

"""
function Get_Drive_Coef(ν::T1,
    ε;
    envelope = Envelopes.Square_Envelope,
    drive_time = 0,
    return_ϕ = false,
    return_ℂ = false
    )where T1<:Number
    
    ϕ(t) = 2π*ν*t
    if return_ϕ
        return ϕ
    end

    function drive_coef(t, params...; return_ℂ = return_ℂ, kwargs...)
        if return_ℂ
            return 2*π*ε*envelope(t)*exp(-1im*ϕ(t))
        else 
            return 2*π*ε*envelope(t)*sin(ϕ(t))
        end
    end

    return drive_coef
end


"""
    Get_Drive_Coef(ν::T1, ε; envelope = Envelopes.Square_Envelope, drive_time = 0, t0 = 0, ϕ = nothing, return_ϕ = false, return_ℂ = false) where T1<:Function

Creates a drive coefficient function based on the provided parameters.

# Arguments
- `ν::T1`: A function that returns the frequency based on the envelope-modulated coefficient.
- `ε`: A scalar coefficient.
- `envelope`: A function that defines the envelope shape of the drive signal. Defaults to `Envelopes.Square_Envelope`.
- `drive_time`: A scalar representing the drive time. Defaults to `0`.
- `t0`: The initial time. Defaults to `0`.
- `ϕ`: A precomputed phase function. If `nothing`, the phase will be computed internally. Defaults to `nothing`.
- `return_ϕ`: A boolean flag indicating whether to return the computed phase function. Defaults to `false`.
- `return_ℂ`: A boolean flag indicating whether to return the complex exponential form. Defaults to `false`.

# Returns
- If `return_ϕ` is `true`, returns the computed phase function `ϕ`.
- Otherwise, returns a function `drive_coef(t, params...; return_ℂ = return_ℂ, kwargs...)` that computes the drive coefficient at time `t`.

# Inner Function `drive_coef`
- `t`: The time at which to evaluate the drive coefficient.
- `params...`: Additional parameters.
- `return_ℂ`: A boolean flag indicating whether to return the complex exponential form.
- `kwargs...`: Additional keyword arguments.


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
    Get_Drive_Coef(ν::T1, ε; envelope = Envelopes.Square_Envelope, drive_time = 0, return_ℂ = false) where T1<:AbstractArray

Creates a drive coefficient function based on the provided parameters.

# Arguments
- `ν::T1`: An array of frequencies.
- `ε`: A scalar coefficient.
- `envelope`: A function that defines the envelope shape of the drive signal. Defaults to `Envelopes.Square_Envelope`.
- `drive_time`: A scalar representing the drive time. Defaults to `0`.
- `return_ℂ`: A boolean flag indicating whether to return the complex exponential form. Defaults to `false`.

# Returns
- A function `drive_coef(t, params...; return_ℂ = return_ℂ, kwargs...)` that computes the drive coefficient at time `t`.

# Inner Function `drive_coef`
- `t`: The time at which to evaluate the drive coefficient.
- `params...`: Additional parameters.
- `return_ℂ`: A boolean flag indicating whether to return the complex exponential form.
- `kwargs...`: Additional keyword arguments.

"""
function Get_Drive_Coef(ν::T1,
    ε;
    envelope = Envelopes.Square_Envelope,
    drive_time = 0,
    return_ℂ = false
    )where T1<:AbstractArray
    if typeof(ε)<:Number
        ε = [ε for i in 1:length(ν)]
    end
    
    function drive_coef(t, params...; return_ℂ = return_ℂ, kwargs...)
        if return_ℂ
            return sum(2*π.*ε.*envelope(t).*exp.(-1im.*(2π.*ν.*t)))
        else 
            return sum(2*π*ε.*envelope(t).*sin.(2π.*ν.*t))
        end
    end

    return drive_coef
end



"""
    Get_Low_Pass_Filtered_Drive_Coef(drive_coef, drive_time; freq_cutoff=4, poles=20, fs=1e2, re_or_im=real, return_T=false)

Generates a function to apply a low-pass filter to the output of a given drive coefficient function.

# Arguments
- `drive_coef`: A callable that computes the drive coefficient as a function of time.
- `drive_time`: The total duration of the drive signal in seconds.

# Keyword Arguments
- `freq_cutoff=4`: The cutoff frequency of the low-pass filter (GHz).
- `poles=20`: The number of poles for the filter
- `fs=1e2`: The sampling frequency used to discretize the drive signal.
- `re_or_im=real`: A callable that extracts either the real or imaginary part of the complex signal. Default is `real`.
- `return_T=false`: If true, returns the low-pass filter function T(ν) instead of the filtered drive coefficient.

# Returns
A function that computes the filtered drive coefficient based on the provided parameters. If `return_T` is true,
it returns a function for the low-pass filter T(ν).

"""
function Get_Low_Pass_Filtered_Drive_Coef(drive_coef, drive_time; freq_cutoff = 4, poles = 20, fs = 1e2, re_or_im = real, return_T = false)
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
    Get_Ĥ_D(op::qt.QuantumObject, drive_coef::Union{Nothing, Function}; TDOS=true, params=nothing, init_time=0.0)

Generates a time-dependent Hamiltonian operator for a quantum object driven by an external field.

# Arguments
- `op::qt.QuantumObject`: The quantum operator to which the driving term will be added.
- `drive_coef::Union{Nothing, Function}`: A callable that computes the drive coefficient as a function of time. If `nothing`, defaults to a constant drive coefficient of 1.

# Keyword Arguments
- `TDOS=true`: If true, returns a time-dependent operator sum with two terms representing the driving Hamiltonian in and out. If false, returns a single term representing the driven Hamiltonian.
- `params=nothing`: Any additional parameters required by the drive_coef function.
- `init_time=0.0`: The initial time at which the system starts being driven.

# Returns
A Time depended operator sum object from qt.TimeDependentOperatorSum

"""
function Get_Ĥ_D(op::qt.QuantumObject, drive_coef::Union{Nothing, Function}; TDOS = true, params = nothing, init_time = 0.0)
    drive_coef = (drive_coef isa Nothing) ? (t)->1 : drive_coef
    if TDOS
        drive_coef_half(t,params...) = drive_coef(t, params...)/2
        drive_coef_half_dag(t, params...) = conj(drive_coef(t, params...))/2
        return qt.TimeDependentOperatorSum([drive_coef_half, drive_coef_half_dag], [op, op'], params = params, init_time = init_time)
    else
        return t -> 0.5*(drive_coef(t)*op+(drive_coef(t)*op)')
    end
end




"""
    Get_Lₜ(op::qt.QuantumObject, drive_coef::Function; params=nothing, init_time=0.0)

Generates a time-dependent Lindblad superoperator for a quantum object driven by an external field.

# Arguments
- `op::qt.QuantumObject`: The quantum operator to which the driving term will be added.
- `drive_coef::Function`: A callable that computes the drive coefficient as a function of time.

# Keyword Arguments
- `params=nothing`: Any additional parameters required by the drive_coef function.
- `init_time=0.0`: The initial time at which the system starts being driven.

# Returns
A Time depended operator sum object from qt.TimeDependentOperatorSum

"""
function Get_Lₜ(op, drive_coef; params = nothing, init_time = 0.0)
    drive_coef_half(t,params...) = drive_coef(t, params...)/2
    drive_coef_half_dag(t, params...) = conj(drive_coef(t, params...))/2
    return qt.TimeDependentOperatorSum([drive_coef_half, drive_coef_half_dag], [qt.liouvillian(op), qt.liouvillian(op')], params = params, init_time = init_time)
end

"""
    Propagator(hilbertspace::qt.HilbertSpace, Ĥₜ::Function, tf; progress_meter=false, ti=0)

Computes the propagator matrix for a quantum system driven by a time-dependent Hamiltonian.

# Arguments
- `hilbertspace::qt.HilbertSpace`: The Hilbert space of the quantum system.
- `Ĥₜ::Function`: A callable that computes the time-dependent Hamiltonian as a function of time.

# Keyword Arguments
- `progress_meter=false`: Whether to display a progress meter during computation.
- `ti=0`: The initial time at which the propagator starts being computed.

# Returns
The propagator matrix U, where U is an operator such that ψ(t_f) = U * ψ(t_i), with ψ being the state vector of the system.

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




