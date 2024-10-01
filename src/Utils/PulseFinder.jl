import QuantumToolbox as qt
using YAXArrays
using DimensionalData
import LsqFit as LF
import CairoMakie as cm

function FindStarkShift(model, drive_op, state1, state2, ε, starkshift_list; make_plot = true)
    ν = model.dressed_energies[state2]- model.dressed_energies[state1]

    νs = ν .+ starkshift_list

    arg_list = []

    for i in 1:length(νs)
        arg_dict = Dict{Any, Any}()
        arg_dict["ε"] = ε
        arg_dict["ν"] = νs[i]
        push!(arg_list, arg_dict)
    end

    floq_sweep_res = Floquet_0_Sweep(model, drive_op, arg_list)

    states_to_track = Dict{Any, Any}()

    states_to_track[tostr(state1)] = model.dressed_states[state1];
    states_to_track[tostr(state2)] = model.dressed_states[state2];

    other_sorts = Dict("F_Energies" => floq_sweep_res["F_Energies"])
    tracking_res = State_Tracker(floq_sweep_res["F_Modes"], states_to_track, other_sorts = other_sorts);

    ys = []
    for state in dims(tracking_res, :State)
        ytemp = []
        for step in dims(tracking_res, :Step)
            val = Real(tracking_res[State = At(state), Step = At(step)]["F_Energies"]/pi)
            if val < 0
                val += 2*abs(νs[step])
            end
            push!(ytemp, val)
        end
        push!(ys, ytemp)
    end

    x = collect(νs .- ν)
    difs1 = abs.(ys[1]-ys[2]);
    difs2 = 2*abs.(νs) - abs.(ys[1]-ys[2])

    difs = []
    for i in 1:length(difs1)
        push!(difs, min(difs1[i], difs2[i]))
    end

    to_fit(t, p) = p[3].*sqrt.((t.-p[1]).^2 .+ p[2].^2)
    p0 = zeros(3)
    p0[1] = x[argmin(difs)]
    p0[2] = minimum(difs)
    p0[3] = abs((maximum(difs)-minimum(difs))/(x[argmax(difs)]-x[argmin(difs)]))
    fit = LF.curve_fit(to_fit, x, difs, p0)

    @info "Fit Stuff: "*tostr(fit.param)
    
    if make_plot 
        f = cm.Figure(size = (800, 500), px_per_unit = 3)

        ax1 = cm.Axis(f[1,1], title = "Floquet Energies", xlabel = "Stark Shifts", ylabel = "Quasienergies (GHz)")
        ax2 = cm.Axis(f[2,1], title = "Difference", xlabel = "Stark Shifts", ylabel = "Dif (GHz)")


        x = Real.(collect(νs.-ν))
        y = []

        colorlist = [:forestgreen, :coral]
        markers = ['+', '×']

        for i in 1:length(dims(tracking_res, :State))
            state = dims(tracking_res, :State)[i]
            cm.scatterlines!(ax1, x, ys[i], label = state, color = colorlist[i], marker = markers[i], linewidth = 0.5, markersize = 20)
            #cm.lines!(ax1, x, y, color = colorlist[i], linewidth = 0.5)
        end
        cm.axislegend(ax1)


        x2 = collect(LinRange(x[1], x[end], 101))
        y2 = to_fit(x2, fit.param)

        fitted_shift = round(fit.param[1], sigdigits = 3)
        cm.lines!(ax2, x2, y2, color = :forestgreen, alpha = 0.25, linewidth= 8, label = "Fitted Stark Shift: $fitted_shift GHz")
        cm.scatterlines!(ax2, x, difs, label = "Difs", marker = '+', markersize = 20, color = :black, linewidth = 0.5)
        cm.axislegend(ax2)
        cm.display(f)
    end
    
    return [fit.param[1], 1/(fit.param[2]*fit.param[3])]
    

end

function OptimizePulse(model, 
    ψ1,
    ψ2,
    ε,
    freq_d,
    stark_shift,
    t_range,
    envelope,
    envelope_args;
    levels = 5,
    samples_per_level = 5,
    solver_kwargs = Dict{Any, Any}(),
    spns = 5,
    envelope_params = Dict{Any, Any}(),
    check_op = nothing)
    
    drive_args = Dict{Any, Any}("pulse_time" => 0.0, "epsilon" => ε, "Envelope" => envelope, "shift"=>stark_shift, "freq_d"=>freq_d)
    drive_args["Envelope Args"] = envelope_args
    
    ti = t_range[1]
    tf = t_range[2]
    tspan = LinRange(ti, tf, samples_per_level)
    level_res = []
    list_of_drive_args = []
    for level in 1:levels
        @info "On Level $level"
        tspan = LinRange(ti, tf, samples_per_level)

        level_res = []
        list_of_drive_args = []
        for i in 1:length(tspan)
            t = tspan[i]
            @info "On Step $i: t = $t"
            drive_args["pulse_time"] = t
            # if "pulse_time" in keys(drive_args["Envelope Args"])
            #     drive_args["Envelope Args"]["pulse_time"] = t
            # end
            # if "sigma" in keys(drive_args["Envelope Args"])
            #     drive_args["Envelope Args"]["sigma"] = t/sigma_factor
            # end
            # if "mu" in keys(drive_args["Envelope Args"])
            #     drive_args["Envelope Args"]["mu"] = t/2
            # end

            drive_args["Envelope Args"] = Envelope_Dict_Cal[drive_args["Envelope"]](t, drive_args["Envelope Args"], envelope_params)

            
            run_res = RunSingleOperator(model, ψ1, drive_args, to_return = "Last", save_step = false, solver_kwargs = solver_kwargs, spns = spns, step_name = "Level_"*string(level)*"_step_"*string(i))
            
            goodness = abs(run_res'*ψ2)^2

            if typeof(check_op)<:qt.QuantumObject
                goodness = abs(run_res'*check_op*run_res)
            end

            drive_args["accuracy"] = goodness
            push!(list_of_drive_args, deepcopy(drive_args))
            push!(level_res, goodness)
            @info "Accurace: $goodness"
            @info "------------------------------------------------------------------------------"
        end
        
        max_loc = argmax(level_res)
        
        @info "Highest Accuracy: $(level_res[max_loc])"
        if max_loc == 1
            ti = tspan[1]
        else
            ti = tspan[max_loc - 1]
        end
        
        if max_loc == length(tspan)
            tf = tspan[end]
        else
            tf = tspan[max_loc+1]
        end
        @info "New ti: $ti, new tf: $tf"
        @info "===================================================================================================\n"
    end

    #opt_drive_time = tspan[argmax(level_res)]

    #drive_args["pulse_time"] = opt_drive_time
    #if "pulse_time" in keys(drive_args["Envelope Args"])
    #    drive_args["Envelope Args"]["pulse_time"] = opt_drive_time
    #end

    return list_of_drive_args[argmax(level_res)]

end