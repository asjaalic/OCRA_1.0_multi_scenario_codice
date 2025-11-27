# SOLVE OPTIMIZATION PROBLEM

function solveWithoutDeg_3(InputParameters::InputParam, SolverParameters::SolverParam, Battery::BatteryParam, NScen, NSteps, Pp, Steps_stages)

    @unpack (NYears, NMonths, NStages, Big, NHoursStep, bin) = InputParameters;                #NSteps, NHoursStage
    @unpack (min_SOC, max_SOC, Eff_charge, Eff_discharge, min_P, max_P, max_SOH, min_SOH, Nfull) = Battery;

    println("Solving Optimization Problem Without Degradation")

    objective_no_deg = zeros(NScen)                   
    revenues_per_stage_no_deg = zeros(NScen, NStages)

    charge_no_deg = zeros(NScen, NSteps)
    discharge_no_deg = zeros(NScen, NSteps)
    soc_no_deg = zeros(NScen, NSteps+1)
    soc_quad_no_deg = zeros(NScen, NSteps+1)
    e_no_deg = zeros(NScen, NStages)

    x_no_deg = zeros(NScen, NSteps+1)
    y_no_deg = zeros(NScen, NSteps+1)
    z_no_deg = zeros(NScen, NSteps+1)

    w_xx_no_deg = zeros(NScen, NSteps+1)
    w_yy_no_deg = zeros(NScen, NSteps+1)
    w_zz_no_deg = zeros(NScen, NSteps+1)
    w_xy_no_deg = zeros(NScen, NSteps+1)
    w_xz_no_deg = zeros(NScen, NSteps+1)
    w_zy_no_deg = zeros(NScen, NSteps+1)

    sim = BuildStageProblemNoDeg_3(InputParameters, SolverParameters, Battery)

    for iScen=1:NScen

        println("Solving scenario $iScen")

        for iStep=1:NSteps      #Update power prices for each scenario
            set_objective_coefficient(sim.M_sim,
            sim.discharge[iStep],
            NHoursStep*Pp[iStep, iScen])

            set_objective_coefficient(sim.M_sim,
            sim.charge[iStep],
            -NHoursStep*Pp[iStep, iScen])
        end

        @timeit to "Solve optimization" optimize!(sim.M_sim)

        if termination_status(sim.M_sim) != MOI.OPTIMAL
            println("NOT OPTIMAL: ", termination_status(sim.M_sim))
        else
            println("Optimization finished")
        end

        @timeit to "Collecting results" begin
            objective_no_deg[iScen] = JuMP.objective_value(sim.M_sim)
            
            for iStep=1:NSteps
                soc_no_deg[iScen, iStep] = JuMP.value(sim.soc[iStep])
                charge_no_deg[iScen, iStep] = JuMP.value(sim.charge[iStep])
                discharge_no_deg[iScen, iStep] = JuMP.value(sim.discharge[iStep])
                soc_quad_no_deg[iScen, iStep] = JuMP.value(sim.soc_quad[iStep])

                x_no_deg[iScen, iStep] = JuMP.value(sim.x[iStep])
                y_no_deg[iScen, iStep] = JuMP.value(sim.y[iStep])
                z_no_deg[iScen, iStep] = JuMP.value(sim.z[iStep])

                w_xx_no_deg[iScen, iStep] = JuMP.value(sim.w_xx[iStep])
                w_yy_no_deg[iScen, iStep] = JuMP.value(sim.w_yy[iStep])
                w_zz_no_deg[iScen, iStep] = JuMP.value(sim.w_zz[iStep])
                w_xy_no_deg[iScen, iStep] = JuMP.value(sim.w_xy[iStep])
                w_xz_no_deg[iScen, iStep] = JuMP.value(sim.w_xz[iStep])
                w_zy_no_deg[iScen, iStep] = JuMP.value(sim.w_zy[iStep])

            end

            soc_no_deg[iScen, end] = JuMP.value(sim.soc[end])
            soc_quad_no_deg[iScen, end] = JuMP.value(sim.soc_quad[end])
            x_no_deg[iScen, end] = JuMP.value(sim.x[end])
            y_no_deg[iScen, end] = JuMP.value(sim.y[end])
            z_no_deg[iScen, end] = JuMP.value(sim.z[end])

            w_xx_no_deg[iScen, end] = JuMP.value(sim.w_xx[end])
            w_yy_no_deg[iScen, end] = JuMP.value(sim.w_yy[end])
            w_zz_no_deg[iScen, end] = JuMP.value(sim.w_zz[end])
            w_xy_no_deg[iScen, end] = JuMP.value(sim.w_xy[end])
            w_xz_no_deg[iScen, end] = JuMP.value(sim.w_xz[end])
            w_zy_no_deg[iScen, end] = JuMP.value(sim.w_zy[end])

        
            for iStage=2:(NStages-1)
                revenues_per_stage_no_deg[iScen, iStage] = sum(Pp[iStep, iScen]*NHoursStep*(discharge_no_deg[iScen, iStep]-charge_no_deg[iScen, iStep]) for iStep=(Steps_stages[iStage]+1):(Steps_stages[iStage+1])) 
            end
            revenues_per_stage_no_deg[iScen, 1] = sum(Pp[iStep, iScen]*NHoursStep*(discharge_no_deg[iScen, iStep] - charge_no_deg[iScen, iStep]) for iStep=(Steps_stages[1]+1):(Steps_stages[2])) 
            revenues_per_stage_no_deg[iScen, NStages] = sum(Pp[iStep, iScen]*NHoursStep*(discharge_no_deg[iScen, iStep] - charge_no_deg[iScen,iStep]) for iStep=(Steps_stages[NStages]+1):(Steps_stages[NStages+1])) 
            
        end
    end
    
    println("Collected results")

    return ResultWithoutDeg_3(
        objective_no_deg,
        revenues_per_stage_no_deg,
        soc_no_deg,
        charge_no_deg,
        discharge_no_deg,
        soc_quad_no_deg,
        x_no_deg,
        y_no_deg,
        z_no_deg,
        w_xx_no_deg,
        w_yy_no_deg,
        w_zz_no_deg,
        w_xy_no_deg,
        w_xz_no_deg,
        w_zy_no_deg,
        e_no_deg,
    )

end


function solveWithoutDeg_4(InputParameters::InputParam, SolverParameters::SolverParam, Battery::BatteryParam, NScen, NSteps, Pp, Steps_stages)

    @unpack (NYears, NMonths, NStages, Big, NHoursStep, bin) = InputParameters;                #NSteps, NHoursStage
    @unpack (min_SOC, max_SOC, Eff_charge, Eff_discharge, min_P, max_P, max_SOH, min_SOH, Nfull) = Battery;

    println("Solving Optimization Problem Without Degradation")

    objective_no_deg = zeros(NScen)                   
    revenues_per_stage_no_deg = zeros(NScen, NStages)

    charge_no_deg = zeros(NScen, NSteps)
    discharge_no_deg = zeros(NScen, NSteps)
    soc_no_deg = zeros(NScen, NSteps+1)
    soc_quad_no_deg = zeros(NScen, NSteps+1)
    e_no_deg = zeros(NScen, NStages)

    x_no_deg = zeros(NScen, NSteps+1)
    y_no_deg = zeros(NScen, NSteps+1)
    z_no_deg = zeros(NScen, NSteps+1)

    w_xx_no_deg = zeros(NScen, NSteps+1)
    w_yy_no_deg = zeros(NScen, NSteps+1)
    w_zz_no_deg = zeros(NScen, NSteps+1)
    w_xy_no_deg = zeros(NScen, NSteps+1)
    w_xz_no_deg = zeros(NScen, NSteps+1)
    w_zy_no_deg = zeros(NScen, NSteps+1)


    sim = BuildStageProblemNoDeg(InputParameters, SolverParameters, Battery)

    for iScen=1:NScen

        println("Solving scenario $iScen")

        for iStep=1:NSteps      #Update power prices for each scenario
            set_objective_coefficient(sim.M_sim,
            sim.discharge[iStep],
            NHoursStep*Pp[iStep, iScen])

            set_objective_coefficient(sim.M_sim,
            sim.charge[iStep],
            -NHoursStep*Pp[iStep, iScen])
        end

        @timeit to "Solve optimization" optimize!(sim.M_sim)

        if termination_status(sim.M_sim) != MOI.OPTIMAL
            println("NOT OPTIMAL: ", termination_status(sim.M_sim))
        else
            println("Optimization finished")
        end

        @timeit to "Collecting results" begin
            objective_no_deg[iScen] = JuMP.objective_value(sim.M_sim)
            
            for iStep=1:NSteps
                soc_no_deg[iScen, iStep] = JuMP.value(sim.soc[iStep])
                charge_no_deg[iScen, iStep] = JuMP.value(sim.charge[iStep])
                discharge_no_deg[iScen, iStep] = JuMP.value(sim.discharge[iStep])
                soc_quad_no_deg[iScen, iStep] = JuMP.value(sim.soc_quad[iStep])

                x_no_deg[iScen, iStep] = JuMP.value(sim.x[iStep])
                y_no_deg[iScen, iStep] = JuMP.value(sim.y[iStep])
                z_no_deg[iScen, iStep] = JuMP.value(sim.z[iStep])
                u_no_deg[iScen, iStep] = JuMP.value(sim.u[iStep])

                w_xx_no_deg[iScen, iStep] = JuMP.value(sim.w_xx[iStep])
                w_yy_no_deg[iScen, iStep] = JuMP.value(sim.w_yy[iStep])
                w_zz_no_deg[iScen, iStep] = JuMP.value(sim.w_zz[iStep])
                w_uu_no_deg[iScen, iStep] = JuMP.value(sim.w_uu[iStep])
                
                w_xy_no_deg[iScen, iStep] = JuMP.value(sim.w_xy[iStep])
                w_xz_no_deg[iScen, iStep] = JuMP.value(sim.w_xz[iStep])
                w_zy_no_deg[iScen, iStep] = JuMP.value(sim.w_zy[iStep])

            end

            soc_no_deg[iScen, end] = JuMP.value(sim.soc[end])
            soc_quad_no_deg[iScen, end] = JuMP.value(sim.soc_quad[end])
            x_no_deg[iScen, end] = JuMP.value(sim.x[end])
            y_no_deg[iScen, end] = JuMP.value(sim.y[end])
            z_no_deg[iScen, end] = JuMP.value(sim.z[end])

            w_xx_no_deg[iScen, end] = JuMP.value(sim.w_xx[end])
            w_yy_no_deg[iScen, end] = JuMP.value(sim.w_yy[end])
            w_zz_no_deg[iScen, end] = JuMP.value(sim.w_zz[end])
            w_xy_no_deg[iScen, end] = JuMP.value(sim.w_xy[end])
            w_xz_no_deg[iScen, end] = JuMP.value(sim.w_xz[end])
            w_zy_no_deg[iScen, end] = JuMP.value(sim.w_zy[end])

        
            for iStage=2:(NStages-1)
                revenues_per_stage_no_deg[iScen, iStage] = sum(Pp[iStep, iScen]*NHoursStep*(discharge_no_deg[iScen, iStep]-charge_no_deg[iScen, iStep]) for iStep=(Steps_stages[iStage]+1):(Steps_stages[iStage+1])) 
            end
            revenues_per_stage_no_deg[iScen, 1] = sum(Pp[iStep, iScen]*NHoursStep*(discharge_no_deg[iScen, iStep] - charge_no_deg[iScen, iStep]) for iStep=(Steps_stages[1]+1):(Steps_stages[2])) 
            revenues_per_stage_no_deg[iScen, NStages] = sum(Pp[iStep, iScen]*NHoursStep*(discharge_no_deg[iScen, iStep] - charge_no_deg[iScen,iStep]) for iStep=(Steps_stages[NStages]+1):(Steps_stages[NStages+1])) 
            
        end
    end
    
    println("Collected results")

    return ResultWithoutDeg_4(
        objective_no_deg,
        revenues_per_stage_no_deg,
        soc_no_deg,
        charge_no_deg,
        discharge_no_deg,
        soc_quad_no_deg,
        x_no_deg,
        y_no_deg,
        z_no_deg,
        u_no_deg,
        w_xx_no_deg,
        w_yy_no_deg,
        w_zz_no_deg,
        w_uu_no_deg,
        w_xy_no_deg,
        w_xz_no_deg,
        w_zy_no_deg,
        w_xu_no_deg,
        w_yu_no_deg,
        w_zu_no_deg,
        e_no_deg,
    )

end


function solveWithoutDeg_5(InputParameters::InputParam, SolverParameters::SolverParam, Battery::BatteryParam, NScen, NSteps, Pp, Steps_stages)

    @unpack (NYears, NMonths, NStages, Big, NHoursStep, bin) = InputParameters;                #NSteps, NHoursStage
    @unpack (min_SOC, max_SOC, Eff_charge, Eff_discharge, min_P, max_P, max_SOH, min_SOH, Nfull) = Battery;

    println("Solving Optimization Problem Without Degradation")

    objective_no_deg = zeros(NScen)                   
    revenues_per_stage_no_deg = zeros(NScen, NStages)

    charge_no_deg = zeros(NScen, NSteps)
    discharge_no_deg = zeros(NScen, NSteps)
    soc_no_deg = zeros(NScen, NSteps+1)
    soc_quad_no_deg = zeros(NScen, NSteps+1)
    e_no_deg = zeros(NScen, NStages)

    x_no_deg = zeros(NScen, NSteps+1)
    y_no_deg = zeros(NScen, NSteps+1)
    z_no_deg = zeros(NScen, NSteps+1)
    u_no_deg = zeros(NScen, NSteps+1)
    t_no_deg = zeros(NScen, NSteps+1)

    w_xx_no_deg = zeros(NScen, NSteps+1)
    w_yy_no_deg = zeros(NScen, NSteps+1)
    w_zz_no_deg = zeros(NScen, NSteps+1)
    w_uu_no_deg = zeros(NScen, NSteps+1)
    w_tt_no_deg = zeros(NScen, NSteps+1)
    w_xy_no_deg = zeros(NScen, NSteps+1)
    w_xz_no_deg = zeros(NScen, NSteps+1)
    w_zy_no_deg = zeros(NScen, NSteps+1)
    w_xu_no_deg = zeros(NScen, NSteps+1)
    w_yu_no_deg = zeros(NScen, NSteps+1)
    w_zu_no_deg = zeros(NScen, NSteps+1)
    w_tx_no_deg = zeros(NScen, NSteps+1)
    w_ty_no_deg = zeros(NScen, NSteps+1)
    w_tz_no_deg = zeros(NScen, NSteps+1)
    w_tu_no_deg = zeros(NScen, NSteps+1)


    sim = BuildStageProblemNoDeg_5(InputParameters, SolverParameters, Battery)

    for iScen=1:NScen

        println("Solving scenario $iScen")

        for iStep=1:NSteps      #Update power prices for each scenario
            set_objective_coefficient(sim.M_sim,
            sim.discharge[iStep],
            NHoursStep*Pp[iStep, iScen])

            set_objective_coefficient(sim.M_sim,
            sim.charge[iStep],
            -NHoursStep*Pp[iStep, iScen])
        end

        @timeit to "Solve optimization" optimize!(sim.M_sim)

        if termination_status(sim.M_sim) != MOI.OPTIMAL
            println("NOT OPTIMAL: ", termination_status(sim.M_sim))
        else
            println("Optimization finished")
        end

        @timeit to "Collecting results" begin
            objective_no_deg[iScen] = JuMP.objective_value(sim.M_sim)
            
            for iStep=1:NSteps
                soc_no_deg[iScen, iStep] = JuMP.value(sim.soc[iStep])
                charge_no_deg[iScen, iStep] = JuMP.value(sim.charge[iStep])
                discharge_no_deg[iScen, iStep] = JuMP.value(sim.discharge[iStep])
                soc_quad_no_deg[iScen, iStep] = JuMP.value(sim.soc_quad[iStep])

                x_no_deg[iScen, iStep] = JuMP.value(sim.x[iStep])
                y_no_deg[iScen, iStep] = JuMP.value(sim.y[iStep])
                z_no_deg[iScen, iStep] = JuMP.value(sim.z[iStep])
                u_no_deg[iScen, iStep] = JuMP.value(sim.u[iStep])
                t_no_deg[iScen, iStep] = JuMP.value(sim.t[iStep])

                w_xx_no_deg[iScen, iStep] = JuMP.value(sim.w_xx[iStep])
                w_yy_no_deg[iScen, iStep] = JuMP.value(sim.w_yy[iStep])
                w_zz_no_deg[iScen, iStep] = JuMP.value(sim.w_zz[iStep])
                w_uu_no_deg[iScen, iStep] = JuMP.value(sim.w_uu[iStep])
                w_tt_no_deg[iScen, iStep] = JuMP.value(sim.w_tt[iStep])

                w_xy_no_deg[iScen, iStep] = JuMP.value(sim.w_xy[iStep])
                w_xz_no_deg[iScen, iStep] = JuMP.value(sim.w_xz[iStep])
                w_zy_no_deg[iScen, iStep] = JuMP.value(sim.w_zy[iStep])
                w_xu_no_deg[iScen, iStep] = JuMP.value(sim.w_xu[iStep])
                w_yu_no_deg[iScen, iStep] = JuMP.value(sim.w_yu[iStep])
                w_zu_no_deg[iScen, iStep] = JuMP.value(sim.w_zu[iStep])
                w_tx_no_deg[iScen, iStep] = JuMP.value(sim.w_tx[iStep])
                w_ty_no_deg[iScen, iStep] = JuMP.value(sim.w_ty[iStep])
                w_tz_no_deg[iScen, iStep] = JuMP.value(sim.w_tz[iStep])
                w_tu_no_deg[iScen, iStep] = JuMP.value(sim.w_tu[iStep])

            end

            soc_no_deg[iScen, end] = JuMP.value(sim.soc[end])
            soc_quad_no_deg[iScen, end] = JuMP.value(sim.soc_quad[end])
            x_no_deg[iScen, end] = JuMP.value(sim.x[end])
            y_no_deg[iScen, end] = JuMP.value(sim.y[end])
            z_no_deg[iScen, end] = JuMP.value(sim.z[end])
            u_no_deg[iScen, end] = JuMP.value(sim.u[end])
            t_no_deg[iScen, end] = JuMP.value(sim.t[end])

            w_xx_no_deg[iScen, end] = JuMP.value(sim.w_xx[end])
            w_yy_no_deg[iScen, end] = JuMP.value(sim.w_yy[end])
            w_zz_no_deg[iScen, end] = JuMP.value(sim.w_zz[end])
            w_uu_no_deg[iScen, end] = JuMP.value(sim.w_uu[end])
            w_tt_no_deg[iScen, end] = JuMP.value(sim.w_tt[end])

            w_xy_no_deg[iScen, end] = JuMP.value(sim.w_xy[end])
            w_xz_no_deg[iScen, end] = JuMP.value(sim.w_xz[end])
            w_zy_no_deg[iScen, end] = JuMP.value(sim.w_zy[end])
            w_xu_no_deg[iScen, end] = JuMP.value(sim.w_xu[end])
            w_yu_no_deg[iScen, end] = JuMP.value(sim.w_yu[end])
            w_zu_no_deg[iScen, end] = JuMP.value(sim.w_zu[end])
            w_tx_no_deg[iScen, end] = JuMP.value(sim.w_tx[end])
            w_ty_no_deg[iScen, end] = JuMP.value(sim.w_ty[end])
            w_tz_no_deg[iScen, end] = JuMP.value(sim.w_tz[end])
            w_tu_no_deg[iScen, end] = JuMP.value(sim.w_tu[end])

        
            for iStage=2:(NStages-1)
                revenues_per_stage_no_deg[iScen, iStage] = sum(Pp[iStep, iScen]*NHoursStep*(discharge_no_deg[iScen, iStep]-charge_no_deg[iScen, iStep]) for iStep=(Steps_stages[iStage]+1):(Steps_stages[iStage+1])) 
            end
            revenues_per_stage_no_deg[iScen, 1] = sum(Pp[iStep, iScen]*NHoursStep*(discharge_no_deg[iScen, iStep] - charge_no_deg[iScen, iStep]) for iStep=(Steps_stages[1]+1):(Steps_stages[2])) 
            revenues_per_stage_no_deg[iScen, NStages] = sum(Pp[iStep, iScen]*NHoursStep*(discharge_no_deg[iScen, iStep] - charge_no_deg[iScen,iStep]) for iStep=(Steps_stages[NStages]+1):(Steps_stages[NStages+1])) 
            
        end
    end
    
    println("Collected results")

    return ResultWithoutDeg_5(
        objective_no_deg,
        revenues_per_stage_no_deg,
        soc_no_deg,
        charge_no_deg,
        discharge_no_deg,
        soc_quad_no_deg,
        x_no_deg,
        y_no_deg,
        z_no_deg,
        u_no_deg,
        w_xx_no_deg,
        w_yy_no_deg,
        w_zz_no_deg,
        w_uu_no_deg,
        w_xy_no_deg,
        w_xz_no_deg,
        w_zy_no_deg,
        w_xu_no_deg,
        w_yu_no_deg,
        w_zu_no_deg,
        e_no_deg,
        t_no_deg,
        w_tt_no_deg,
        w_tx_no_deg,
        w_ty_no_deg,
        w_tz_no_deg,
        w_tu_no_deg,
    )

end