# STAGE MAXIMIZATION PROBLEM FORMULATION

function Problem_OCRA_2(InputParameters::InputParam, SolverParameters::SolverParam, Battery::BatteryParam, NSteps_sim, vector_stages_index, iScen, vector_prices, vector_downtime_stages)       

    @unpack (MIPGap, MIPFocus, Method, Cuts, Heuristics) = SolverParameters;
  
    @unpack (NYears, NMonths, NStages, Big, NHoursStep, conv, bin) = InputParameters;     #NSteps,NHoursStage
    @unpack (min_SOC, max_SOC, min_P, max_P, Eff_charge, Eff_discharge, max_SOH, min_SOH, Nfull) = Battery ;         
    disc = 7

    k = min_SOH/(2*Nfull)
    Small = 0.64

    M = Model(Gurobi.Optimizer)
    set_optimizer_attribute(M, "MIPGap", 0.05)

    # DEFINE VARIABLES

    @variable(M, min_SOC <= soc[iStep=1:NSteps_sim+1] <= max_SOC, base_name = "Energy")                # MWh   energy_Capacity NSteps
    @variable(M, min_SOC^2 <= soc_quad[iStep=1:NSteps_sim+1] <= max_SOC^2, base_name = "Square energy")

    @variable(M, min_P <= charge[iStep=1:NSteps_sim] <= max_P, base_name= "Charge")      #max_disc   0<=discharge<=1
    @variable(M, min_P <= discharge[iStep=1:NSteps_sim] <= max_P, base_name= "Discharge")
    
    @variable(M, 0 <= deg[iStep=1:NSteps_sim] <= Small, base_name = "Degradation")

    @variable(M, 0 <= revamping[iStage=1:NStages] <= (max_SOH-min_SOH), base_name = "Revamping")
    @variable(M, min_SOH <= capacity[iStep=1:NSteps_sim+1] <= max_SOH, base_name = "Energy_Capacity")        #energy_Capacity     [iStage=1:NStages]
    @variable(M, e[iStage=1:NStages], Bin, base_name = "Binary Revamp")
    
    @variable(M, 0<= rev_vendita[iStage=1:NStages] <= max_SOH, base_name = "Vendita rev")
    @variable(M, -max_SOH<= rev_acquisto[iStage=1:NStages] <= 0, base_name = "Acquisto rev")

    #VARIABLES FOR DISCRETIZATION of Stored Energy

    @variable(M, x[iStep=1:NSteps_sim+1], Bin, base_name = "Binary_1")
    @variable(M, y[iStep=1:NSteps_sim+1], Bin, base_name = "Binary_2")
    @variable(M, z[iStep=1:NSteps_sim+1], Bin, base_name = "Binary_3")
    
    @variable(M, 0<= w_xx[iStep=1:NSteps_sim+1] <= 1, base_name = "xx")
    @variable(M, 0<= w_yy[iStep=1:NSteps_sim+1] <= 1, base_name = "yy")
    @variable(M, 0<= w_zz[iStep=1:NSteps_sim+1] <= 1, base_name = "zz")
    @variable(M, 0<= w_xy[iStep=1:NSteps_sim+1] <= 1, base_name = "xy")
    @variable(M, 0<= w_xz[iStep=1:NSteps_sim+1] <= 1, base_name = "xz")
    @variable(M, 0<= w_zy[iStep=1:NSteps_sim+1] <= 1, base_name = "yz")

    # DEFINE OBJECTIVE function - length(Battery_price) = NStages+1=21
    @objective(
      M,
      MathOptInterface.MAX_SENSE, 
      sum(vector_prices[iScen][iStep]*NHoursStep*(discharge[iStep]-charge[iStep]) for iStep=1:NSteps_sim) -               #-sum(Battery_price[iStage]*(revamping[iStage]) for iStage=1:NStages) + 
      Battery_price_purchase[1]*(revamping[1]) -
      sum(Battery_price_purchase[iStage]*(capacity[vector_stages_index[iStage]+2] + rev_acquisto[iStage]) for iStage=2:NStages) +
      sum(Battery_price_sale[iStage]*(capacity[vector_stages_index[iStage]+1] - rev_vendita[iStage]) for iStage=2:NStages) +
      Battery_price_sale[NStages+1]*(capacity[end]- min_SOH) - 
      sum(fix*e[iStage] for iStage=1:NStages) + 2300    
      )
         
    # DEFINE CONSTRAINTS

    #@constraint(M, Charge_op[iStep=1:NSteps_sim], charge[iStep] <= max_P*e[iStep])
    #@constraint(M, Disch_op[iStep=1:NSteps_sim], discharge[iStep] <= max_P*(1-e[iStep]))

    @constraint(M,energy[iStep=1:NSteps_sim], soc[iStep] + (charge[iStep]*Eff_charge-discharge[iStep]/Eff_discharge)*NHoursStep == soc[iStep+1] )

    @constraint(M, en_bal[iStep=1:NSteps_sim+1], min_SOC + ((max_SOC-min_SOC)/disc)*(x[iStep]+2*y[iStep]+4*z[iStep]) == soc[iStep])
    @constraint(M, en_square[iStep=1:NSteps_sim+1], soc_quad[iStep] == min_SOC^2+ 2*min_SOC*((max_SOC-min_SOC)/disc)*(x[iStep]+2*y[iStep]+4*z[iStep])+(w_xx[iStep]+4*w_xy[iStep]+8*w_xz[iStep]+4*w_yy[iStep]+16*w_zz[iStep]+16*w_zy[iStep])*((max_SOC-min_SOC)/disc)^2)

    # INEQUALITIES CONSTRAINTS
    @constraint(M, xx_1[iStep=1:NSteps_sim+1], w_xx[iStep] <= x[iStep])
    @constraint(M, xx_2[iStep=1:NSteps_sim+1], w_xx[iStep] >= 2*x[iStep]-1)

    @constraint(M, xy_1[iStep=1:NSteps_sim+1], w_xy[iStep] <= x[iStep])
    @constraint(M, xy_2[iStep=1:NSteps_sim+1], w_xy[iStep] <= y[iStep])
    @constraint(M, xy_3[iStep=1:NSteps_sim+1], w_xy[iStep] >= x[iStep]+y[iStep]-1)

    @constraint(M, xz_1[iStep=1:NSteps_sim+1], w_xz[iStep] <= x[iStep])
    @constraint(M, xz_2[iStep=1:NSteps_sim+1], w_xz[iStep] <= z[iStep])
    @constraint(M, xz_3[iStep=1:NSteps_sim+1], w_xz[iStep] >= x[iStep]+z[iStep]-1)

    @constraint(M, yy_1[iStep=1:NSteps_sim+1], w_yy[iStep] <= y[iStep])
    @constraint(M, yy_2[iStep=1:NSteps_sim+1], w_yy[iStep] >= 2*y[iStep]-1)

    @constraint(M, zz_1[iStep=1:NSteps_sim+1], w_zz[iStep] <= z[iStep])
    @constraint(M, zz_2[iStep=1:NSteps_sim+1], w_zz[iStep] >= 2*z[iStep]-1)

    @constraint(M, zy_1[iStep=1:NSteps_sim+1], w_zy[iStep] <= z[iStep])
    @constraint(M, zy_2[iStep=1:NSteps_sim+1], w_zy[iStep] <= y[iStep])
    @constraint(M, zy_3[iStep=1:NSteps_sim+1], w_zy[iStep] >= z[iStep]+y[iStep]-1) 

    # CONSTRAINTS ON DEGRADATION
    @constraint(M, deg_1[iStep=1:NSteps_sim], deg[iStep] >= soc_quad[iStep]/max_SOC^2 - soc_quad[iStep+1]/max_SOC^2 + (2/max_SOC)*(soc[iStep+1]-soc[iStep]))
    @constraint(M, deg_2[iStep=1:NSteps_sim], deg[iStep] >= soc_quad[iStep+1]/max_SOC^2 - soc_quad[iStep]/max_SOC^2 + (2/max_SOC)*(soc[iStep]-soc[iStep+1]))

    #CONSTRAINT ON REVAMPING
    @constraint(M, energy_capacity[iStage=1:NStages], capacity[vector_stages_index[iScen, iStage]+2] == capacity[vector_stages_index[iScen,iStage]+1]-deg[vector_stages_index[iScen, iStage]+1]*k+revamping[iStage])
   
    @constraint(M, initial_e[iStep=1], capacity[iStep] == min_SOH) #max_SOH

    @constraint(M,en_cap[iStage in 1:NStages, iStep in ((vector_stages_index[iScen,iStage]+2):vector_stages_index[iScen,iStage+1])], capacity[iStep+1]== capacity[iStep]-deg[iStep]*k) 

    @constraint(M, stop_charge[iStage in 2:NStages, iStep in (vector_stages_index[iScen,iStage]:(vector_stages_index[iScen,iStage]+vector_downtime_stages[iScen,iStage-1]))], charge[iStep] <= (1-e[iStage])*max_P)
    
    @constraint(M, stop_discharge[iStage in 2:NStages, iStep in (vector_stages_index[iScen,iStage]:(vector_stages_index[iScen,iStage]+vector_downtime_stages[iScen, iStage-1]))], discharge[iStep] <= (1-e[iStage])*max_P) 

    @constraint(M, rev[iStage=1:NStages], revamping[iStage] <= (max_SOH-min_SOH)*e[iStage])

    # CONSTRAINT SU VARIABILI AUSILIARIE

    @constraint(M, vendita[iStage=1], rev_vendita[iStage] == 0)
    @constraint(M, vendita_1[iStage=2:NStages], rev_vendita[iStage] >= 0)
    @constraint(M, vendita_2[iStage=2:NStages], rev_vendita[iStage] >= capacity[vector_stages_index[iScen,iStage]+1]- e[iStage]*Big)

    @constraint(M, acquisto_1[iStage=2:NStages], rev_acquisto[iStage] >= -(1-e[iStage])*Big)
    @constraint(M, acquisto_2[iStage=2:NStages], rev_acquisto[iStage] >= -capacity[vector_stages_index[iScen,iStage]+2] )
    @constraint(M, acquisto_3[iStage=1], rev_acquisto[iStage] == 0)

    return Problem_OCRA2(
        M,
        soc,
        soc_quad,
        charge,
        discharge,
        deg,
        x,
        y,
        z,
        w_xx,
        w_yy,
        w_zz,
        w_xy,
        w_xz,
        w_zy,
        capacity,
        revamping,
        e,
        rev_vendita,
        rev_acquisto,
      )
end

