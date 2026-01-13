# STAGE MAXIMIZATION PROBLEM FORMULATION

function BuildStageProblem_3(InputParameters::InputParam, SolverParameters::SolverParam, Battery::BatteryParam, NSteps_sim, vector_stages_index, iScen, vector_prices)       

    @unpack (MIPGap, MIPFocus, Method, Cuts, Heuristics) = SolverParameters;
  
    @unpack (NYears, NMonths, NStages, Big, NHoursStep, conv, bin) = InputParameters;     #NSteps,NHoursStage
    @unpack (min_SOC, max_SOC, min_P, max_P, Eff_charge, Eff_discharge, max_SOH, min_SOH, Nfull) = Battery ;         
    disc = 7

    k = min_SOH/(2*Nfull)
    Small = 0.64

    M = Model(Gurobi.Optimizer)
    set_optimizer_attribute(M, "MIPGap", 0.01)

    # DEFINE VARIABLES

    @variable(M, min_SOC <= soc[iStep=1:NSteps_sim+1] <= max_SOC, base_name = "Energy")                # MWh   energy_Capacity NSteps
    @variable(M, min_SOC^2 <= soc_quad[iStep=1:NSteps_sim+1] <= max_SOC^2, base_name = "Square energy")

    @variable(M, min_P <= charge[iStep=1:NSteps_sim] <= max_P, base_name= "Charge")      #max_disc   0<=discharge<=1
    @variable(M, min_P <= discharge[iStep=1:NSteps_sim] <= max_P, base_name= "Discharge")
    
    @variable(M, 0 <= deg[iStep=1:NSteps_sim] <= Small, base_name = "Degradation")

    @variable(M, 0 <= revamping[iStage=1:NStages] <= (max_SOH-min_SOH), base_name = "Revamping")
    @variable(M, min_SOH <= capacity[iStep=1:NSteps_sim+1] <= max_SOH, base_name = "Energy_Capacity")        #energy_Capacity     [iStage=1:NStages]
    
    @variable(M, e[iStep=1:NSteps_sim], Bin, base_name ="Binary operation")

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
      sum(vector_prices[iScen][iStep]*discharge[iStep]-vector_prices[iScen][iStep]*charge[iStep] for iStep=1:NSteps_sim) 
      -sum(Battery_price_purchase[iStage]*(revamping[iStage]) for iStage=1:NStages)  
      +Battery_price_purchase[NStages+1]*(capacity[end]-min_SOH)+2300   
      )
         
    # DEFINE CONSTRAINTS

    @constraint(M, Charge_op[iStep=1:NSteps_sim], charge[iStep] <= max_P*e[iStep])
    @constraint(M, Disch_op[iStep=1:NSteps_sim], discharge[iStep] <= max_P*(1-e[iStep]))

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

    @constraint(M,en_cap[iStage in 1:NStages, iStep in (vector_stages_index[iScen,iStage]+2:vector_stages_index[iScen,iStage+1])], capacity[iStep+1]== capacity[iStep]-deg[iStep]*k) 
    

    return BuildStageProblem_3(
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
      )
end


function BuildStageProblem_4(InputParameters::InputParam, SolverParameters::SolverParam, Battery::BatteryParam, NSteps_sim, vector_stages_index, iScen, vector_prices)       

    @unpack (MIPGap, MIPFocus, Method, Cuts, Heuristics) = SolverParameters;
  
    @unpack (NYears, NMonths, NStages, Big, NHoursStep, conv, bin) = InputParameters;     #NSteps,NHoursStage
    @unpack (min_SOC, max_SOC, min_P, max_P, Eff_charge, Eff_discharge, max_SOH, min_SOH, Nfull) = Battery ;         
    disc = 15

    k = min_SOH/(2*Nfull)
    Small = 1

    M = Model(Gurobi.Optimizer)
    set_optimizer_attribute(M, "MIPGap", 0.01)

    # DEFINE VARIABLES

    @variable(M, min_SOC <= soc[iStep=1:NSteps_sim+1] <= max_SOC, base_name = "Energy")                # MWh   energy_Capacity NSteps
    @variable(M, min_SOC^2 <= soc_quad[iStep=1:NSteps_sim+1] <= max_SOC^2, base_name = "Square energy")

    @variable(M, min_P <= charge[iStep=1:NSteps_sim] <= max_P, base_name= "Charge")      #max_disc   0<=discharge<=1
    @variable(M, min_P <= discharge[iStep=1:NSteps_sim] <= max_P, base_name= "Discharge")
    
    @variable(M, 0 <= deg[iStep=1:NSteps_sim] <= Small, base_name = "Degradation")

    @variable(M, 0 <= revamping[iStage=1:NStages] <= (max_SOH-min_SOH), base_name = "Revamping")
    @variable(M, min_SOH <= capacity[iStep=1:NSteps_sim+1] <= max_SOH, base_name = "Energy_Capacity")        #energy_Capacity     [iStage=1:NStages]
    
    @variable(M, e[iStep=1:NSteps_sim], Bin, base_name ="Binary operation")

    #VARIABLES FOR DISCRETIZATION of Stored Energy

    @variable(M, x[iStep=1:NSteps_sim+1], Bin, base_name = "Binary_1")
    @variable(M, y[iStep=1:NSteps_sim+1], Bin, base_name = "Binary_2")
    @variable(M, z[iStep=1:NSteps_sim+1], Bin, base_name = "Binary_3")
    @variable(M, u[iStep=1:NSteps_sim+1], Bin, base_name = "Binary_4")
    
    @variable(M, 0<= w_xx[iStep=1:NSteps_sim+1] <= 1, base_name = "xx")
    @variable(M, 0<= w_yy[iStep=1:NSteps_sim+1] <= 1, base_name = "yy")
    @variable(M, 0<= w_zz[iStep=1:NSteps_sim+1] <= 1, base_name = "zz")
    @variable(M, 0<= w_xy[iStep=1:NSteps_sim+1] <= 1, base_name = "xy")
    @variable(M, 0<= w_xz[iStep=1:NSteps_sim+1] <= 1, base_name = "xz")
    @variable(M, 0<= w_zy[iStep=1:NSteps_sim+1] <= 1, base_name = "yz")

  #=  @variable(M, 0 <= w_uu[iStep=1:NSteps_sim+1] <=1, base_name = "uu")
    @variable(M, 0 <= w_xu[iStep=1:NSteps_sim+1] <=1, base_name = "xu")
    @variable(M, 0 <= w_yu[iStep=1:NSteps_sim+1] <=1, base_name = "yu")
    @variable(M, 0 <= w_zu[iStep=1:NSteps_sim+1] <=1, base_name = "zu")=#

    # DEFINE OBJECTIVE function - length(Battery_price) = NStages+1=21

    @objective(
      M,
      MathOptInterface.MAX_SENSE, 
      sum(vector_prices[iScen][iStep]*discharge[iStep]-vector_prices[iScen][iStep]*charge[iStep] for iStep=1:NSteps_sim) 
      -sum(Battery_price_purchase[iStage]*(revamping[iStage]) for iStage=1:NStages)  
      +Battery_price_purchase[NStages+1]*(capacity[end]-min_SOH)+2300   
      )
         
    # DEFINE CONSTRAINTS

    @constraint(M, Charge_op[iStep=1:NSteps_sim], charge[iStep] <= max_P*e[iStep])
    @constraint(M, Disch_op[iStep=1:NSteps_sim], discharge[iStep] <= max_P*(1-e[iStep]))

    @constraint(M,energy[iStep=1:NSteps_sim], soc[iStep] + (charge[iStep]*Eff_charge-discharge[iStep]/Eff_discharge)*NHoursStep == soc[iStep+1] )

    @constraint(M, en_bal[iStep=1:NSteps_sim+1], min_SOC + ((max_SOC-min_SOC)/disc)*(x[iStep]+2*y[iStep]+4*z[iStep]+8*u[iStep]) == soc[iStep])    
    @constraint(M, en_square[iStep=1:NSteps_sim+1], soc_quad[iStep] == min_SOC^2+ 2*min_SOC*((max_SOC-min_SOC)/disc)*(x[iStep]+2*y[iStep]+4*z[iStep]+8*u[iStep])+(w_xx[iStep]+4*w_xy[iStep]+8*w_xz[iStep]+16*w_xu[iStep]+4*w_yy[iStep]+16*w_zz[iStep]+16*w_zy[iStep]+32*w_yu[iStep]+64*w_zu[iStep]+64*w_uu[iStep])*((max_SOC-min_SOC)/disc)^2)

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

  #=  @constraint(M, uu_1[iStep=1:NSteps_sim+1], w_uu[iStep] <= u[iStep])
    @constraint(M, uu_2[iStep=1:NSteps_sim+1], w_uu[iStep] >= 2*u[iStep]-1)

    @constraint(M, xu_1[iStep=1:NSteps_sim+1], w_xu[iStep] <= x[iStep])
    @constraint(M, xu_2[iStep=1:NSteps_sim+1], w_xu[iStep] <= u[iStep])
    @constraint(M, xu_3[iStep=1:NSteps_sim+1], w_xu[iStep] >= x[iStep]+u[iStep]-1)

    @constraint(M, yu_1[iStep=1:NSteps_sim+1], w_yu[iStep] <= y[iStep])
    @constraint(M, yu_2[iStep=1:NSteps_sim+1], w_yu[iStep] <= u[iStep])
    @constraint(M, yu_3[iStep=1:NSteps_sim+1], w_yu[iStep] >= y[iStep]+u[iStep]-1)

    @constraint(M, zu_1[iStep=1:NSteps_sim+1], w_zu[iStep] <= z[iStep])
    @constraint(M, zu_2[iStep=1:NSteps_sim+1], w_zu[iStep] <= u[iStep])
    @constraint(M, zu_3[iStep=1:NSteps_sim+1], w_zu[iStep] >= z[iStep]+u[iStep]-1)=#

    # CONSTRAINTS ON DEGRADATION
    @constraint(M, deg_1[iStep=1:NSteps_sim], deg[iStep] >= soc_quad[iStep]/max_SOC^2 - soc_quad[iStep+1]/max_SOC^2 + (2/max_SOC)*(soc[iStep+1]-soc[iStep]))
    @constraint(M, deg_2[iStep=1:NSteps_sim], deg[iStep] >= soc_quad[iStep+1]/max_SOC^2 - soc_quad[iStep]/max_SOC^2 + (2/max_SOC)*(soc[iStep]-soc[iStep+1]))

    #CONSTRAINT ON REVAMPING
    @constraint(M, energy_capacity[iStage=1:NStages], capacity[vector_stages_index[iScen, iStage]+2] == capacity[vector_stages_index[iScen,iStage]+1]-deg[vector_stages_index[iScen, iStage]+1]*k+revamping[iStage])
   
    @constraint(M, initial_e[iStep=1], capacity[iStep] == min_SOH) #max_SOH

    #@constraint(M, en_cap[iStep=1:NSteps], capacity[iStep+1]== capacity[iStep]-deg[iStep]*k)

    @constraint(M,en_cap[iStage in 1:NStages, iStep in (vector_stages_index[iScen,iStage]+2:vector_stages_index[iScen,iStage+1])], capacity[iStep+1]== capacity[iStep]-deg[iStep]*k) 
    

    return BuildStageProblem_4(
        M,
        soc,
        soc_quad,
        charge,
        discharge,
        deg,
        x,
        y,
        z,
        #u,
        w_xx,
        w_yy,
        w_zz,
        w_xy,
        w_xz,
        w_zy,
        #w_uu,
        #w_xu,
        #w_yu,
        #w_zu,
        capacity,
        revamping,
        e,
      )
end


function BuildStageProblem_5(InputParameters::InputParam, SolverParameters::SolverParam, Battery::BatteryParam, NSteps_sim, vector_stages_index, iScen, vector_prices)       

    @unpack (MIPGap, MIPFocus, Method, Cuts, Heuristics) = SolverParameters;
  
    @unpack (NYears, NMonths, NStages, Big, NHoursStep, conv, bin) = InputParameters;     #NSteps,NHoursStage
    @unpack (min_SOC, max_SOC, min_P, max_P, Eff_charge, Eff_discharge, max_SOH, min_SOH, Nfull) = Battery ;         
    disc = 31

    k = min_SOH/(2*Nfull)
    Small = 0.64

    M = Model(Gurobi.Optimizer)
    set_optimizer_attribute(M, "MIPGap", 0.01)

    # DEFINE VARIABLES

    @variable(M, min_SOC <= soc[iStep=1:NSteps_sim+1] <= max_SOC, base_name = "Energy")                # MWh   energy_Capacity NSteps
    @variable(M, min_SOC^2 <= soc_quad[iStep=1:NSteps_sim+1] <= max_SOC^2, base_name = "Square energy")

    @variable(M, min_P <= charge[iStep=1:NSteps_sim] <= max_P, base_name= "Charge")      #max_disc   0<=discharge<=1
    @variable(M, min_P <= discharge[iStep=1:NSteps_sim] <= max_P, base_name= "Discharge")
    
    @variable(M, 0 <= deg[iStep=1:NSteps_sim] <= Small, base_name = "Degradation")

    @variable(M, 0 <= revamping[iStage=1:NStages] <= (max_SOH-min_SOH), base_name = "Revamping")
    @variable(M, min_SOH <= capacity[iStep=1:NSteps_sim+1] <= max_SOH, base_name = "Energy_Capacity")        #energy_Capacity     [iStage=1:NStages]
    
    @variable(M, e[iStep=1:NSteps_sim], Bin, base_name ="Binary operation")

    #VARIABLES FOR DISCRETIZATION of Stored Energy

    @variable(M, x[iStep=1:NSteps_sim+1], Bin, base_name = "Binary_1")
    @variable(M, y[iStep=1:NSteps_sim+1], Bin, base_name = "Binary_2")
    @variable(M, z[iStep=1:NSteps_sim+1], Bin, base_name = "Binary_3")
    @variable(M, u[iStep=1:NSteps_sim+1], Bin, base_name = "Binary_4")
    @variable(M, t[iStep=1:NSteps_sim+1], Bin, base_name = "Binary_5")
    
    @variable(M, 0<= w_xx[iStep=1:NSteps_sim+1] <= 1, base_name = "xx")
    @variable(M, 0<= w_yy[iStep=1:NSteps_sim+1] <= 1, base_name = "yy")
    @variable(M, 0<= w_zz[iStep=1:NSteps_sim+1] <= 1, base_name = "zz")
    @variable(M, 0<= w_xy[iStep=1:NSteps_sim+1] <= 1, base_name = "xy")
    @variable(M, 0<= w_xz[iStep=1:NSteps_sim+1] <= 1, base_name = "xz")
    @variable(M, 0<= w_zy[iStep=1:NSteps_sim+1] <= 1, base_name = "yz")

    @variable(M, 0 <= w_uu[iStep=1:NSteps_sim+1] <=1, base_name = "uu")
    @variable(M, 0 <= w_xu[iStep=1:NSteps_sim+1] <=1, base_name = "xu")
    @variable(M, 0 <= w_yu[iStep=1:NSteps_sim+1] <=1, base_name = "yu")
    @variable(M, 0 <= w_zu[iStep=1:NSteps_sim+1] <=1, base_name = "zu")

    @variable(M, 0 <= w_tt[iStep=1:NSteps_sim+1] <=1, base_name = "tt")
    @variable(M, 0 <= w_tx[iStep=1:NSteps_sim+1] <=1, base_name = "tx")
    @variable(M, 0 <= w_ty[iStep=1:NSteps_sim+1] <=1, base_name = "ty")
    @variable(M, 0 <= w_tz[iStep=1:NSteps_sim+1] <=1, base_name = "tz")
    @variable(M, 0 <= w_tu[iStep=1:NSteps_sim+1] <=1, base_name = "tu")

    # DEFINE OBJECTIVE function - length(Battery_price) = NStages+1=21

    @objective(
      M,
      MathOptInterface.MAX_SENSE, 
      sum(vector_prices[iScen][iStep]*discharge[iStep]-vector_prices[iScen][iStep]*charge[iStep] for iStep=1:NSteps_sim) 
      -sum(Battery_price_purchase[iStage]*(revamping[iStage]) for iStage=1:NStages)  
      +Battery_price_purchase[NStages+1]*(capacity[end]-min_SOH)+2300   
      )
         
    # DEFINE CONSTRAINTS

    @constraint(M, Charge_op[iStep=1:NSteps_sim], charge[iStep] <= max_P*e[iStep])
    @constraint(M, Disch_op[iStep=1:NSteps_sim], discharge[iStep] <= max_P*(1-e[iStep]))

    @constraint(M,energy[iStep=1:NSteps_sim], soc[iStep] + (charge[iStep]*Eff_charge-discharge[iStep]/Eff_discharge)*NHoursStep == soc[iStep+1] )

    @constraint(M, en_bal[iStep=1:NSteps_sim+1], min_SOC + ((max_SOC-min_SOC)/disc)*(x[iStep]+2*y[iStep]+4*z[iStep]+8*u[iStep]+16*t[iStep]) == soc[iStep])
    @constraint(M, en_square[iStep=1:NSteps_sim+1], soc_quad[iStep] == min_SOC^2+ 2*min_SOC*((max_SOC-min_SOC)/disc)*(x[iStep]+2*y[iStep]+4*z[iStep]+8*u[iStep]+16*t[iStep])+(w_xx[iStep]+4*w_xy[iStep]+8*w_xz[iStep]+16*w_xu[iStep]+32*w_tx[iStep]+4*w_yy[iStep]+16*w_zy[iStep]+32*w_yu[iStep]+64*w_ty[iStep]+16*w_zz[iStep]+64*w_zu[iStep]+128*w_tz[iStep]+64*w_uu[iStep]+256*w_tu[iStep]+256*w_tt[iStep])*((max_SOC-min_SOC)/disc)^2)

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

    @constraint(M, uu_1[iStep=1:NSteps_sim+1], w_uu[iStep] <= u[iStep])
    @constraint(M, uu_2[iStep=1:NSteps_sim+1], w_uu[iStep] >= 2*u[iStep]-1)

    @constraint(M, xu_1[iStep=1:NSteps_sim+1], w_xu[iStep] <= x[iStep])
    @constraint(M, xu_2[iStep=1:NSteps_sim+1], w_xu[iStep] <= u[iStep])
    @constraint(M, xu_3[iStep=1:NSteps_sim+1], w_xu[iStep] >= x[iStep]+u[iStep]-1)

    @constraint(M, yu_1[iStep=1:NSteps_sim+1], w_yu[iStep] <= y[iStep])
    @constraint(M, yu_2[iStep=1:NSteps_sim+1], w_yu[iStep] <= u[iStep])
    @constraint(M, yu_3[iStep=1:NSteps_sim+1], w_yu[iStep] >= y[iStep]+u[iStep]-1)

    @constraint(M, zu_1[iStep=1:NSteps_sim+1], w_zu[iStep] <= z[iStep])
    @constraint(M, zu_2[iStep=1:NSteps_sim+1], w_zu[iStep] <= u[iStep])
    @constraint(M, zu_3[iStep=1:NSteps_sim+1], w_zu[iStep] >= z[iStep]+u[iStep]-1)

    @constraint(M, tt_1[iStep=1:NSteps_sim+1], w_tt[iStep] <= t[iStep])
    @constraint(M, tt_2[iStep=1:NSteps_sim+1], w_tt[iStep] >= 2*t[iStep]-1)

    @constraint(M, tx_1[iStep=1:NSteps_sim+1], w_tx[iStep] <= x[iStep])
    @constraint(M, tx_2[iStep=1:NSteps_sim+1], w_tx[iStep] <= t[iStep])
    @constraint(M, tx_3[iStep=1:NSteps_sim+1], w_tx[iStep] >= x[iStep]+t[iStep]-1)

    @constraint(M, ty_1[iStep=1:NSteps_sim+1], w_ty[iStep] <= y[iStep])
    @constraint(M, ty_2[iStep=1:NSteps_sim+1], w_ty[iStep] <= t[iStep])
    @constraint(M, ty_3[iStep=1:NSteps_sim+1], w_ty[iStep] >= y[iStep]+t[iStep]-1)

    @constraint(M, tz_1[iStep=1:NSteps_sim+1], w_tz[iStep] <= z[iStep])
    @constraint(M, tz_2[iStep=1:NSteps_sim+1], w_tz[iStep] <= t[iStep])
    @constraint(M, tz_3[iStep=1:NSteps_sim+1], w_tz[iStep] >= z[iStep]+t[iStep]-1)

    @constraint(M, tu_1[iStep=1:NSteps_sim+1], w_tu[iStep] <= t[iStep])
    @constraint(M, tu_2[iStep=1:NSteps_sim+1], w_tu[iStep] <= u[iStep])
    @constraint(M, tu_3[iStep=1:NSteps_sim+1], w_tu[iStep] >= t[iStep]+u[iStep]-1) 
    

    # CONSTRAINTS ON DEGRADATION
    @constraint(M, deg_1[iStep=1:NSteps_sim], deg[iStep] >= soc_quad[iStep]/max_SOC^2 - soc_quad[iStep+1]/max_SOC^2 + (2/max_SOC)*(soc[iStep+1]-soc[iStep]))
    @constraint(M, deg_2[iStep=1:NSteps_sim], deg[iStep] >= soc_quad[iStep+1]/max_SOC^2 - soc_quad[iStep]/max_SOC^2 + (2/max_SOC)*(soc[iStep]-soc[iStep+1]))

    #CONSTRAINT ON REVAMPING
    @constraint(M, energy_capacity[iStage=1:NStages], capacity[vector_stages_index[iScen, iStage]+2] == capacity[vector_stages_index[iScen,iStage]+1]-deg[vector_stages_index[iScen, iStage]+1]*k+revamping[iStage])
   
    @constraint(M, initial_e[iStep=1], capacity[iStep] == min_SOH) #max_SOH

    #@constraint(M, en_cap[iStep=1:NSteps], capacity[iStep+1]== capacity[iStep]-deg[iStep]*k)

    @constraint(M,en_cap[iStage in 1:NStages, iStep in (vector_stages_index[iScen,iStage]+2:vector_stages_index[iScen,iStage+1])], capacity[iStep+1]== capacity[iStep]-deg[iStep]*k) 
    

    return BuildStageProblem_5(
        M,
        soc,
        soc_quad,
        charge,
        discharge,
        deg,
        x,
        y,
        z,
        u,
        w_xx,
        w_yy,
        w_zz,
        w_xy,
        w_xz,
        w_zy,
        w_uu,
        w_xu,
        w_yu,
        w_zu,
        capacity,
        revamping,
        e,
        t,
        w_tt,
        w_tx,
        w_ty,
        w_tz,
        w_tu,
      )
end


