# STAGE MAXIMIZATION PROBLEM FORMULATION

function BuildStageProblemNoDeg(InputParameters::InputParam, SolverParameters::SolverParam, Battery::BatteryParam)       #, state_variables::states When we have 2 hydropower plants- 2 turbines

    @unpack (MIPGap, MIPFocus, Method, Cuts, Heuristics) = SolverParameters;
  
    @unpack (NYears, NMonths, NStages, Big, NHoursStep, conv, disc) = InputParameters;     #NSteps,NHoursStage
    @unpack (min_SOC, max_SOC, min_P, max_P, Eff_charge, Eff_discharge, max_SOH, min_SOH, Nfull) = Battery ;         

    M_sim = Model(Gurobi.Optimizer)
    set_optimizer_attribute(M_sim, "MIPGap", 0.01)

    # DEFINE VARIABLES
    @variable(M_sim, min_SOC <= soc[iStep=1:NSteps+1] <= max_SOC, base_name = "Energy")                # MWh   energy_Capacity NSteps
    @variable(M_sim, min_SOC^2 <= soc_quad[iStep=1:NSteps+1] <= max_SOC^2, base_name = "Square energy")

    @variable(M_sim, min_P <= charge[iStep=1:NSteps] <= max_P, base_name= "Charge")      #max_disc   0<=discharge<=1
    @variable(M_sim, min_P <= discharge[iStep=1:NSteps] <= max_P, base_name= "Discharge")

    @variable(M_sim, e[iStep=1:NSteps], Bin, base_name ="Binary operation")

    #VARIABLES FOR DISCRETIZATION of Stored Energy

    @variable(M_sim, x[iStep=1:NSteps+1], Bin, base_name = "Binary_1")
    @variable(M_sim, y[iStep=1:NSteps+1], Bin, base_name = "Binary_2")
    @variable(M_sim, z[iStep=1:NSteps+1], Bin, base_name = "Binary_3")
    #@variable(M_sim, u[iStep=1:NSteps+1], Bin, base_name = "Binary_4")
    #@variable(M_sim, t[iStep=1:NSteps+1], Bin, base_name = "Binary_5")
    
    @variable(M_sim, 0<= w_xx[iStep=1:NSteps+1] <= 1, base_name = "xx")
    @variable(M_sim, 0<= w_yy[iStep=1:NSteps+1] <= 1, base_name = "yy")
    @variable(M_sim, 0<= w_zz[iStep=1:NSteps+1] <= 1, base_name = "zz")
    @variable(M_sim, 0<= w_xy[iStep=1:NSteps+1] <= 1, base_name = "xy")
    @variable(M_sim, 0<= w_xz[iStep=1:NSteps+1] <= 1, base_name = "xz")
    @variable(M_sim, 0<= w_zy[iStep=1:NSteps+1] <= 1, base_name = "yz")

  #=  @variable(M_sim, 0 <= w_uu[iStep=1:NSteps+1] <=1, base_name = "uu")
    @variable(M_sim, 0 <= w_xu[iStep=1:NSteps+1] <=1, base_name = "xu")
    @variable(M_sim, 0 <= w_yu[iStep=1:NSteps+1] <=1, base_name = "yu")
    @variable(M_sim, 0 <= w_zu[iStep=1:NSteps+1] <=1, base_name = "zu")=#

  #=  @variable(M_sim, 0 <= w_tt[iStep=1:NSteps+1] <=1, base_name = "tt")
    @variable(M_sim, 0 <= w_tx[iStep=1:NSteps+1] <=1, base_name = "tx")
    @variable(M_sim, 0 <= w_ty[iStep=1:NSteps+1] <=1, base_name = "ty")
    @variable(M_sim, 0 <= w_tz[iStep=1:NSteps+1] <=1, base_name = "tz")
    @variable(M_sim, 0 <= w_tu[iStep=1:NSteps+1] <=1, base_name = "tu")=#
  
    # DEFINE OBJECTIVE function - length(Battery_price) = NStages+1=21

    @objective(
      M_sim,
      MathOptInterface.MAX_SENSE, 
      sum(1*discharge[iStep]+1*charge[iStep] for iStep=1:NSteps) )
         
    # DEFINE CONSTRAINTS

    @constraint(M_sim, Charge_op[iStep=1:NSteps], charge[iStep] <= max_P*e[iStep])
    @constraint(M_sim, Disch_op[iStep=1:NSteps], discharge[iStep] <= max_P*(1-e[iStep]))

    @constraint(M_sim, energy[iStep=1:NSteps], soc[iStep] + (charge[iStep]*Eff_charge-discharge[iStep]/Eff_discharge)*NHoursStep == soc[iStep+1] )

    @constraint(M_sim, en_bal[iStep=1:NSteps+1], min_SOC + ((max_SOC-min_SOC)/disc)*(x[iStep]+2*y[iStep]+4*z[iStep]) == soc[iStep])
   # @constraint(M_sim, en_bal[iStep=1:NSteps+1], min_SOC + ((max_SOC-min_SOC)/disc)*(x[iStep]+2*y[iStep]+4*z[iStep]+8*u[iStep]) == soc[iStep])    
    #@constraint(M_sim, en_bal[iStep=1:NSteps], min_SOC + ((max_SOC-min_SOC)/disc)*(x[iStep]+2*y[iStep]+4*z[iStep]+8*u[iStep]+16*t[iStep]) == soc[iStep])

    @constraint(M_sim, en_square[iStep=1:NSteps+1], soc_quad[iStep] == min_SOC^2+ 2*min_SOC*((max_SOC-min_SOC)/disc)*(x[iStep]+2*y[iStep]+4*z[iStep])+(w_xx[iStep]+4*w_xy[iStep]+8*w_xz[iStep]+4*w_yy[iStep]+16*w_zz[iStep]+16*w_zy[iStep])*((max_SOC-min_SOC)/disc)^2)
    #@constraint(M_sim, en_square[iStep=1:NSteps+1], soc_quad[iStep] == min_SOC^2+ 2*min_SOC*((max_SOC-min_SOC)/disc)*(x[iStep]+2*y[iStep]+4*z[iStep]+8*u[iStep])+(w_xx[iStep]+4*w_xy[iStep]+8*w_xz[iStep]+16*w_xu[iStep]+4*w_yy[iStep]+16*w_zz[iStep]+16*w_zy[iStep]+32*w_yu[iStep]+64*w_zu[iStep]+64*w_uu[iStep])*((max_SOC-min_SOC)/disc)^2)
    #@constraint(M_sim, en_square[iStep=1:NSteps+1], soc_quad[iStep] == min_SOC^2+ 2*min_SOC*((max_SOC-min_SOC)/disc)*(x[iStep]+2*y[iStep]+4*z[iStep]+8*u[iStep]+16*t[iStep])+(w_xx[iStep]+4*w_xy[iStep]+8*w_xz[iStep]+16*w_xu[iStep]+32*w_tx[iStep]+4*w_yy[iStep]+16*w_zy[iStep]+32*w_yu[iStep]+64*w_ty[iStep]+16*w_zz[iStep]+64*w_zu[iStep]+128*w_tz[iStep]+64*w_uu[iStep]+256*w_tu[iStep]+256*w_tt[iStep])*((max_SOC-min_SOC)/disc)^2)

    # INEQUALITIES CONSTRAINTS
    @constraint(M_sim, xx_1[iStep=1:NSteps+1], w_xx[iStep] <= x[iStep])
    @constraint(M_sim, xx_2[iStep=1:NSteps+1], w_xx[iStep] >= 2*x[iStep]-1)

    @constraint(M_sim, xy_1[iStep=1:NSteps+1], w_xy[iStep] <= x[iStep])
    @constraint(M_sim, xy_2[iStep=1:NSteps+1], w_xy[iStep] <= y[iStep])
    @constraint(M_sim, xy_3[iStep=1:NSteps+1], w_xy[iStep] >= x[iStep]+y[iStep]-1)

    @constraint(M_sim, xz_1[iStep=1:NSteps+1], w_xz[iStep] <= x[iStep])
    @constraint(M_sim, xz_2[iStep=1:NSteps+1], w_xz[iStep] <= z[iStep])
    @constraint(M_sim, xz_3[iStep=1:NSteps+1], w_xz[iStep] >= x[iStep]+z[iStep]-1)

    @constraint(M_sim, yy_1[iStep=1:NSteps+1], w_yy[iStep] <= y[iStep])
    @constraint(M_sim, yy_2[iStep=1:NSteps+1], w_yy[iStep] >= 2*y[iStep]-1)

    @constraint(M_sim, zz_1[iStep=1:NSteps+1], w_zz[iStep] <= z[iStep])
    @constraint(M_sim, zz_2[iStep=1:NSteps+1], w_zz[iStep] >= 2*z[iStep]-1)

    @constraint(M_sim, zy_1[iStep=1:NSteps+1], w_zy[iStep] <= z[iStep])
    @constraint(M_sim, zy_2[iStep=1:NSteps+1], w_zy[iStep] <= y[iStep])
    @constraint(M_sim, zy_3[iStep=1:NSteps+1], w_zy[iStep] >= z[iStep]+y[iStep]-1)

   #= @constraint(M_sim, uu_1[iStep=1:NSteps+1], w_uu[iStep] <= u[iStep])
    @constraint(M_sim, uu_2[iStep=1:NSteps+1], w_uu[iStep] >= 2*u[iStep]-1)

    @constraint(M_sim, xu_1[iStep=1:NSteps+1], w_xu[iStep] <= x[iStep])
    @constraint(M_sim, xu_2[iStep=1:NSteps+1], w_xu[iStep] <= u[iStep])
    @constraint(M_sim, xu_3[iStep=1:NSteps+1], w_xu[iStep] >= x[iStep]+u[iStep]-1)

    @constraint(M_sim, yu_1[iStep=1:NSteps+1], w_yu[iStep] <= y[iStep])
    @constraint(M_sim, yu_2[iStep=1:NSteps+1], w_yu[iStep] <= u[iStep])
    @constraint(M_sim, yu_3[iStep=1:NSteps+1], w_yu[iStep] >= y[iStep]+u[iStep]-1)

    @constraint(M_sim, zu_1[iStep=1:NSteps+1], w_zu[iStep] <= z[iStep])
    @constraint(M_sim, zu_2[iStep=1:NSteps+1], w_zu[iStep] <= u[iStep])
    @constraint(M_sim, zu_3[iStep=1:NSteps+1], w_zu[iStep] >= z[iStep]+u[iStep]-1)=#

    #= PER 5 VARIABILI binarie

    @constraint(M_sim, tt_1[iStep=1:NSteps+1], w_tt[iStep] <= t[iStep])
    @constraint(M_sim, tt_2[iStep=1:NSteps+1], w_tt[iStep] >= 2*t[iStep]-1)

    @constraint(M_sim, tx_1[iStep=1:NSteps+1], w_tx[iStep] <= x[iStep])
    @constraint(M_sim, tx_2[iStep=1:NSteps+1], w_tx[iStep] <= t[iStep])
    @constraint(M_sim, tx_3[iStep=1:NSteps+1], w_tx[iStep] >= x[iStep]+t[iStep]-1)

    @constraint(M_sim, ty_1[iStep=1:NSteps+1], w_ty[iStep] <= y[iStep])
    @constraint(M_sim, ty_2[iStep=1:NSteps+1], w_ty[iStep] <= t[iStep])
    @constraint(M_sim, ty_3[iStep=1:NSteps+1], w_ty[iStep] >= y[iStep]+t[iStep]-1)

    @constraint(M_sim, tz_1[iStep=1:NSteps+1], w_tz[iStep] <= z[iStep])
    @constraint(M_sim, tz_2[iStep=1:NSteps+1], w_tz[iStep] <= t[iStep])
    @constraint(M_sim, tz_3[iStep=1:NSteps+1], w_tz[iStep] >= z[iStep]+t[iStep]-1)

    @constraint(M_sim, tu_1[iStep=1:NSteps+1], w_tu[iStep] <= t[iStep])
    @constraint(M_sim, tu_2[iStep=1:NSteps+1], w_tu[iStep] <= u[iStep])
    @constraint(M_sim, tu_3[iStep=1:NSteps+1], w_tu[iStep] >= t[iStep]+u[iStep]-1)=#
    
    return BuildStageProblem_no_deg(
        M_sim,
        soc,
        soc_quad,
        charge,
        discharge,
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
       #= w_uu,
        w_xu,
        w_yu,
        w_zu,=#
        e,
        #=t,
        w_tt,
        w_tx,
        w_ty,
        w_tz,
        w_tu,=#
      )
end



