# STRUCTURES USED IN THE PROBLEM

# Input data
#-----------------------------------------------

# Input parameters 
@with_kw struct InputParam{F<:Float64,I<:Int}
    NYears::F                                     # Number of years
    NMonths::I
    NStages::I                                    # Number of stages of N months in the problem FORMULATION-- calcolato come NYears/NMonths*12
    NHoursStep::F                                 # Number of hours in each time step
    Big::F                                        # A big number
    conv::F                                       # A small number for degradation convergence
    disc::I                                       # Discretization points
end

# Battery's characteristics
@with_kw struct BatteryParam{F<:Float64,I<:Int}
    min_SOC::F                                     # Battery's minimum energy storage capacity
    max_SOC::F                                     # Batter's maximum energy storage capacity
    Eff_charge::F                                  # Battery's efficiency for charging operation
    Eff_discharge::F                               # Battery's efficiency for discharging operation
    min_P::F
    max_P::F
    max_SOH::F                                     # Maximum SOH that can be achieved because of volume issues
    min_SOH::F                                     # Minimum SOH to be respected by contract
    Nfull::I                                       # Maximum number of full cycles for DoD=100%     
    downtime::I
end
  
# solver parameters
@with_kw struct SolverParam{F<:Float64,I<:Int}
    MIPGap::F 
    MIPFocus::I
    Method::F
    Cuts::F
    Heuristics::F
end
  
# Indirizzi cartelle
@with_kw struct caseData{S<:String}
    DataPath::S
    InputPath::S
    ResultPath::S
    CaseName::S
end

# runMode Parameters
@with_kw mutable struct runModeParam{B<:Bool}

    # Solver settings
    solveMIP::B     #If using SOS2

    batterySystemFromFile::B 

    #runMode self defined reading of input 
    setInputParameters::B             #from .in file
 
    excel_savings::B 

end

# Optimization problem
struct BuildStageProblem
    M::Any
    soc::Any
    soc_quad::Any
    charge::Any 
    discharge::Any
    deg::Any
    x::Any
    y::Any
    z::Any
    #u::Any
    w_xx::Any
    w_yy::Any
    w_zz::Any
    w_xy::Any
    w_xz::Any
    w_zy::Any
    #w_uu::Any
    #w_xu::Any
    #w_yu::Any
    #w_zu::Any
    capacity::Any
    revamping::Any
    e::Any
    #=t::Any
    w_tt::Any
    w_tx::Any
    w_ty::Any
    w_tz::Any
    w_tu::Any=#
end

struct Results
    objective::Any
    net_revenues_per_stage::Any
    WEM_stage::Any
    cost_rev::Any
    deg_stage::Any
    soc::Any
    charge::Any
    discharge::Any
    deg::Any
    soc_quad::Any
    x::Any
    y::Any
    z::Any
    #u::Any
    w_xx::Any
    w_yy::Any
    w_zz::Any
   # w_uu::Any
    w_xy::Any
    w_xz::Any
    w_zy::Any
    #w_xu::Any
    #w_yu::Any
    #w_zu::Any
    rev::Any
    cap::Any
end

# Optimization problem
struct BuildStageProblem_no_deg
    M_sim::Any
    soc::Any
    soc_quad::Any
    charge::Any 
    discharge::Any
    x::Any
    y::Any
    z::Any
    #u::Any
    w_xx::Any
    w_yy::Any
    w_zz::Any
    w_xy::Any
    w_xz::Any
    w_zy::Any
    #w_uu::Any
    #w_xu::Any
    #w_yu::Any
    #w_zu::Any
    e::Any
    #=t::Any
    w_tt::Any
    w_tx::Any
    w_ty::Any
    w_tz::Any
    w_tu::Any=#
end

struct ResultWithoutDeg
    objective_no_deg::Any
    revenues_per_stage_no_deg::Any
    soc_no_deg::Any
    charge_no_deg::Any
    discharge_no_deg::Any
    soc_quad_no_deg::Any
    x_no_deg::Any
    y_no_deg::Any
    z_no_deg::Any
    #u_no_deg::Any
    w_xx_no_deg::Any
    w_yy_no_deg::Any
    w_zz_no_deg::Any
    #w_uu_no_deg::Any
    w_xy_no_deg::Any
    w_xz_no_deg::Any
    w_zy_no_deg::Any
    #w_xu_no_deg::Any
    #w_yu_no_deg::Any
    #w_zu_no_deg::Any
    e_no_deg::Any
end

struct results_ex_post
    deg_singola::Any
    deg_stage::Any
    costo_stage::Any
    net_revenues::Any
    vector_prices::Any
    vector_stages_index::Any
    NSteps_scenario::Any
    vector_downtime_stages::Any
end

struct Data_analysed
    WEM_rev_stat_frame::Any
    cost_stat_frame::Any
    net_rev_stat_frame::Any
    deg_stat_frame::Any
    Initial_cap_stat_frame::Any
    Final_cap_stat_frame::Any
end