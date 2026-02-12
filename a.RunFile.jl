#RUN FILE

# Calls the Packages used for the optimization problem
using JuMP
using Printf
using Gurobi
using MathOptInterface
using JLD
using TimerOutputs
using DataFrames
using XLSX
using Parameters
using Dates
using CSV
using Statistics

# Calls the other Julia files
include("Structures.jl")
include("SetInputParameters.jl") 

include("Saving in xlsx.jl")

include("formulation_without_deg.jl")
include("problem_without_deg.jl")

include("solveOptimizationAlgorithm.jl")      
include("ProblemFormulationInequalities.jl") 

include("solve_OCRA_2.jl")
include("OCRA_2_Stage_problem.jl")

include("Ex_post analysis.jl")

date = string(today())

# PREPARE INPUT DATA
to = TimerOutput()

@timeit to "Set input data" begin

  #Set run case - indirizzi delle cartelle di input ed output
  case = set_runCase()

  @unpack (DataPath,InputPath,ResultPath,CaseName, Results_JLD) = case;

  # Set run mode (how and what to run) and Input parameters
  runMode = read_runMode_file()
  @unpack (OCRA_1) = runMode

  InputParameters = set_parameters(runMode, case)
  @unpack (NYears, NMonths, NHoursStep, NStages, Big, conv, bin)= InputParameters;    #NSteps, NHoursStage

  # Upload battery's characteristics
  Battery = set_battery_system(runMode, case)
  @unpack (min_SOC, max_SOC, Eff_charge, Eff_discharge, min_P, max_P, max_SOH, min_SOH, Nfull, downtime, fix, cost) = Battery; 

  # Set solver parameters (Gurobi etc)
  SolverParameters = set_solverParameters()

  # Set BESS revamping costs 
  Battery_price_purchase = read_csv("Yearly_battery_mid_cost.csv",case.DataPath) 
  Battery_price_sale = set_price(Battery_price_purchase,cost);
 
  # Set scenarios to simulate
  Pp = read_csv("100_scenarios_non_negative.csv", case.DataPath); #100_scenarios_non_negative
  NScen = size(Pp)[2]
  NSteps = size(Pp)[1]
  Steps_stages = [0 8760 17520 26280 35040 43800 52560 61320 70080 78840 87600]
  #[0 4380 8760 13140 17520 21900 26280 30660 35040 39420 43800 48180 52560 56940 61320 65700 70080 74460 78840 83220 87600]
 #[0 4344 8760 13104 17520 21864 26280 30624 35040 39384 43800 48144 52560 56904 61320 65664 70080 74424 78840 83184 87600]
 
end

@timeit to "Solve problem WITHOUT degradation" begin
  if bin ==3
    Results_No_Deg_3 = solveWithoutDeg_3(InputParameters, SolverParameters, Battery, NScen, NSteps, Pp, Steps_stages)
    #save(joinpath(FinalResPath, "optimization_results_no_deg.jld"), "optimization_results_no_deg", Results_No_Deg_3)
  elseif bin == 4
    Results_No_Deg_4 = solveWithoutDeg_4(InputParameters, SolverParameters, Battery, NScen, NSteps, Pp, Steps_stages)
    #save(joinpath(FinalResPath, "optimization_results_no_deg.jld"), "optimization_results_no_deg", Results_No_Deg_4)
  elseif bin == 5
    Results_No_Deg_5 = solveWithoutDeg_5(InputParameters, SolverParameters, Battery, NScen, NSteps, Pp, Steps_stages)
    #save(joinpath(FinalResPath, "optimization_results_no_deg.jld"), "optimization_results_no_deg", Results_No_Deg_5)
  else
    println("Inserire un numero tra 3, 4 e 5")
    throw(error())
  end
end

@timeit to "Ex-post costs evaluation and price filtering" begin
  if bin == 3
    Results_ex_post = ex_post_analysis_3(Results_No_Deg_3, InputParameters, Battery, NScen, NSteps, Pp, Steps_stages)
  elseif bin == 4
    Results_ex_post = ex_post_analysis_4(Results_No_Deg_4, InputParameters, Battery, NScen, NSteps, Pp, Steps_stages)
  else
    Results_ex_post = ex_post_analysis_5(Results_No_Deg_5, InputParameters, Battery, NScen, NSteps, Pp, Steps_stages)
  end
end

PathJLDFiles = set_run_name(run_mode, Results_JLD, Configuration)

main=0
@timeit to "Save results WITHOUT degradation" begin
  #cartella = "C:\\GitHub_Asja\\Ocra_1.0_multi_scenario\\Results_multi_scenario"
  cd(ResultsPath)
  if bin == 3
    main = data_saving_without_deg_3(Results_No_Deg_3, Results_ex_post, NSteps, NStages)
    save(joinpath(PathJLDFiles, "Results_No_deg_3.jld"), "Results_No_Deg_3", Results_No_Deg_3)
    save(joinpath(PathJLDFiles, "Results_ex_post.jld"), "Results_ex_post", Results_ex_post)
  elseif bin == 4
    main = data_saving_without_deg_4(Results_No_Deg_4, Results_ex_post, NSteps, NStages)
    save(joinpath(PathJLDFiles, "Results_No_Deg_4.jld"), "Results_No_Deg_4", Results_No_Deg_4)
    save(joinpath(PathJLDFiles, "Results_ex_post.jld"), "Results_ex_post", Results_ex_post)
  else
    main = data_saving_without_deg_5(Results_No_Deg_5, Results_ex_post, NSteps, NStages)
    save(joinpath(PathJLDFiles, "Results_No_Deg_5.jld"), "Results_No_Deg_5", Results_No_Deg_5)
    save(joinpath(PathJLDFiles, "Results_ex_post.jld"), "Results_ex_post", Results_ex_post)
  end
end


@timeit to "Solve optimization problem WITH degradation" begin
  if OCRA_1 #if this is true, we solve the OCRA 1.0
    if bin == 3
      ResultsOpt_3 = solveOptimizationProblem_3(InputParameters, SolverParameters, Battery, NScen, Results_ex_post);
    elseif bin == 4
      ResultsOpt_4 = solveOptimizationProblem_4(InputParameters, SolverParameters, Battery, NScen, Results_ex_post);
    else
      ResultsOpt_5 = solveOptimizationProblem_5(InputParameters, SolverParameters, Battery, NScen, Results_ex_post);
    end
  else  #OTHERWISE WE SOLVE THE OCRA 2 BUT ONLY WITH THREE VARIABLES
    if bin !=3
      error("You have to re-run the normal problem with only 3 binary variables!")
    else
      Results_OCRA_2 = solveOCRA_2(InputParameters, SolverParameters, Battery, NScen, Results_ex_post)
    end
  end
end

@timeit to "Evaluate results WITH degradation: statistical analysis" begin
  if OCRA_1
    if bin == 3
      Results_statistics = statistical_analysis_3(ResultsOpt_3, Results_ex_post, NScen, NStages)
    elseif bin == 4
      Results_statistics = statistical_analysis_4(ResultsOpt_4, Results_ex_post, NScen, NStages)
    else
      Results_statistics = statistical_analysis_5(ResultsOpt_5, Results_ex_post, NScen, NStages)
    end
  else
      Results_statistics_OCRA_2 = analysis_OCRA_2(Results_OCRA_2, Results_ex_post, NScen, NStages)
  end
end

# SAVE DATA IN EXCEL FILES
if runMode.excel_savings
  cd(main)
  if OCRA_1
    if bin == 3
      Saving = data_saving_3(InputParameters, ResultsOpt_3, Results_ex_post, Results_statistics, Battery, NScen)
    elseif bin == 4
      Saving = data_saving_4(InputParameters, ResultsOpt_4, Results_ex_post, Results_statistics, Battery, NScen)
    else
      Saving = data_saving_5(InputParameters, ResultsOpt_5, Results_ex_post, Results_statistics, Battery, NScen)
    end
  else
      Saving_OCRA_2 = saving_OCRA_2(InputParameters, Results_OCRA_2, Results_ex_post, Results_statistics_OCRA_2, Battery, NScen)
  end
else
  println("Solved without saving results in xlsx format.")
end


#end
print(to)




