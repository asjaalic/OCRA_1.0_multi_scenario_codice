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

include("solveOptimizationAlgorithm.jl")        #solveOptimizationAlgorithm_3cuts
include("ProblemFormulationInequalities.jl")      #ProblemFormulationCutsTaylor_3

include("Saving in xlsx.jl")

include("formulation_without_deg.jl")
include("problem_without_deg.jl")
include("Ex_post analysis.jl")

date = string(today())

# PREPARE INPUT DATA
to = TimerOutput()

@timeit to "Set input data" begin

  #Set run case - indirizzi delle cartelle di input ed output
  case = set_runCase()

  @unpack (DataPath,InputPath,ResultPath,CaseName) = case;

  # Set run mode (how and what to run) and Input parameters
  runMode = read_runMode_file()
  InputParameters = set_parameters(runMode, case)
  @unpack (NYears, NMonths, NHoursStep, NStages, Big, conv, disc)= InputParameters;    #NSteps, NHoursStage

  # Upload battery's characteristics
  Battery = set_battery_system(runMode, case)
  @unpack (min_SOC, max_SOC, Eff_charge, Eff_discharge, min_P, max_P, max_SOH, min_SOH, Nfull) = Battery; 

  # Set solver parameters (Gurobi etc)
  SolverParameters = set_solverParameters()

  # Set BESS revamping costs 
  Battery_price_purchase = read_csv("Battery_decreasing_prices_high.csv",case.DataPath) 
 
  # Set scenarios to simulate
  Pp = read_csv("5_scenari.csv", case.DataPath);
  NScen = size(Pp)[2]
  NSteps = size(Pp)[1]
  Steps_stages = [0 4380 8760 13140 17520 21900 26280 30660 35040 39420 43800 48180 52560 56940 61320 65700 70080 74460 78840 83220 87600]
 #[0 4344 8760 13104 17520 21864 26280 30624 35040 39384 43800 48144 52560 56904 61320 65664 70080 74424 78840 83184 87600]
 
  # Where and how to save the results
  FinalResPath= set_run_name(case, ResultPath, NSteps)

end

@timeit to "Solve problem WITHOUT degradation" begin
  Results_No_Deg = solveWithoutDeg(InputParameters, SolverParameters, Battery, NScen, NSteps, Pp, Steps_stages)
  save(joinpath(FinalResPath, "optimization_results_no_deg.jld"), "optimization_results_no_deg", Results_No_Deg)
end

@timeit to "Ex-post costs evaluation and price filtering" begin
  Results_ex_post = ex_post_analysis(Results_No_Deg, InputParameters, Battery, NScen, NSteps, Pp, Steps_stages)
  save(joinpath(FinalResPath, "results_ex_post.jld"), "results_ex_post", Results_ex_post)
end

main=0
@timeit to "Save results WITHOUT degradation" begin
  cartella = "C:\\Users\\Utente\\Desktop\\ASJA\\OCRA_1.0_UNI_PAVIA\\Risultati_completo"
  cd(cartella)
  main = data_saving_without_deg(Results_No_Deg, Results_ex_post, NSteps, NStages)
end


@timeit to "Solve optimization problem WITH degradation" begin
  ResultsOpt = solveOptimizationProblem(InputParameters, SolverParameters, Battery, NScen, Results_ex_post);
  save(joinpath(FinalResPath, "optimization_results.jld"), "optimization_results", ResultsOpt)
end

@timeit to "Evaluate results WITH degradation: statistical analysis" begin
  Results_statistics = statistical_analysis(ResultsOpt, Results_ex_post, NScen, NStages)
end

# SAVE DATA IN EXCEL FILES
if runMode.excel_savings
  cd(main)
  Saving = data_saving(InputParameters, ResultsOpt, Results_ex_post, Results_statistics, Battery, NScen)
else
  println("Solved without saving results in xlsx format.")
end


#end
print(to)




