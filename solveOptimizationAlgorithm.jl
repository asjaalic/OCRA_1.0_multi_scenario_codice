# SOLVE OPTIMIZATION PROBLEM

function solveOptimizationProblem_3(InputParameters::InputParam, SolverParameters::SolverParam, Battery::BatteryParam, NScen, Results_ex_post::results_ex_post)

    @unpack (NYears, NMonths, NStages, Big, NHoursStep, bin) = InputParameters;                #NSteps, NHoursStage
    @unpack (min_SOC, max_SOC, Eff_charge, Eff_discharge, min_P, max_P, max_SOH, min_SOH, Nfull ) = Battery;
    @unpack (vector_prices, vector_stages_index, NSteps_scenario) = Results_ex_post;

    println("Solving Optimization Problem")

    k = min_SOH/(2*Nfull)

    # PER OGNI SCENARIO, C'E' SEMPRE LO STESSO NUMERO DI STAGES(semestri)
    objective = zeros(NScen);                 
    net_revenues_per_stage = zeros(NScen, NStages);
    WEM_stage = zeros(NScen, NStages);
    cost_rev = zeros(NScen, NStages);
    deg_stage = zeros(NScen, NStages);

    # PER OGNI SCENARIO, CAMBIA IL NUMERO DI STEPS
    charge = Vector{Vector{Float64}}(undef,NScen)
    discharge = Vector{Vector{Float64}}(undef,NScen)
    soc = Vector{Vector{Float64}}(undef,NScen)
    deg = Vector{Vector{Float64}}(undef,NScen)
    cap = Vector{Vector{Float64}}(undef,NScen)
    rev= zeros(NScen, NStages);

    soc_quad = Vector{Vector{Float64}}(undef,NScen)
    x = Vector{Vector{Float64}}(undef,NScen)
    y = Vector{Vector{Float64}}(undef,NScen)
    z = Vector{Vector{Float64}}(undef,NScen)

    w_xx = Vector{Vector{Float64}}(undef,NScen)
    w_yy = Vector{Vector{Float64}}(undef,NScen)
    w_zz = Vector{Vector{Float64}}(undef,NScen)
    w_xy = Vector{Vector{Float64}}(undef,NScen)
    w_xz = Vector{Vector{Float64}}(undef,NScen)
    w_zy = Vector{Vector{Float64}}(undef,NScen)

    problem = 0 #BuildStageProblem(InputParameters, SolverParameters, Battery)
    NSteps_sim= 0

    for iScen=1:NScen

        println("Solving scenario $iScen")
        NSteps_sim = Int(NSteps_scenario[iScen])
        problem = BuildStageProblem_3(InputParameters, SolverParameters, Battery, NSteps_sim, vector_stages_index, iScen, vector_prices)

        charge[iScen] = Float64[]
        discharge[iScen] = Float64[]
        soc[iScen] = Float64[]
        deg[iScen] = Float64[]
        cap[iScen] = Float64[]
        soc_quad[iScen] = Float64[]
        x[iScen] = Float64[]
        y[iScen] = Float64[]
        z[iScen] = Float64[]

        w_xx[iScen] = Float64[]
        w_yy[iScen] = Float64[]
        w_zz[iScen] = Float64[]
        w_xy[iScen] = Float64[]
        w_xz[iScen] = Float64[]
        w_zy[iScen] = Float64[]

        @timeit to "Solve optimization" optimize!(problem.M)

        if termination_status(problem.M) != MOI.OPTIMAL
            println("NOT OPTIMAL: ", termination_status(problem.M))
        else
            println("Optimization finished")
        end

        @timeit to "Collecting results" begin
            objective = JuMP.objective_value(problem.M)
            
            for iStep=1:NSteps_sim
                
                push!(soc[iScen],JuMP.value(problem.soc[iStep]))
                push!(soc_quad[iScen],JuMP.value(problem.soc_quad[iStep]))
                push!(charge[iScen],JuMP.value(problem.charge[iStep]))
                push!(discharge[iScen],JuMP.value(problem.discharge[iStep]))
                push!(deg[iScen],JuMP.value(problem.deg[iStep]))
                push!(cap[iScen],JuMP.value(problem.capacity[iStep]))

                push!(x[iScen],JuMP.value(problem.x[iStep]))
                push!(y[iScen],JuMP.value(problem.y[iStep]))
                push!(z[iScen],JuMP.value(problem.z[iStep]))

                push!(w_xx[iScen],JuMP.value(problem.w_xx[iStep]))
                push!(w_yy[iScen],JuMP.value(problem.w_yy[iStep]))
                push!(w_zz[iScen],JuMP.value(problem.w_zz[iStep]))

                push!(w_xy[iScen],JuMP.value(problem.w_xy[iStep]))
                push!(w_xz[iScen],JuMP.value(problem.w_xz[iStep]))
                push!(w_zy[iScen],JuMP.value(problem.w_zy[iStep]))

            end
                    
            push!(soc[iScen],JuMP.value(problem.soc[end]))
            push!(soc_quad[iScen],JuMP.value(problem.soc_quad[end]))
            push!(cap[iScen],JuMP.value(problem.capacity[end]))

            push!(x[iScen],JuMP.value(problem.x[end]))
            push!(y[iScen],JuMP.value(problem.y[end]))
            push!(z[iScen],JuMP.value(problem.z[end]))

            push!(w_xx[iScen],JuMP.value(problem.w_xx[end]))
            push!(w_yy[iScen],JuMP.value(problem.w_yy[end]))
            push!(w_zz[iScen],JuMP.value(problem.w_zz[end]))

            push!(w_xy[iScen],JuMP.value(problem.w_xy[end]))
            push!(w_xz[iScen],JuMP.value(problem.w_xz[end]))
            push!(w_zy[iScen],JuMP.value(problem.w_zy[end]))        
            
            for iStage=1:NStages
                rev[iScen,iStage] = JuMP.value(problem.revamping[iStage])
                deg_stage[iScen, iStage] = sum(deg[iScen][iStep] for iStep=(vector_stages_index[iScen, iStage]+1):(vector_stages_index[iScen,iStage+1]))*k
            end
        
            for iStage=2:(NStages-1)
                net_revenues_per_stage[iScen, iStage] = sum(vector_prices[iScen][iStep]*NHoursStep*(discharge[iScen][iStep] - charge[iScen][iStep]) for iStep=(vector_stages_index[iScen, iStage]+1):(vector_stages_index[iScen,iStage+1])) - Battery_price_purchase[iStage]*rev[iScen, iStage] 
                WEM_stage[iScen, iStage] = sum(vector_prices[iScen][iStep]*NHoursStep*(discharge[iScen][iStep] - charge[iScen][iStep]) for iStep=(vector_stages_index[iScen, iStage]+1):(vector_stages_index[iScen, iStage+1]))
                cost_rev[iScen, iStage] = Battery_price_purchase[iStage]*(rev[iScen, iStage]) 

            end

            net_revenues_per_stage[iScen, 1] = sum(vector_prices[iScen][iStep]*NHoursStep*(discharge[iScen][iStep]-charge[iScen][iStep]) for iStep=(vector_stages_index[iScen, 1]+1):(vector_stages_index[iScen, 2])) - Battery_price_purchase[1]*rev[iScen, 1]  
            WEM_stage[iScen, 1] = sum(vector_prices[iScen][iStep]*NHoursStep*(discharge[iScen][iStep]-charge[iScen][iStep]) for iStep=(vector_stages_index[iScen, 1]+1):(vector_stages_index[iScen, 2]))
            cost_rev[iScen, 1] = Battery_price_purchase[1]*rev[iScen, 1] 

            net_revenues_per_stage[iScen, NStages] = sum(vector_prices[iScen][iStep]*NHoursStep*(discharge[iScen][iStep]-charge[iScen][iStep]) for iStep=(vector_stages_index[iScen, NStages]+1):(vector_stages_index[iScen, NStages+1])) + Battery_price_purchase[NStages+1]*(cap[iScen][end]-min_SOH) - Battery_price_purchase[NStages]*rev[iScen, NStages] 
            WEM_stage[iScen, NStages]= sum(vector_prices[iScen][iStep]*NHoursStep*(discharge[iScen][iStep]-charge[iScen][iStep]) for iStep=(vector_stages_index[iScen, NStages]+1):(vector_stages_index[iScen, NStages+1]))
            cost_rev[iScen, NStages] = Battery_price_purchase[NStages]*rev[iScen, NStages] - Battery_price_purchase[NStages+1]*(cap[iScen][end]-min_SOH)  
        end
        
        problem = 0
    end
    
    println("Collected results")

    return Results_3(
        objective,
        net_revenues_per_stage,
        WEM_stage,
        cost_rev,
        deg_stage,
        soc,
        charge,
        discharge,
        deg,
        soc_quad,
        x,
        y,
        z,
        w_xx,
        w_yy,
        w_zz,
        w_xy,
        w_xz,
        w_zy,
        rev,
        cap,  
    )

end

# SOLVE OPTIMIZATION PROBLEM

function solveOptimizationProblem_4(InputParameters::InputParam, SolverParameters::SolverParam, Battery::BatteryParam, NScen, Results_ex_post::results_ex_post)

    @unpack (NYears, NMonths, NStages, Big, NHoursStep, bin) = InputParameters;                #NSteps, NHoursStage
    @unpack (min_SOC, max_SOC, Eff_charge, Eff_discharge, min_P, max_P, max_SOH, min_SOH, Nfull ) = Battery;
    @unpack (vector_prices, vector_stages_index, NSteps_scenario) = Results_ex_post;

    println("Solving Optimization Problem")

    k = min_SOH/(2*Nfull)

    # PER OGNI SCENARIO, C'E' SEMPRE LO STESSO NUMERO DI STAGES(semestri)
    objective = zeros(NScen);                 
    net_revenues_per_stage = zeros(NScen, NStages);
    WEM_stage = zeros(NScen, NStages);
    cost_rev = zeros(NScen, NStages);
    deg_stage = zeros(NScen, NStages);

    # PER OGNI SCENARIO, CAMBIA IL NUMERO DI STEPS
    charge = Vector{Vector{Float64}}(undef,NScen)
    discharge = Vector{Vector{Float64}}(undef,NScen)
    soc = Vector{Vector{Float64}}(undef,NScen)
    deg = Vector{Vector{Float64}}(undef,NScen)
    cap = Vector{Vector{Float64}}(undef,NScen)
    rev= zeros(NScen, NStages);

    soc_quad = Vector{Vector{Float64}}(undef,NScen)
    x = Vector{Vector{Float64}}(undef,NScen)
    y = Vector{Vector{Float64}}(undef,NScen)
    z = Vector{Vector{Float64}}(undef,NScen)
   # u = Vector{Vector{Float64}}(undef,NScen)
    w_xx = Vector{Vector{Float64}}(undef,NScen)
    w_yy = Vector{Vector{Float64}}(undef,NScen)
    w_zz = Vector{Vector{Float64}}(undef,NScen)
    w_xy = Vector{Vector{Float64}}(undef,NScen)
    w_xz = Vector{Vector{Float64}}(undef,NScen)
    w_zy = Vector{Vector{Float64}}(undef,NScen)
    #w_uu = Vector{Vector{Float64}}(undef,NScen)
    #w_xu = Vector{Vector{Float64}}(undef,NScen)
    #w_yu = Vector{Vector{Float64}}(undef,NScen)
    #w_zu = Vector{Vector{Float64}}(undef,NScen)

    problem = 0 #BuildStageProblem(InputParameters, SolverParameters, Battery)
    NSteps_sim= 0

    for iScen=1:NScen

        println("Solving scenario $iScen")
        NSteps_sim = Int(NSteps_scenario[iScen])
        problem = BuildStageProblem_4(InputParameters, SolverParameters, Battery, NSteps_sim, vector_stages_index, iScen, vector_prices)

        charge[iScen] = Float64[]
        discharge[iScen] = Float64[]
        soc[iScen] = Float64[]
        deg[iScen] = Float64[]
        cap[iScen] = Float64[]
        soc_quad[iScen] = Float64[]
        x[iScen] = Float64[]
        y[iScen] = Float64[]
        z[iScen] = Float64[]
        #u[iScen] = Float64[]
        w_xx[iScen] = Float64[]
        w_yy[iScen] = Float64[]
        w_zz[iScen] = Float64[]
        w_xy[iScen] = Float64[]
        w_xz[iScen] = Float64[]
        w_zy[iScen] = Float64[]
        #w_uu[iScen] = Float64[]
        #w_xu[iScen] = Float64[]
        #w_yu[iScen] = Float64[]
        #w_zu[iScen] = Float64[]

        @timeit to "Solve optimization" optimize!(problem.M)

        if termination_status(problem.M) != MOI.OPTIMAL
            println("NOT OPTIMAL: ", termination_status(problem.M))
        else
            println("Optimization finished")
        end

        @timeit to "Collecting results" begin
            objective = JuMP.objective_value(problem.M)
            
            for iStep=1:NSteps_sim
                
                push!(soc[iScen],JuMP.value(problem.soc[iStep]))
                push!(soc_quad[iScen],JuMP.value(problem.soc_quad[iStep]))
                push!(charge[iScen],JuMP.value(problem.charge[iStep]))
                push!(discharge[iScen],JuMP.value(problem.discharge[iStep]))
                push!(deg[iScen],JuMP.value(problem.deg[iStep]))
                push!(cap[iScen],JuMP.value(problem.capacity[iStep]))

                push!(x[iScen],JuMP.value(problem.x[iStep]))
                push!(y[iScen],JuMP.value(problem.y[iStep]))
                push!(z[iScen],JuMP.value(problem.z[iStep]))
                #push!(u[iScen],JuMP.value(problem.u[iStep]))

                push!(w_xx[iScen],JuMP.value(problem.w_xx[iStep]))
                push!(w_yy[iScen],JuMP.value(problem.w_yy[iStep]))
                push!(w_zz[iScen],JuMP.value(problem.w_zz[iStep]))
                #push!(w_uu[iScen],JuMP.value(problem.w_uu[iStep]))

                push!(w_xy[iScen],JuMP.value(problem.w_xy[iStep]))
                push!(w_xz[iScen],JuMP.value(problem.w_xz[iStep]))
                push!(w_zy[iScen],JuMP.value(problem.w_zy[iStep]))
                #push!(w_xu[iScen],JuMP.value(problem.w_xu[iStep]))
                #push!(w_yu[iScen],JuMP.value(problem.w_yu[iStep]))
                #push!(w_zu[iScen],JuMP.value(problem.w_zu[iStep]))

            end
                    
            push!(soc[iScen],JuMP.value(problem.soc[end]))
            push!(soc_quad[iScen],JuMP.value(problem.soc_quad[end]))
            push!(cap[iScen],JuMP.value(problem.capacity[end]))

            push!(x[iScen],JuMP.value(problem.x[end]))
            push!(y[iScen],JuMP.value(problem.y[end]))
            push!(z[iScen],JuMP.value(problem.z[end]))
            #push!(u[iScen],JuMP.value(problem.u[end]))

            push!(w_xx[iScen],JuMP.value(problem.w_xx[end]))
            push!(w_yy[iScen],JuMP.value(problem.w_yy[end]))
            push!(w_zz[iScen],JuMP.value(problem.w_zz[end]))
           # push!(w_uu[iScen],JuMP.value(problem.w_uu[end]))

            push!(w_xy[iScen],JuMP.value(problem.w_xy[end]))
            push!(w_xz[iScen],JuMP.value(problem.w_xz[end]))
            push!(w_zy[iScen],JuMP.value(problem.w_zy[end]))
            #=push!(w_xu[iScen],JuMP.value(problem.w_xu[end]))
            push!(w_yu[iScen],JuMP.value(problem.w_yu[end]))
            push!(w_zu[iScen],JuMP.value(problem.w_zu[end]))=#
        
            
            for iStage=1:NStages
                rev[iScen,iStage] = JuMP.value(problem.revamping[iStage])
                deg_stage[iScen, iStage] = sum(deg[iScen][iStep] for iStep=(vector_stages_index[iScen, iStage]+1):(vector_stages_index[iScen,iStage+1]))*k
            end
        
            for iStage=2:(NStages-1)
                net_revenues_per_stage[iScen, iStage] = sum(vector_prices[iScen][iStep]*NHoursStep*(discharge[iScen][iStep] - charge[iScen][iStep]) for iStep=(vector_stages_index[iScen, iStage]+1):(vector_stages_index[iScen,iStage+1])) - Battery_price_purchase[iStage]*rev[iScen, iStage] 
                WEM_stage[iScen, iStage] = sum(vector_prices[iScen][iStep]*NHoursStep*(discharge[iScen][iStep] - charge[iScen][iStep]) for iStep=(vector_stages_index[iScen, iStage]+1):(vector_stages_index[iScen, iStage+1]))
                cost_rev[iScen, iStage] = Battery_price_purchase[iStage]*(rev[iScen, iStage]) 

            end

            net_revenues_per_stage[iScen, 1] = sum(vector_prices[iScen][iStep]*NHoursStep*(discharge[iScen][iStep]-charge[iScen][iStep]) for iStep=(vector_stages_index[iScen, 1]+1):(vector_stages_index[iScen, 2])) - Battery_price_purchase[1]*rev[iScen, 1]  
            WEM_stage[iScen, 1] = sum(vector_prices[iScen][iStep]*NHoursStep*(discharge[iScen][iStep]-charge[iScen][iStep]) for iStep=(vector_stages_index[iScen, 1]+1):(vector_stages_index[iScen, 2]))
            cost_rev[iScen, 1] = Battery_price_purchase[1]*rev[iScen, 1] 

            net_revenues_per_stage[iScen, NStages] = sum(vector_prices[iScen][iStep]*NHoursStep*(discharge[iScen][iStep]-charge[iScen][iStep]) for iStep=(vector_stages_index[iScen, NStages]+1):(vector_stages_index[iScen, NStages+1])) + Battery_price_purchase[NStages+1]*(cap[iScen][end]-min_SOH) - Battery_price_purchase[NStages]*rev[iScen, NStages] 
            WEM_stage[iScen, NStages]= sum(vector_prices[iScen][iStep]*NHoursStep*(discharge[iScen][iStep]-charge[iScen][iStep]) for iStep=(vector_stages_index[iScen, NStages]+1):(vector_stages_index[iScen, NStages+1]))
            cost_rev[iScen, NStages] = Battery_price_purchase[NStages]*rev[iScen, NStages] - Battery_price_purchase[NStages+1]*(cap[iScen][end]-min_SOH)  
        end
        
        problem = 0
    end
    
    println("Collected results")

    return Results_4(
        objective,
        net_revenues_per_stage,
        WEM_stage,
        cost_rev,
        deg_stage,
        soc,
        charge,
        discharge,
        deg,
        soc_quad,
        x,
        y,
        z,
        #u,
        w_xx,
        w_yy,
        w_zz,
        #w_uu,
        w_xy,
        w_xz,
        w_zy,
        #w_xu,
        #w_yu,
        #w_zu,
        rev,
        cap,  
    )

end

# SOLVE OPTIMIZATION PROBLEM

function solveOptimizationProblem_5(InputParameters::InputParam, SolverParameters::SolverParam, Battery::BatteryParam, NScen, Results_ex_post::results_ex_post)

    @unpack (NYears, NMonths, NStages, Big, NHoursStep, bin) = InputParameters;                #NSteps, NHoursStage
    @unpack (min_SOC, max_SOC, Eff_charge, Eff_discharge, min_P, max_P, max_SOH, min_SOH, Nfull ) = Battery;
    @unpack (vector_prices, vector_stages_index, NSteps_scenario) = Results_ex_post;

    println("Solving Optimization Problem")

    k = min_SOH/(2*Nfull)

    # PER OGNI SCENARIO, C'E' SEMPRE LO STESSO NUMERO DI STAGES(semestri)
    objective = zeros(NScen);                 
    net_revenues_per_stage = zeros(NScen, NStages);
    WEM_stage = zeros(NScen, NStages);
    cost_rev = zeros(NScen, NStages);
    deg_stage = zeros(NScen, NStages);

    # PER OGNI SCENARIO, CAMBIA IL NUMERO DI STEPS
    charge = Vector{Vector{Float64}}(undef,NScen)
    discharge = Vector{Vector{Float64}}(undef,NScen)
    soc = Vector{Vector{Float64}}(undef,NScen)
    deg = Vector{Vector{Float64}}(undef,NScen)
    cap = Vector{Vector{Float64}}(undef,NScen)
    rev= zeros(NScen, NStages);

    soc_quad = Vector{Vector{Float64}}(undef,NScen)
    x = Vector{Vector{Float64}}(undef,NScen)
    y = Vector{Vector{Float64}}(undef,NScen)
    z = Vector{Vector{Float64}}(undef,NScen)
    u = Vector{Vector{Float64}}(undef,NScen)
    t = Vector{Vector{Float64}}(undef,NScen)

    w_xx = Vector{Vector{Float64}}(undef,NScen)
    w_yy = Vector{Vector{Float64}}(undef,NScen)
    w_zz = Vector{Vector{Float64}}(undef,NScen)
    w_uu = Vector{Vector{Float64}}(undef,NScen)
    w_tt = Vector{Vector{Float64}}(undef,NScen)

    w_xy = Vector{Vector{Float64}}(undef,NScen)
    w_xz = Vector{Vector{Float64}}(undef,NScen)
    w_zy = Vector{Vector{Float64}}(undef,NScen)
    w_xu = Vector{Vector{Float64}}(undef,NScen)
    w_yu = Vector{Vector{Float64}}(undef,NScen)
    w_zu = Vector{Vector{Float64}}(undef,NScen)
    w_tx = Vector{Vector{Float64}}(undef,NScen)
    w_ty = Vector{Vector{Float64}}(undef,NScen)
    w_tz = Vector{Vector{Float64}}(undef,NScen)
    w_tu = Vector{Vector{Float64}}(undef,NScen)

    problem = 0 #BuildStageProblem(InputParameters, SolverParameters, Battery)
    NSteps_sim= 0

    for iScen=1:NScen

        println("Solving scenario $iScen")
        NSteps_sim = Int(NSteps_scenario[iScen])
        problem = BuildStageProblem_5(InputParameters, SolverParameters, Battery, NSteps_sim, vector_stages_index, iScen, vector_prices)

        charge[iScen] = Float64[]
        discharge[iScen] = Float64[]
        soc[iScen] = Float64[]
        deg[iScen] = Float64[]
        cap[iScen] = Float64[]
        soc_quad[iScen] = Float64[]
        x[iScen] = Float64[]
        y[iScen] = Float64[]
        z[iScen] = Float64[]
        u[iScen] = Float64[]
        t[iScen] = Float64[]

        w_xx[iScen] = Float64[]
        w_yy[iScen] = Float64[]
        w_zz[iScen] = Float64[]
        w_uu[iScen] = Float64[]
        w_tt[iScen] = Float64[]

        w_xy[iScen] = Float64[]
        w_xz[iScen] = Float64[]
        w_zy[iScen] = Float64[]
        w_xu[iScen] = Float64[]
        w_yu[iScen] = Float64[]
        w_zu[iScen] = Float64[]
        w_tx[iScen] = Float64[]
        w_ty[iScen] = Float64[]
        w_tz[iScen] = Float64[]
        w_tu[iScen] = Float64[]

        @timeit to "Solve optimization" optimize!(problem.M)

        if termination_status(problem.M) != MOI.OPTIMAL
            println("NOT OPTIMAL: ", termination_status(problem.M))
        else
            println("Optimization finished")
        end

        @timeit to "Collecting results" begin
            objective = JuMP.objective_value(problem.M)
            
            for iStep=1:NSteps_sim
                
                push!(soc[iScen],JuMP.value(problem.soc[iStep]))
                push!(soc_quad[iScen],JuMP.value(problem.soc_quad[iStep]))
                push!(charge[iScen],JuMP.value(problem.charge[iStep]))
                push!(discharge[iScen],JuMP.value(problem.discharge[iStep]))
                push!(deg[iScen],JuMP.value(problem.deg[iStep]))
                push!(cap[iScen],JuMP.value(problem.capacity[iStep]))

                push!(x[iScen],JuMP.value(problem.x[iStep]))
                push!(y[iScen],JuMP.value(problem.y[iStep]))
                push!(z[iScen],JuMP.value(problem.z[iStep]))
                push!(u[iScen],JuMP.value(problem.u[iStep]))
                push!(t[iScen],JuMP.value(problem.t[iStep]))

                push!(w_xx[iScen],JuMP.value(problem.w_xx[iStep]))
                push!(w_yy[iScen],JuMP.value(problem.w_yy[iStep]))
                push!(w_zz[iScen],JuMP.value(problem.w_zz[iStep]))
                push!(w_uu[iScen],JuMP.value(problem.w_uu[iStep]))
                push!(w_tt[iScen],JuMP.value(problem.w_tt[iStep]))

                push!(w_xy[iScen],JuMP.value(problem.w_xy[iStep]))
                push!(w_xz[iScen],JuMP.value(problem.w_xz[iStep]))
                push!(w_zy[iScen],JuMP.value(problem.w_zy[iStep]))
                push!(w_xu[iScen],JuMP.value(problem.w_xu[iStep]))
                push!(w_yu[iScen],JuMP.value(problem.w_yu[iStep]))
                push!(w_zu[iScen],JuMP.value(problem.w_zu[iStep]))
                push!(w_tx[iScen],JuMP.value(problem.w_tx[iStep]))
                push!(w_ty[iScen],JuMP.value(problem.w_ty[iStep]))
                push!(w_tz[iScen],JuMP.value(problem.w_tz[iStep]))
                push!(w_tu[iScen],JuMP.value(problem.w_tu[iStep]))

            end
                    
            push!(soc[iScen],JuMP.value(problem.soc[end]))
            push!(soc_quad[iScen],JuMP.value(problem.soc_quad[end]))
            push!(cap[iScen],JuMP.value(problem.capacity[end]))

            push!(x[iScen],JuMP.value(problem.x[end]))
            push!(y[iScen],JuMP.value(problem.y[end]))
            push!(z[iScen],JuMP.value(problem.z[end]))
            push!(u[iScen],JuMP.value(problem.u[end]))
            push!(t[iScen],JuMP.value(problem.t[end]))

            push!(w_xx[iScen],JuMP.value(problem.w_xx[end]))
            push!(w_yy[iScen],JuMP.value(problem.w_yy[end]))
            push!(w_zz[iScen],JuMP.value(problem.w_zz[end]))
            push!(w_uu[iScen],JuMP.value(problem.w_uu[end]))
            push!(w_tt[iScen],JuMP.value(problem.w_tt[end]))

            push!(w_xy[iScen],JuMP.value(problem.w_xy[end]))
            push!(w_xz[iScen],JuMP.value(problem.w_xz[end]))
            push!(w_zy[iScen],JuMP.value(problem.w_zy[end]))
            push!(w_xu[iScen],JuMP.value(problem.w_xu[end]))
            push!(w_yu[iScen],JuMP.value(problem.w_yu[end]))
            push!(w_zu[iScen],JuMP.value(problem.w_zu[end]))
            push!(w_tx[iScen],JuMP.value(problem.w_tx[end]))
            push!(w_ty[iScen],JuMP.value(problem.w_ty[end]))
            push!(w_tz[iScen],JuMP.value(problem.w_tz[end]))
            push!(w_tu[iScen],JuMP.value(problem.w_tu[end]))
        
            
            for iStage=1:NStages
                rev[iScen,iStage] = JuMP.value(problem.revamping[iStage])
                deg_stage[iScen, iStage] = sum(deg[iScen][iStep] for iStep=(vector_stages_index[iScen, iStage]+1):(vector_stages_index[iScen,iStage+1]))*k
            end
        
            for iStage=2:(NStages-1)
                net_revenues_per_stage[iScen, iStage] = sum(vector_prices[iScen][iStep]*NHoursStep*(discharge[iScen][iStep] - charge[iScen][iStep]) for iStep=(vector_stages_index[iScen, iStage]+1):(vector_stages_index[iScen,iStage+1])) - Battery_price_purchase[iStage]*rev[iScen, iStage] 
                WEM_stage[iScen, iStage] = sum(vector_prices[iScen][iStep]*NHoursStep*(discharge[iScen][iStep] - charge[iScen][iStep]) for iStep=(vector_stages_index[iScen, iStage]+1):(vector_stages_index[iScen, iStage+1]))
                cost_rev[iScen, iStage] = Battery_price_purchase[iStage]*(rev[iScen, iStage]) 

            end

            net_revenues_per_stage[iScen, 1] = sum(vector_prices[iScen][iStep]*NHoursStep*(discharge[iScen][iStep]-charge[iScen][iStep]) for iStep=(vector_stages_index[iScen, 1]+1):(vector_stages_index[iScen, 2])) - Battery_price_purchase[1]*rev[iScen, 1]  
            WEM_stage[iScen, 1] = sum(vector_prices[iScen][iStep]*NHoursStep*(discharge[iScen][iStep]-charge[iScen][iStep]) for iStep=(vector_stages_index[iScen, 1]+1):(vector_stages_index[iScen, 2]))
            cost_rev[iScen, 1] = Battery_price_purchase[1]*rev[iScen, 1] 

            net_revenues_per_stage[iScen, NStages] = sum(vector_prices[iScen][iStep]*NHoursStep*(discharge[iScen][iStep]-charge[iScen][iStep]) for iStep=(vector_stages_index[iScen, NStages]+1):(vector_stages_index[iScen, NStages+1])) + Battery_price_purchase[NStages+1]*(cap[iScen][end]-min_SOH) - Battery_price_purchase[NStages]*rev[iScen, NStages] 
            WEM_stage[iScen, NStages]= sum(vector_prices[iScen][iStep]*NHoursStep*(discharge[iScen][iStep]-charge[iScen][iStep]) for iStep=(vector_stages_index[iScen, NStages]+1):(vector_stages_index[iScen, NStages+1]))
            cost_rev[iScen, NStages] = Battery_price_purchase[NStages]*rev[iScen, NStages] - Battery_price_purchase[NStages+1]*(cap[iScen][end]-min_SOH)  
        end
        
        problem = 0
    end
    
    println("Collected results")

    return Results_5(
        objective,
        net_revenues_per_stage,
        WEM_stage,
        cost_rev,
        deg_stage,
        soc,
        charge,
        discharge,
        deg,
        soc_quad,
        x,
        y,
        z,
        u,
        w_xx,
        w_yy,
        w_zz,
        w_uu,
        w_xy,
        w_xz,
        w_zy,
        w_xu,
        w_yu,
        w_zu,
        rev,
        cap,  
        t,
        w_tt,
        w_tx,
        w_ty,
        w_tz,
        w_tu,
    )

end

