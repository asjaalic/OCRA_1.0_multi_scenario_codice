#EX-POST ANALYSIS

function ex_post_analysis_3(Results_No_Deg_3::ResultWithoutDeg_3, InputParameters::InputParam, Battery::BatteryParam, NScen, NSteps, Pp, Steps_stages)

    @unpack (charge_no_deg, discharge_no_deg, soc_no_deg, soc_quad_no_deg, revenues_per_stage_no_deg) = Results_No_Deg_3
    @unpack (NYears, NMonths, NHoursStep, NStages, Big, conv, bin) = InputParameters
    @unpack (min_SOC, max_SOC, Eff_charge, Eff_discharge, min_P, max_P, max_SOH, min_SOH, Nfull, downtime) = Battery

    deg_singola = zeros(NScen, NSteps)
    deg_stage = zeros(NScen, NStages)
    costo_stage = zeros(NScen, NStages)
    net_revenues = zeros(NScen, NStages)
    NSteps_scenario = zeros(NScen)

    #CALCOLO DEGRADAZIONE SINGOLA, TOTALE, COSTI e NET REVENUES
    for iScen=1:NScen
        for iStep=1:NSteps
            deg_singola[iScen, iStep] = min_SOH/(2*Nfull)*abs((soc_quad_no_deg[iScen, iStep]-soc_quad_no_deg[iScen, iStep+1])/max_SOC^2+2/max_SOC*(soc_no_deg[iScen, iStep+1]-soc_no_deg[iScen, iStep]))
        end
        for iStage=1:NStages
            deg_stage[iScen, iStage] = sum(deg_singola[iScen,iStep] for iStep = Steps_stages[iStage]+1:Steps_stages[iStage+1] )
            costo_stage[iScen, iStage] = deg_stage[iScen, iStage]*Battery_price_purchase[iStage]
            net_revenues[iScen, iStage] = revenues_per_stage_no_deg[iScen, iStage] - costo_stage[iScen, iStage]
        end
    end

    # PRICE-FILTERING STRATEGY
    vector_prices = Vector{Vector{Float64}}(undef,NScen)
    vector_stages_index = fill(0, NScen, NStages+1)
    vector_downtime_stages = fill(0, NScen, NStages) 
    cont_two= 0

    for iScen=1:NScen
        vector_prices[iScen] = Float64[]
        cont = 0
        
        # PER OGNI STAGE (SEMESTRE O ANNO)
        for iStage=1:NStages
            
            # FACCIO LA PRICE FILTERING
            for iStep = Steps_stages[iStage]+1:Steps_stages[iStage+1]

                if charge_no_deg[iScen, iStep]!=0 || discharge_no_deg[iScen, iStep] !=0         # Salvo i prezzi dell'energia per cui la batteria effettivamente fa un ciclo di carica o scarica
                    cont+=1
                    push!(vector_prices[iScen], Pp[iStep, iScen])
                end
            end
            vector_stages_index[iScen,iStage+1] = Int(cont)

            # CREO IL VETTORE CON IL NUMERO DI MOMENTI CON CUI LA BATTERIA DEVE STARE FERMA -> downtime period non è valido per il primo semestre
            for iStep=Steps_stages[iStage]+1:Steps_stages[iStage]+downtime
                if charge_no_deg[iScen, iStep]!=0 || discharge_no_deg[iScen, iStep] !=0
                    cont_two +=1                  
                end
            end
            println("Downtime scenario $iScen stage $iStage: $cont_two")
            vector_downtime_stages[iScen,iStage]=Int(cont_two)
            cont_two= 0
        end
        NSteps_scenario[iScen]= vector_stages_index[iScen,end]
    end

   
    
    return results_ex_post(
        deg_singola,
        deg_stage,
        costo_stage,
        net_revenues,
        vector_prices,
        vector_stages_index,
        NSteps_scenario,
        vector_downtime_stages,
    )
end

function statistical_analysis_3(ResultsOpt_3::Results_3, Results_ex_post::results_ex_post, NScen, NStages)

    @unpack (net_revenues_per_stage, WEM_stage, cost_rev, deg_stage, rev, cap) = ResultsOpt_3;
    @unpack (vector_prices, vector_stages_index) = Results_ex_post;

    Initial_capacity = zeros(NScen, NStages);
    Final_capacity = zeros(NScen, NStages);

    for iScen=1:NScen
        for iStage=1:NStages
            Initial_capacity[iScen, iStage] = cap[iScen][2]
            Final_capacity[iScen, iStage] = cap[iScen][end]
        end
    end

    #1=max, 2= min, 3= avg, 4 =median
    WEM_rev_stat = zeros(NStages, 4)
    cost_stat = zeros(NStages, 4)
    net_rev_stat = zeros(NStages, 4)
    deg_stat = zeros(NStages, 4)
    Initial_cap_stat = zeros(NStages, 4)
    Final_cap_stat = zeros(NStages, 4)

    for iStage=1:NStages

        # Calcolo valori massimi su NScen scenari
        WEM_rev_stat[iStage,1] = findmax(WEM_stage[:,iStage])[1]
        cost_stat[iStage,1] = findmax(cost_rev[:,iStage])[1]
        net_rev_stat[iStage,1] = findmax(net_revenues_per_stage[:,iStage])[1]
        deg_stat[iStage,1] = findmax(deg_stage[:,iStage])[1]
        Initial_cap_stat[iStage,1] = findmax(Initial_capacity[:,iStage])[1]
        Final_cap_stat[iStage,1] = findmax(Final_capacity[:,iStage])[1]

        # Calcolo valori MINIMI su NScen scenari
        WEM_rev_stat[iStage,2] = findmin(WEM_stage[:,iStage])[1]
        cost_stat[iStage,2] = findmin(cost_rev[:,iStage])[1]
        net_rev_stat[iStage,2] = findmin(net_revenues_per_stage[:,iStage])[1]
        deg_stat[iStage,2] = findmin(deg_stage[:,iStage])[1]
        Initial_cap_stat[iStage,2] = findmin(Initial_capacity[:,iStage])[1]
        Final_cap_stat[iStage,2] = findmin(Final_capacity[:,iStage])[1]

        # Calcolo valori MEDI su NScen scenari
        WEM_rev_stat[iStage,3] = mean(WEM_stage[:,iStage])
        cost_stat[iStage,3] = mean(cost_rev[:,iStage])
        net_rev_stat[iStage,3] = mean(net_revenues_per_stage[:,iStage])
        deg_stat[iStage,3] = mean(deg_stage[:,iStage])
        Initial_cap_stat[iStage,3] = mean(Initial_capacity[:,iStage])
        Final_cap_stat[iStage,3] = mean(Final_capacity[:,iStage])

        # Calcolo valori MEDIANI su NScen scenari
        WEM_rev_stat[iStage,4] = median(WEM_stage[:,iStage])
        cost_stat[iStage,4] = median(cost_rev[:,iStage])
        net_rev_stat[iStage,4] = median(net_revenues_per_stage[:,iStage])
        deg_stat[iStage,4] = median(deg_stage[:,iStage])
        Initial_cap_stat[iStage,4] = median(Initial_capacity[:,iStage])
        Final_cap_stat[iStage,4] = median(Final_capacity[:,iStage])
    end

    WEM_rev_stat_frame = DataFrame()
    cost_stat_frame = DataFrame()
    net_rev_stat_frame = DataFrame()
    deg_stat_frame = DataFrame()
    Initial_cap_stat_frame = DataFrame()
    Final_cap_stat_frame = DataFrame()

    WEM_rev_stat_frame[!,"Stages"] = 1:1:NStages
    WEM_rev_stat_frame[!,"Max"] = WEM_rev_stat[:,1]
    WEM_rev_stat_frame[!,"Min"] = WEM_rev_stat[:,2]
    WEM_rev_stat_frame[!,"Average"] = WEM_rev_stat[:,3]
    WEM_rev_stat_frame[!,"Median"] = WEM_rev_stat[:,4]

    cost_stat_frame[!,"Stages"] = 1:1:NStages
    cost_stat_frame[!,"Max"] = cost_stat[:,1]
    cost_stat_frame[!,"Min"] = cost_stat[:,2]
    cost_stat_frame[!,"Average"] = cost_stat[:,3]
    cost_stat_frame[!,"Median"] = cost_stat[:,4]

    net_rev_stat_frame[!,"Stages"] = 1:1:NStages
    net_rev_stat_frame[!,"Max"] = net_rev_stat[:,1]
    net_rev_stat_frame[!,"Min"] = net_rev_stat[:,2]
    net_rev_stat_frame[!,"Average"] = net_rev_stat[:,3]
    net_rev_stat_frame[!,"Median"] = net_rev_stat[:,4]

    deg_stat_frame[!,"Stages"] = 1:1:NStages
    deg_stat_frame[!, "Max"] = deg_stat[:,1]
    deg_stat_frame[!, "Min"] = deg_stat[:,2]
    deg_stat_frame[!, "Average"] = deg_stat[:,3]
    deg_stat_frame[!, "Median"] = deg_stat[:,4]

    Initial_cap_stat_frame[!, "Stages"] =1:1:NStages
    Initial_cap_stat_frame[!, "Max"] = Initial_cap_stat[:,1]
    Initial_cap_stat_frame[!, "Min"] = Initial_cap_stat[:,2]
    Initial_cap_stat_frame[!, "Average"] = Initial_cap_stat[:,3]
    Initial_cap_stat_frame[!, "Median"] = Initial_cap_stat[:,4]

    Final_cap_stat_frame[!, "Stages"] =1:1:NStages
    Final_cap_stat_frame[!, "Max"] = Final_cap_stat[:,1]
    Final_cap_stat_frame[!, "Min"] = Final_cap_stat[:,2]
    Final_cap_stat_frame[!, "Average"] = Final_cap_stat[:,3]
    Final_cap_stat_frame[!, "Median"] = Final_cap_stat[:,4]

    return Data_analysed(
        WEM_rev_stat_frame,
        cost_stat_frame,
        net_rev_stat_frame,
        deg_stat_frame,
        Initial_cap_stat_frame,
        Final_cap_stat_frame,
    )
end

function ex_post_analysis_4(Results_No_Deg_4::ResultWithoutDeg_4, InputParameters::InputParam, Battery::BatteryParam, NScen, NSteps, Pp, Steps_stages)

    @unpack (charge_no_deg, discharge_no_deg, soc_no_deg, soc_quad_no_deg, revenues_per_stage_no_deg) = Results_No_Deg_4
    @unpack (NYears, NMonths, NHoursStep, NStages, Big, conv, bin) = InputParameters
    @unpack (min_SOC, max_SOC, Eff_charge, Eff_discharge, min_P, max_P, max_SOH, min_SOH, Nfull, downtime) = Battery

    deg_singola = zeros(NScen, NSteps)
    deg_stage = zeros(NScen, NStages)
    costo_stage = zeros(NScen, NStages)
    net_revenues = zeros(NScen, NStages)
    NSteps_scenario = zeros(NScen)

    #CALCOLO DEGRADAZIONE SINGOLA, TOTALE, COSTI e NET REVENUES
    for iScen=1:NScen
        for iStep=1:NSteps
            deg_singola[iScen, iStep] = min_SOH/(2*Nfull)*abs((soc_quad_no_deg[iScen, iStep]-soc_quad_no_deg[iScen, iStep+1])/max_SOC^2+2/max_SOC*(soc_no_deg[iScen, iStep+1]-soc_no_deg[iScen, iStep]))
        end
        for iStage=1:NStages
            deg_stage[iScen, iStage] = sum(deg_singola[iScen,iStep] for iStep = Steps_stages[iStage]+1:Steps_stages[iStage+1] )
            costo_stage[iScen, iStage] = deg_stage[iScen, iStage]*Battery_price_purchase[iStage]
            net_revenues[iScen, iStage] = revenues_per_stage_no_deg[iScen, iStage] - costo_stage[iScen, iStage]
        end
    end

    # PRICE-FILTERING STRATEGY
    vector_prices = Vector{Vector{Float64}}(undef,NScen)
    vector_stages_index = fill(0, NScen, NStages+1)
    vector_downtime_stages = fill(0, NScen, NStages) 
    cont_two= 0

    for iScen=1:NScen
        vector_prices[iScen] = Float64[]
        cont = 0
        
        # PER OGNI STAGE (SEMESTRE O ANNO)
        for iStage=1:NStages
            
            # FACCIO LA PRICE FILTERING
            for iStep = Steps_stages[iStage]+1:Steps_stages[iStage+1]

                if charge_no_deg[iScen, iStep]!=0 || discharge_no_deg[iScen, iStep] !=0         # Salvo i prezzi dell'energia per cui la batteria effettivamente fa un ciclo di carica o scarica
                    cont+=1
                    push!(vector_prices[iScen], Pp[iStep, iScen])
                end
            end
            vector_stages_index[iScen,iStage+1] = Int(cont)

            # CREO IL VETTORE CON IL NUMERO DI MOMENTI CON CUI LA BATTERIA DEVE STARE FERMA -> downtime period non è valido per il primo semestre
            for iStep=Steps_stages[iStage]+1:Steps_stages[iStage]+downtime
                if charge_no_deg[iScen, iStep]!=0 || discharge_no_deg[iScen, iStep] !=0
                    cont_two +=1                  
                end
            end
            println("Downtime scenario $iScen stage $iStage: $cont_two")
            vector_downtime_stages[iScen,iStage]=Int(cont_two)
            cont_two= 0
        end
        NSteps_scenario[iScen]= vector_stages_index[iScen,end]
    end

    return results_ex_post(
        deg_singola,
        deg_stage,
        costo_stage,
        net_revenues,
        vector_prices,
        vector_stages_index,
        NSteps_scenario,
        vector_downtime_stages,
    )
end

function statistical_analysis_4(ResultsOpt_4::Results_4, Results_ex_post::results_ex_post, NScen, NStages)

    @unpack (net_revenues_per_stage, WEM_stage, cost_rev, deg_stage, rev, cap) = ResultsOpt_4;
    @unpack (vector_prices, vector_stages_index) = Results_ex_post;

    Initial_capacity = zeros(NScen, NStages);
    Final_capacity = zeros(NScen, NStages);

    for iScen=1:NScen
        for iStage=1:NStages
            Initial_capacity[iScen, iStage] = cap[iScen][2]
            Final_capacity[iScen, iStage] = cap[iScen][end]
        end
    end

    #1=max, 2= min, 3= avg, 4 =median
    WEM_rev_stat = zeros(NStages, 4)
    cost_stat = zeros(NStages, 4)
    net_rev_stat = zeros(NStages, 4)
    deg_stat = zeros(NStages, 4)
    Initial_cap_stat = zeros(NStages, 4)
    Final_cap_stat = zeros(NStages, 4)

    for iStage=1:NStages

        # Calcolo valori massimi su NScen scenari
        WEM_rev_stat[iStage,1] = findmax(WEM_stage[:,iStage])[1]
        cost_stat[iStage,1] = findmax(cost_rev[:,iStage])[1]
        net_rev_stat[iStage,1] = findmax(net_revenues_per_stage[:,iStage])[1]
        deg_stat[iStage,1] = findmax(deg_stage[:,iStage])[1]
        Initial_cap_stat[iStage,1] = findmax(Initial_capacity[:,iStage])[1]
        Final_cap_stat[iStage,1] = findmax(Final_capacity[:,iStage])[1]

        # Calcolo valori MINIMI su NScen scenari
        WEM_rev_stat[iStage,2] = findmin(WEM_stage[:,iStage])[1]
        cost_stat[iStage,2] = findmin(cost_rev[:,iStage])[1]
        net_rev_stat[iStage,2] = findmin(net_revenues_per_stage[:,iStage])[1]
        deg_stat[iStage,2] = findmin(deg_stage[:,iStage])[1]
        Initial_cap_stat[iStage,2] = findmin(Initial_capacity[:,iStage])[1]
        Final_cap_stat[iStage,2] = findmin(Final_capacity[:,iStage])[1]

        # Calcolo valori MEDI su NScen scenari
        WEM_rev_stat[iStage,3] = mean(WEM_stage[:,iStage])
        cost_stat[iStage,3] = mean(cost_rev[:,iStage])
        net_rev_stat[iStage,3] = mean(net_revenues_per_stage[:,iStage])
        deg_stat[iStage,3] = mean(deg_stage[:,iStage])
        Initial_cap_stat[iStage,3] = mean(Initial_capacity[:,iStage])
        Final_cap_stat[iStage,3] = mean(Final_capacity[:,iStage])

        # Calcolo valori MEDIANI su NScen scenari
        WEM_rev_stat[iStage,4] = median(WEM_stage[:,iStage])
        cost_stat[iStage,4] = median(cost_rev[:,iStage])
        net_rev_stat[iStage,4] = median(net_revenues_per_stage[:,iStage])
        deg_stat[iStage,4] = median(deg_stage[:,iStage])
        Initial_cap_stat[iStage,4] = median(Initial_capacity[:,iStage])
        Final_cap_stat[iStage,4] = median(Final_capacity[:,iStage])
    end

    WEM_rev_stat_frame = DataFrame()
    cost_stat_frame = DataFrame()
    net_rev_stat_frame = DataFrame()
    deg_stat_frame = DataFrame()
    Initial_cap_stat_frame = DataFrame()
    Final_cap_stat_frame = DataFrame()

    WEM_rev_stat_frame[!,"Stages"] = 1:1:NStages
    WEM_rev_stat_frame[!,"Max"] = WEM_rev_stat[:,1]
    WEM_rev_stat_frame[!,"Min"] = WEM_rev_stat[:,2]
    WEM_rev_stat_frame[!,"Average"] = WEM_rev_stat[:,3]
    WEM_rev_stat_frame[!,"Median"] = WEM_rev_stat[:,4]

    cost_stat_frame[!,"Stages"] = 1:1:NStages
    cost_stat_frame[!,"Max"] = cost_stat[:,1]
    cost_stat_frame[!,"Min"] = cost_stat[:,2]
    cost_stat_frame[!,"Average"] = cost_stat[:,3]
    cost_stat_frame[!,"Median"] = cost_stat[:,4]

    net_rev_stat_frame[!,"Stages"] = 1:1:NStages
    net_rev_stat_frame[!,"Max"] = net_rev_stat[:,1]
    net_rev_stat_frame[!,"Min"] = net_rev_stat[:,2]
    net_rev_stat_frame[!,"Average"] = net_rev_stat[:,3]
    net_rev_stat_frame[!,"Median"] = net_rev_stat[:,4]

    deg_stat_frame[!,"Stages"] = 1:1:NStages
    deg_stat_frame[!, "Max"] = deg_stat[:,1]
    deg_stat_frame[!, "Min"] = deg_stat[:,2]
    deg_stat_frame[!, "Average"] = deg_stat[:,3]
    deg_stat_frame[!, "Median"] = deg_stat[:,4]

    Initial_cap_stat_frame[!, "Stages"] =1:1:NStages
    Initial_cap_stat_frame[!, "Max"] = Initial_cap_stat[:,1]
    Initial_cap_stat_frame[!, "Min"] = Initial_cap_stat[:,2]
    Initial_cap_stat_frame[!, "Average"] = Initial_cap_stat[:,3]
    Initial_cap_stat_frame[!, "Median"] = Initial_cap_stat[:,4]

    Final_cap_stat_frame[!, "Stages"] =1:1:NStages
    Final_cap_stat_frame[!, "Max"] = Final_cap_stat[:,1]
    Final_cap_stat_frame[!, "Min"] = Final_cap_stat[:,2]
    Final_cap_stat_frame[!, "Average"] = Final_cap_stat[:,3]
    Final_cap_stat_frame[!, "Median"] = Final_cap_stat[:,4]

    return Data_analysed(
        WEM_rev_stat_frame,
        cost_stat_frame,
        net_rev_stat_frame,
        deg_stat_frame,
        Initial_cap_stat_frame,
        Final_cap_stat_frame,
    )
end


function ex_post_analysis_5(Results_No_Deg_5::ResultWithoutDeg_5, InputParameters::InputParam, Battery::BatteryParam, NScen, NSteps, Pp, Steps_stages)

    @unpack (charge_no_deg, discharge_no_deg, soc_no_deg, soc_quad_no_deg, revenues_per_stage_no_deg) = Results_No_Deg_5
    @unpack (NYears, NMonths, NHoursStep, NStages, Big, conv, bin) = InputParameters
    @unpack (min_SOC, max_SOC, Eff_charge, Eff_discharge, min_P, max_P, max_SOH, min_SOH, Nfull, downtime) = Battery

    deg_singola = zeros(NScen, NSteps)
    deg_stage = zeros(NScen, NStages)
    costo_stage = zeros(NScen, NStages)
    net_revenues = zeros(NScen, NStages)
    NSteps_scenario = zeros(NScen)

    #CALCOLO DEGRADAZIONE SINGOLA, TOTALE, COSTI e NET REVENUES
    for iScen=1:NScen
        for iStep=1:NSteps
            deg_singola[iScen, iStep] = min_SOH/(2*Nfull)*abs((soc_quad_no_deg[iScen, iStep]-soc_quad_no_deg[iScen, iStep+1])/max_SOC^2+2/max_SOC*(soc_no_deg[iScen, iStep+1]-soc_no_deg[iScen, iStep]))
        end
        for iStage=1:NStages
            deg_stage[iScen, iStage] = sum(deg_singola[iScen,iStep] for iStep = Steps_stages[iStage]+1:Steps_stages[iStage+1] )
            costo_stage[iScen, iStage] = deg_stage[iScen, iStage]*Battery_price_purchase[iStage]
            net_revenues[iScen, iStage] = revenues_per_stage_no_deg[iScen, iStage] - costo_stage[iScen, iStage]
        end
    end

    # PRICE-FILTERING STRATEGY
    vector_prices = Vector{Vector{Float64}}(undef,NScen)
    vector_stages_index = fill(0, NScen, NStages+1)
    vector_downtime_stages = fill(0, NScen, NStages) 
    cont_two= 0

    for iScen=1:NScen
        vector_prices[iScen] = Float64[]
        cont = 0
        
        # PER OGNI STAGE (SEMESTRE O ANNO)
        for iStage=1:NStages
            
            # FACCIO LA PRICE FILTERING
            for iStep = Steps_stages[iStage]+1:Steps_stages[iStage+1]

                if charge_no_deg[iScen, iStep]!=0 || discharge_no_deg[iScen, iStep] !=0         # Salvo i prezzi dell'energia per cui la batteria effettivamente fa un ciclo di carica o scarica
                    cont+=1
                    push!(vector_prices[iScen], Pp[iStep, iScen])
                end
            end
            vector_stages_index[iScen,iStage+1] = Int(cont)

            # CREO IL VETTORE CON IL NUMERO DI MOMENTI CON CUI LA BATTERIA DEVE STARE FERMA -> downtime period non è valido per il primo semestre
            for iStep=Steps_stages[iStage]+1:Steps_stages[iStage]+downtime
                if charge_no_deg[iScen, iStep]!=0 || discharge_no_deg[iScen, iStep] !=0
                    cont_two +=1                  
                end
            end
            println("Downtime scenario $iScen stage $iStage: $cont_two")
            vector_downtime_stages[iScen,iStage]=Int(cont_two)
            cont_two= 0
        end
        NSteps_scenario[iScen]= vector_stages_index[iScen,end]
    end

    return results_ex_post(
        deg_singola,
        deg_stage,
        costo_stage,
        net_revenues,
        vector_prices,
        vector_stages_index,
        NSteps_scenario,
        vector_downtime_stages,
    )
end

function statistical_analysis_5(ResultsOpt_5::Results_5, Results_ex_post::results_ex_post, NScen, NStages)

    @unpack (net_revenues_per_stage, WEM_stage, cost_rev, deg_stage, rev, cap) = ResultsOpt_5;
    @unpack (vector_prices, vector_stages_index) = Results_ex_post;

    Initial_capacity = zeros(NScen, NStages);
    Final_capacity = zeros(NScen, NStages);

    for iScen=1:NScen
        for iStage=1:NStages
            Initial_capacity[iScen, iStage] = cap[iScen][2]
            Final_capacity[iScen, iStage] = cap[iScen][end]
        end
    end

    #1=max, 2= min, 3= avg, 4 =median
    WEM_rev_stat = zeros(NStages, 4)
    cost_stat = zeros(NStages, 4)
    net_rev_stat = zeros(NStages, 4)
    deg_stat = zeros(NStages, 4)
    Initial_cap_stat = zeros(NStages, 4)
    Final_cap_stat = zeros(NStages, 4)

    for iStage=1:NStages

        # Calcolo valori massimi su NScen scenari
        WEM_rev_stat[iStage,1] = findmax(WEM_stage[:,iStage])[1]
        cost_stat[iStage,1] = findmax(cost_rev[:,iStage])[1]
        net_rev_stat[iStage,1] = findmax(net_revenues_per_stage[:,iStage])[1]
        deg_stat[iStage,1] = findmax(deg_stage[:,iStage])[1]
        Initial_cap_stat[iStage,1] = findmax(Initial_capacity[:,iStage])[1]
        Final_cap_stat[iStage,1] = findmax(Final_capacity[:,iStage])[1]

        # Calcolo valori MINIMI su NScen scenari
        WEM_rev_stat[iStage,2] = findmin(WEM_stage[:,iStage])[1]
        cost_stat[iStage,2] = findmin(cost_rev[:,iStage])[1]
        net_rev_stat[iStage,2] = findmin(net_revenues_per_stage[:,iStage])[1]
        deg_stat[iStage,2] = findmin(deg_stage[:,iStage])[1]
        Initial_cap_stat[iStage,2] = findmin(Initial_capacity[:,iStage])[1]
        Final_cap_stat[iStage,2] = findmin(Final_capacity[:,iStage])[1]

        # Calcolo valori MEDI su NScen scenari
        WEM_rev_stat[iStage,3] = mean(WEM_stage[:,iStage])
        cost_stat[iStage,3] = mean(cost_rev[:,iStage])
        net_rev_stat[iStage,3] = mean(net_revenues_per_stage[:,iStage])
        deg_stat[iStage,3] = mean(deg_stage[:,iStage])
        Initial_cap_stat[iStage,3] = mean(Initial_capacity[:,iStage])
        Final_cap_stat[iStage,3] = mean(Final_capacity[:,iStage])

        # Calcolo valori MEDIANI su NScen scenari
        WEM_rev_stat[iStage,4] = median(WEM_stage[:,iStage])
        cost_stat[iStage,4] = median(cost_rev[:,iStage])
        net_rev_stat[iStage,4] = median(net_revenues_per_stage[:,iStage])
        deg_stat[iStage,4] = median(deg_stage[:,iStage])
        Initial_cap_stat[iStage,4] = median(Initial_capacity[:,iStage])
        Final_cap_stat[iStage,4] = median(Final_capacity[:,iStage])
    end

    WEM_rev_stat_frame = DataFrame()
    cost_stat_frame = DataFrame()
    net_rev_stat_frame = DataFrame()
    deg_stat_frame = DataFrame()
    Initial_cap_stat_frame = DataFrame()
    Final_cap_stat_frame = DataFrame()

    WEM_rev_stat_frame[!,"Stages"] = 1:1:NStages
    WEM_rev_stat_frame[!,"Max"] = WEM_rev_stat[:,1]
    WEM_rev_stat_frame[!,"Min"] = WEM_rev_stat[:,2]
    WEM_rev_stat_frame[!,"Average"] = WEM_rev_stat[:,3]
    WEM_rev_stat_frame[!,"Median"] = WEM_rev_stat[:,4]

    cost_stat_frame[!,"Stages"] = 1:1:NStages
    cost_stat_frame[!,"Max"] = cost_stat[:,1]
    cost_stat_frame[!,"Min"] = cost_stat[:,2]
    cost_stat_frame[!,"Average"] = cost_stat[:,3]
    cost_stat_frame[!,"Median"] = cost_stat[:,4]

    net_rev_stat_frame[!,"Stages"] = 1:1:NStages
    net_rev_stat_frame[!,"Max"] = net_rev_stat[:,1]
    net_rev_stat_frame[!,"Min"] = net_rev_stat[:,2]
    net_rev_stat_frame[!,"Average"] = net_rev_stat[:,3]
    net_rev_stat_frame[!,"Median"] = net_rev_stat[:,4]

    deg_stat_frame[!,"Stages"] = 1:1:NStages
    deg_stat_frame[!, "Max"] = deg_stat[:,1]
    deg_stat_frame[!, "Min"] = deg_stat[:,2]
    deg_stat_frame[!, "Average"] = deg_stat[:,3]
    deg_stat_frame[!, "Median"] = deg_stat[:,4]

    Initial_cap_stat_frame[!, "Stages"] =1:1:NStages
    Initial_cap_stat_frame[!, "Max"] = Initial_cap_stat[:,1]
    Initial_cap_stat_frame[!, "Min"] = Initial_cap_stat[:,2]
    Initial_cap_stat_frame[!, "Average"] = Initial_cap_stat[:,3]
    Initial_cap_stat_frame[!, "Median"] = Initial_cap_stat[:,4]

    Final_cap_stat_frame[!, "Stages"] =1:1:NStages
    Final_cap_stat_frame[!, "Max"] = Final_cap_stat[:,1]
    Final_cap_stat_frame[!, "Min"] = Final_cap_stat[:,2]
    Final_cap_stat_frame[!, "Average"] = Final_cap_stat[:,3]
    Final_cap_stat_frame[!, "Median"] = Final_cap_stat[:,4]

    return Data_analysed(
        WEM_rev_stat_frame,
        cost_stat_frame,
        net_rev_stat_frame,
        deg_stat_frame,
        Initial_cap_stat_frame,
        Final_cap_stat_frame,
    )
end



function analysis_OCRA_2(Results_OCRA_2::Results_OCRA2, Results_ex_post::results_ex_post, NScen, NStages)

    @unpack (net_revenues_per_stage, WEM_stage, cost_rev, deg_stage, rev, cap) = Results_OCRA_2;
    @unpack (vector_prices, vector_stages_index) = Results_ex_post;

    Initial_capacity = zeros(NScen, NStages);
    Final_capacity = zeros(NScen, NStages);

    for iScen=1:NScen
        for iStage=1:NStages
            Initial_capacity[iScen, iStage] = cap[iScen][2]
            Final_capacity[iScen, iStage] = cap[iScen][end]
        end
    end

    #1=max, 2= min, 3= avg, 4 =median
    WEM_rev_stat = zeros(NStages, 4)
    cost_stat = zeros(NStages, 4)
    net_rev_stat = zeros(NStages, 4)
    deg_stat = zeros(NStages, 4)
    Initial_cap_stat = zeros(NStages, 4)
    Final_cap_stat = zeros(NStages, 4)

    for iStage=1:NStages

        # Calcolo valori massimi su NScen scenari
        WEM_rev_stat[iStage,1] = findmax(WEM_stage[:,iStage])[1]
        cost_stat[iStage,1] = findmax(cost_rev[:,iStage])[1]
        net_rev_stat[iStage,1] = findmax(net_revenues_per_stage[:,iStage])[1]
        deg_stat[iStage,1] = findmax(deg_stage[:,iStage])[1]
        Initial_cap_stat[iStage,1] = findmax(Initial_capacity[:,iStage])[1]
        Final_cap_stat[iStage,1] = findmax(Final_capacity[:,iStage])[1]

        # Calcolo valori MINIMI su NScen scenari
        WEM_rev_stat[iStage,2] = findmin(WEM_stage[:,iStage])[1]
        cost_stat[iStage,2] = findmin(cost_rev[:,iStage])[1]
        net_rev_stat[iStage,2] = findmin(net_revenues_per_stage[:,iStage])[1]
        deg_stat[iStage,2] = findmin(deg_stage[:,iStage])[1]
        Initial_cap_stat[iStage,2] = findmin(Initial_capacity[:,iStage])[1]
        Final_cap_stat[iStage,2] = findmin(Final_capacity[:,iStage])[1]

        # Calcolo valori MEDI su NScen scenari
        WEM_rev_stat[iStage,3] = mean(WEM_stage[:,iStage])
        cost_stat[iStage,3] = mean(cost_rev[:,iStage])
        net_rev_stat[iStage,3] = mean(net_revenues_per_stage[:,iStage])
        deg_stat[iStage,3] = mean(deg_stage[:,iStage])
        Initial_cap_stat[iStage,3] = mean(Initial_capacity[:,iStage])
        Final_cap_stat[iStage,3] = mean(Final_capacity[:,iStage])

        # Calcolo valori MEDIANI su NScen scenari
        WEM_rev_stat[iStage,4] = median(WEM_stage[:,iStage])
        cost_stat[iStage,4] = median(cost_rev[:,iStage])
        net_rev_stat[iStage,4] = median(net_revenues_per_stage[:,iStage])
        deg_stat[iStage,4] = median(deg_stage[:,iStage])
        Initial_cap_stat[iStage,4] = median(Initial_capacity[:,iStage])
        Final_cap_stat[iStage,4] = median(Final_capacity[:,iStage])
    end

    WEM_rev_stat_frame = DataFrame()
    cost_stat_frame = DataFrame()
    net_rev_stat_frame = DataFrame()
    deg_stat_frame = DataFrame()
    Initial_cap_stat_frame = DataFrame()
    Final_cap_stat_frame = DataFrame()

    WEM_rev_stat_frame[!,"Stages"] = 1:1:NStages
    WEM_rev_stat_frame[!,"Max"] = WEM_rev_stat[:,1]
    WEM_rev_stat_frame[!,"Min"] = WEM_rev_stat[:,2]
    WEM_rev_stat_frame[!,"Average"] = WEM_rev_stat[:,3]
    WEM_rev_stat_frame[!,"Median"] = WEM_rev_stat[:,4]

    cost_stat_frame[!,"Stages"] = 1:1:NStages
    cost_stat_frame[!,"Max"] = cost_stat[:,1]
    cost_stat_frame[!,"Min"] = cost_stat[:,2]
    cost_stat_frame[!,"Average"] = cost_stat[:,3]
    cost_stat_frame[!,"Median"] = cost_stat[:,4]

    net_rev_stat_frame[!,"Stages"] = 1:1:NStages
    net_rev_stat_frame[!,"Max"] = net_rev_stat[:,1]
    net_rev_stat_frame[!,"Min"] = net_rev_stat[:,2]
    net_rev_stat_frame[!,"Average"] = net_rev_stat[:,3]
    net_rev_stat_frame[!,"Median"] = net_rev_stat[:,4]

    deg_stat_frame[!,"Stages"] = 1:1:NStages
    deg_stat_frame[!, "Max"] = deg_stat[:,1]
    deg_stat_frame[!, "Min"] = deg_stat[:,2]
    deg_stat_frame[!, "Average"] = deg_stat[:,3]
    deg_stat_frame[!, "Median"] = deg_stat[:,4]

    Initial_cap_stat_frame[!, "Stages"] =1:1:NStages
    Initial_cap_stat_frame[!, "Max"] = Initial_cap_stat[:,1]
    Initial_cap_stat_frame[!, "Min"] = Initial_cap_stat[:,2]
    Initial_cap_stat_frame[!, "Average"] = Initial_cap_stat[:,3]
    Initial_cap_stat_frame[!, "Median"] = Initial_cap_stat[:,4]

    Final_cap_stat_frame[!, "Stages"] =1:1:NStages
    Final_cap_stat_frame[!, "Max"] = Final_cap_stat[:,1]
    Final_cap_stat_frame[!, "Min"] = Final_cap_stat[:,2]
    Final_cap_stat_frame[!, "Average"] = Final_cap_stat[:,3]
    Final_cap_stat_frame[!, "Median"] = Final_cap_stat[:,4]

    return Data_analysed(
        WEM_rev_stat_frame,
        cost_stat_frame,
        net_rev_stat_frame,
        deg_stat_frame,
        Initial_cap_stat_frame,
        Final_cap_stat_frame,
    )
end

