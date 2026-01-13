function data_saving_3(InputParameters::InputParam, ResultsOpt_3::Results_3, Results_ex_post::results_ex_post, Results_statistics::Data_analysed, Battery::BatteryParam,  NScen)

    @unpack (NYears, NMonths, NStages, Big, NHoursStep, bin) = InputParameters;       #NSteps,NHoursStage
    @unpack (charge, discharge, rev, cap, soc, soc_quad, deg, deg_stage, WEM_stage, cost_rev, net_revenues_per_stage, x, y, z, w_xx, w_yy, w_zz, w_xy, w_xz, w_zy) = ResultsOpt_3;

    @unpack (vector_prices, vector_stages_index, NSteps_scenario) = Results_ex_post;
    @unpack (WEM_rev_stat_frame, cost_stat_frame, net_rev_stat_frame, deg_stat_frame, Initial_cap_stat_frame, Final_cap_stat_frame) = Results_statistics;
    @unpack (min_SOC, max_SOC, min_P, max_P, Eff_charge, Eff_discharge, max_SOH, min_SOH, Nfull) = Battery ; 

    hour=string(now())
    a=replace(hour,':'=> '-')

    nameF= "$max_SOH SoH_max, $min_SOC SoC_min, mid-cost projections N=7 new, filtrato"
    folder = "$nameF"
    mkdir(folder)
    cd(folder)

    a = pwd()
    cd(a)

    tot_WEM_rev = zeros(NScen)
    tot_rev_costs = zeros(NScen)
    tot_net_rev = zeros(NScen)
    tot_deg = zeros(NScen)
    tot_rev = zeros(NScen)
    initial_cap = zeros(NScen)

    for iScen=1:NScen
        tot_WEM_rev[iScen] = sum(WEM_stage[iScen,:])
        tot_rev_costs[iScen] = sum(cost_rev[iScen,:])
        tot_net_rev[iScen] = sum(net_revenues_per_stage[iScen,:])
        tot_deg[iScen] = sum(deg_stage[iScen,:])
        tot_rev[iScen] = sum(rev[iScen,:])
        initial_cap[iScen] = cap[iScen][2]
    end

    general = DataFrame()
    general[!,"Scenario"] = 1:1:NScen
    general[!,"Initial capacity MWh"] = initial_cap[:]
    general[!,"Revamping MWh"] = tot_rev[:]
    general[!,"Degradation MWh"] = tot_deg[:]
    general[!,"WEM revenues €"] = tot_WEM_rev[:]
    general[!,"Cost revamping"] = tot_rev_costs[:]
    general[!,"Net_Revenues €"] = tot_net_rev[:]

    box_whiskers_WEM = DataFrame()
    box_whiskers_costs = DataFrame()
    box_whiskers_net = DataFrame()

    for iScen=1:NScen
        box_whiskers_WEM[!,"Scen $iScen"] = WEM_stage[iScen,:]
        box_whiskers_costs[!,"Scen $iScen"] = cost_rev[iScen,:]
        box_whiskers_net[!,"Scen $iScen"] = net_revenues_per_stage[iScen,:]
    end

    XLSX.writetable("General results WITH degradation.xlsx", overwrite=true,                                       #$nameFile
    results_scenarios = (collect(DataFrames.eachcol(general)),DataFrames.names(general)),
    WEM_rev = (collect(DataFrames.eachcol(WEM_rev_stat_frame)),DataFrames.names(WEM_rev_stat_frame)),
    Revamping_costs = (collect(DataFrames.eachcol(cost_stat_frame)),DataFrames.names(cost_stat_frame)),
    Net_revenues = (collect(DataFrames.eachcol(net_rev_stat_frame)),DataFrames.names(net_rev_stat_frame)),
    Degradation = (collect(DataFrames.eachcol(deg_stat_frame)),DataFrames.names(deg_stat_frame)),
    Initial_capacity = (collect(DataFrames.eachcol(Initial_cap_stat_frame)),DataFrames.names(Initial_cap_stat_frame)),
    Final_capacity = (collect(DataFrames.eachcol(Final_cap_stat_frame)),DataFrames.names(Final_cap_stat_frame)),

    arbitrage_rev = (collect(DataFrames.eachcol(box_whiskers_WEM)),DataFrames.names(box_whiskers_WEM)),
    revamping_costs = (collect(DataFrames.eachcol(box_whiskers_costs)),DataFrames.names(box_whiskers_costs)),
    net_rev = (collect(DataFrames.eachcol(box_whiskers_net)),DataFrames.names(box_whiskers_net)),
    )

    initial_cap_stage = zeros(NScen, NStages)
    final_cap_stage = zeros(NScen, NStages)

    indice_uno= 0
    indice_due = 0

    # SALVO I SINGOLI SCENARI
    for iScen=1:NScen
            for iStage=1:NStages
                indice_uno = vector_stages_index[iScen,iStage]+2
                indice_due = vector_stages_index[iScen,iStage+1]+1
                initial_cap_stage[iScen, iStage] =cap[iScen][indice_uno]
                final_cap_stage[iScen, iStage] = cap[iScen][indice_due]
            end

            stage_results = DataFrame()

            stage_results[!, "Stages"] = 1:1:NStages
            stage_results[!, "Initial capacity MWh"] = initial_cap_stage[iScen,:]
            stage_results[!, "Final capacity MWh"] = final_cap_stage[iScen,:]
            stage_results[!, "Revamping MWh"] = rev[iScen, :]
            stage_results[!, "Degradation MWh"] = deg_stage[iScen, :]
            stage_results[!, "WEM revenues €"] = WEM_stage[iScen, :]
            stage_results[!, "Revamping costs €"] = cost_rev[iScen,:]
            stage_results[!, "Net revenues €"] = net_revenues_per_stage[iScen,:]

            hourly_results = DataFrame()

            hourly_results[!, "Steps"] = 1:1:NSteps_scenario[iScen]
            hourly_results[!, "Energy pries €/MWh"] = vector_prices[iScen][:]
            hourly_results[!, "Capacity MWh"] = cap[iScen][1:end-1]
            hourly_results[!, "SOC Mwh"] = soc[iScen][1:end-1]
            hourly_results[!, "SOC_quad Mwh"] = soc_quad[iScen][1:end-1]
            hourly_results[!, "Charge MWh"] = charge[iScen][:]
            hourly_results[!, "Discharge MWh"] = discharge[iScen][:]
            hourly_results[!, "X"] = x[iScen][1:end-1]
            hourly_results[!, "Y"] = y[iScen][1:end-1]
            hourly_results[!, "Z"] = z[iScen][1:end-1]
            hourly_results[!, "XX"] = w_xx[iScen][1:end-1]
            hourly_results[!, "YY"] = w_yy[iScen][1:end-1]
            hourly_results[!, "ZZ"] = w_zz[iScen][1:end-1]
            hourly_results[!, "XY"] = w_xy[iScen][1:end-1]
            hourly_results[!, "XZ"] = w_xz[iScen][1:end-1]
            hourly_results[!, "ZY"] = w_zy[iScen][1:end-1]

            XLSX.writetable("$iScen scenario .xlsx", overwrite=true,                                       #$nameFile
            results_stages = (collect(DataFrames.eachcol(stage_results)),DataFrames.names(stage_results)),
            hourly_values = (collect(DataFrames.eachcol(hourly_results)),DataFrames.names(hourly_results)),
            )

    end

    cd(main)             # ritorno nella cartella di salvataggio dati


    return println("Saved data in xlsx")
end


function data_saving_without_deg_3(Results_No_Deg_3::ResultWithoutDeg_3, Results_ex_post::results_ex_post, NSteps, NStages)

    @unpack (revenues_per_stage_no_deg, soc_no_deg, soc_quad_no_deg, charge_no_deg, discharge_no_deg, x_no_deg, y_no_deg,z_no_deg) = Results_No_Deg_3;
    @unpack (deg_singola, deg_stage, costo_stage, net_revenues, vector_prices, vector_stages_index, vector_downtime_stages) = Results_ex_post;

    cart = "Test 3 BIN 27.11"
    mkdir(cart)
    cd(cart)
    main = pwd()

    c = "Risultati senza degradazione, $max_SOH SoH_max, $min_SOC SoC_min, imd cost"
    mkdir(c)
    cd(c)
    
    # SAVE TOTAL RESULTS IN A GENERAL EXCEL
    general_no_deg = DataFrame()

    tot_wem = zeros(NScen);
    tot_costs = zeros(NScen);
    tot_net = zeros(NScen);
    tot_deg = zeros(NScen);
    
    for iScen=1:NScen
        tot_wem[iScen] = sum(revenues_per_stage_no_deg[iScen,:])
        tot_costs[iScen] = sum(costo_stage[iScen,:])
        tot_net[iScen] = sum(net_revenues[iScen,:])
        tot_deg[iScen] = sum(deg_stage[iScen,:])
    end

    general_no_deg[!, "Scenario"] = 1:1:NScen
    general_no_deg[!, "Degradation MWWh"] = tot_deg[:]
    general_no_deg[!, "WEM revenues €"] = tot_wem[:]
    general_no_deg[!, "Revamping costs €"] = tot_costs[:]
    general_no_deg[!, "Net revenues €"] = tot_net[:]

    XLSX.writetable("general results.xlsx", overwrite=true,                                       #$nameFile
    Total = (collect(DataFrames.eachcol(general_no_deg)),DataFrames.names(general_no_deg)),
    )

    #SAVE SINGLE VALUES FOR EACH SCENARIO
    initial_cap = zeros(NScen, NStages);
    final_cap = zeros(NScen,NStages);

    for iScen=1:NScen

        for iStage=1:NStages
            final_cap[iScen, iStage] = min_SOH
            initial_cap[iScen, iStage] = final_cap[iScen, iStage] + deg_stage[iScen, iStage]
        end

        results_stage = DataFrame()
        results_stage[!, "Stage"] = 1:1:NStages
        results_stage[!, "Initial capacity MWh"] = initial_cap[iScen,:]
        results_stage[!, "Final capacity MWh"] = final_cap[iScen,:]
        results_stage[!, "Degradation MWh"] = deg_stage[iScen,:]
        results_stage[!, "WEM revenues €"] = revenues_per_stage_no_deg[iScen,:]
        results_stage[!, "Revamping costs €"] = costo_stage[iScen,:]
        results_stage[!, "Net revenues €"] = net_revenues[iScen,:]
        results_stage[!, "Battery cost€"] = Battery_price_purchase[1:NStages]

        hourly_no_deg = DataFrame()
        hourly_no_deg[!, "Steps"] = 1:1:NSteps
        hourly_no_deg[!, "Power prices €/MWh"] = Pp[:, iScen]
        hourly_no_deg[!, "SOC MW"] = soc_no_deg[iScen, 1:end-1]
        hourly_no_deg[!, "SOC_quad MW"] = soc_quad_no_deg[iScen, 1:end-1]
        hourly_no_deg[!, "Discharged MW"] = discharge_no_deg[iScen, :]
        hourly_no_deg[!, "Charge MWs"] = charge_no_deg[iScen, :]
        hourly_no_deg[!, "Degradation MWh"] = deg_singola[iScen, :]
        hourly_no_deg[!,"X"]= x_no_deg[iScen,1:end-1]
        hourly_no_deg[!,"Y"]= y_no_deg[iScen,1:end-1]
        hourly_no_deg[!,"Z"]= z_no_deg[iScen,1:end-1]

        filtered_pp = DataFrame()
        filtered_pp[!,"Filtered energy prices €/Mwh"] = vector_prices[iScen]

        vector_indexes = DataFrame()
        vector_indexes[!,"Stage indexes"] = vector_stages_index[iScen,:]

        downtime_indexes = DataFrame()
        downtime_indexes[!,"Filtered downtime"] = vector_downtime_stages[iScen,:]

        XLSX.writetable("Scenario $iScen.xlsx", overwrite=true,                                       #$nameFile
        Total = (collect(DataFrames.eachcol(results_stage)),DataFrames.names(results_stage)),
        Hourly_results = (collect(DataFrames.eachcol(hourly_no_deg)),DataFrames.names(hourly_no_deg)),
        Filtered_prices = (collect(DataFrames.eachcol(filtered_pp)),DataFrames.names(filtered_pp)),
        New_stage_vector = (collect(DataFrames.eachcol(vector_indexes)),DataFrames.names(vector_indexes)),
        Downtime_vector = (collect(DataFrames.eachcol(downtime_indexes)),DataFrames.names(downtime_indexes)),
        )
    end

    return (main)

end

function data_saving_4(InputParameters::InputParam, ResultsOpt_4::Results_4, Results_ex_post::results_ex_post, Results_statistics::Data_analysed, Battery::BatteryParam,  NScen)

    @unpack (NYears, NMonths, NStages, Big, NHoursStep, bin) = InputParameters;       #NSteps,NHoursStage
    @unpack (charge, discharge, rev, cap, soc, soc_quad, deg, deg_stage, WEM_stage, cost_rev, net_revenues_per_stage, x, y, z,u, w_xx, w_yy, w_zz, w_uu, w_xy, w_xz, w_zy, w_xu) = ResultsOpt_4;

    @unpack (vector_prices, vector_stages_index, NSteps_scenario) = Results_ex_post;
    @unpack (WEM_rev_stat_frame, cost_stat_frame, net_rev_stat_frame, deg_stat_frame, Initial_cap_stat_frame, Final_cap_stat_frame) = Results_statistics;
    @unpack (min_SOC, max_SOC, min_P, max_P, Eff_charge, Eff_discharge, max_SOH, min_SOH, Nfull) = Battery ; 

    hour=string(now())
    a=replace(hour,':'=> '-')

    nameF= "OCRA filtrato $NScen scenari disc=16, E_Max=$max_SOH "
    folder = "$nameF"
    mkdir(folder)
    cd(folder)

    a = pwd()
    cd(a)

    tot_WEM_rev = zeros(NScen)
    tot_rev_costs = zeros(NScen)
    tot_net_rev = zeros(NScen)
    tot_deg = zeros(NScen)
    tot_rev = zeros(NScen)
    initial_cap = zeros(NScen)

    for iScen=1:NScen
        tot_WEM_rev[iScen] = sum(WEM_stage[iScen,:])
        tot_rev_costs[iScen] = sum(cost_rev[iScen,:])
        tot_net_rev[iScen] = sum(net_revenues_per_stage[iScen,:])
        tot_deg[iScen] = sum(deg_stage[iScen,:])
        tot_rev[iScen] = sum(rev[iScen,:])
        initial_cap[iScen] = cap[iScen][2]
    end

    general = DataFrame()
    general[!,"Scenario"] = 1:1:NScen
    general[!,"Initial capacity MWh"] = initial_cap[:]
    general[!,"Revamping MWh"] = tot_rev[:]
    general[!,"Degradation MWh"] = tot_deg[:]
    general[!,"WEM revenues €"] = tot_WEM_rev[:]
    general[!,"Cost revamping"] = tot_rev_costs[:]
    general[!,"Net_Revenues €"] = tot_net_rev[:]

    box_whiskers_WEM = DataFrame()
    box_whiskers_costs = DataFrame()
    box_whiskers_net = DataFrame()

    for iScen=1:NScen
        box_whiskers_WEM[!,"Scen $iScen"] = WEM_stage[iScen,:]
        box_whiskers_costs[!,"Scen $iScen"] = cost_rev[iScen,:]
        box_whiskers_net[!,"Scen $iScen"] = net_revenues_per_stage[iScen,:]
    end

    XLSX.writetable("General results WITH degradation.xlsx", overwrite=true,                                       #$nameFile
    results_scenarios = (collect(DataFrames.eachcol(general)),DataFrames.names(general)),
    WEM_rev = (collect(DataFrames.eachcol(WEM_rev_stat_frame)),DataFrames.names(WEM_rev_stat_frame)),
    Revamping_costs = (collect(DataFrames.eachcol(cost_stat_frame)),DataFrames.names(cost_stat_frame)),
    Net_revenues = (collect(DataFrames.eachcol(net_rev_stat_frame)),DataFrames.names(net_rev_stat_frame)),
    Degradation = (collect(DataFrames.eachcol(deg_stat_frame)),DataFrames.names(deg_stat_frame)),
    Initial_capacity = (collect(DataFrames.eachcol(Initial_cap_stat_frame)),DataFrames.names(Initial_cap_stat_frame)),
    Final_capacity = (collect(DataFrames.eachcol(Final_cap_stat_frame)),DataFrames.names(Final_cap_stat_frame)),

    arbitrage_rev = (collect(DataFrames.eachcol(box_whiskers_WEM)),DataFrames.names(box_whiskers_WEM)),
    revamping_costs = (collect(DataFrames.eachcol(box_whiskers_costs)),DataFrames.names(box_whiskers_costs)),
    net_rev = (collect(DataFrames.eachcol(box_whiskers_net)),DataFrames.names(box_whiskers_net)),
    )

    initial_cap_stage = zeros(NScen, NStages)
    final_cap_stage = zeros(NScen, NStages)

    indice_uno= 0
    indice_due = 0

    # SALVO I SINGOLI SCENARI
    for iScen=1:NScen
            for iStage=1:NStages
                indice_uno = vector_stages_index[iScen,iStage]+2
                indice_due = vector_stages_index[iScen,iStage+1]+1
                initial_cap_stage[iScen, iStage] =cap[iScen][indice_uno]
                final_cap_stage[iScen, iStage] = cap[iScen][indice_due]
            end

            stage_results = DataFrame()

            stage_results[!, "Stages"] = 1:1:NStages
            stage_results[!, "Initial capacity MWh"] = initial_cap_stage[iScen,:]
            stage_results[!, "Final capacity MWh"] = final_cap_stage[iScen,:]
            stage_results[!, "Revamping MWh"] = rev[iScen, :]
            stage_results[!, "Degradation MWh"] = deg_stage[iScen, :]
            stage_results[!, "WEM revenues €"] = WEM_stage[iScen, :]
            stage_results[!, "Revamping costs €"] = cost_rev[iScen,:]
            stage_results[!, "Net revenues €"] = net_revenues_per_stage[iScen,:]

            hourly_results = DataFrame()

            hourly_results[!, "Steps"] = 1:1:NSteps_scenario[iScen]
            hourly_results[!, "Energy pries €/MWh"] = vector_prices[iScen][:]
            hourly_results[!, "Capacity MWh"] = cap[iScen][1:end-1]
            hourly_results[!, "SOC Mwh"] = soc[iScen][1:end-1]
            hourly_results[!, "SOC_quad Mwh"] = soc_quad[iScen][1:end-1]
            hourly_results[!, "Charge MWh"] = charge[iScen][:]
            hourly_results[!, "Discharge MWh"] = discharge[iScen][:]
            hourly_results[!, "X"] = x[iScen][1:end-1]
            hourly_results[!, "Y"] = y[iScen][1:end-1]
            hourly_results[!, "Z"] = z[iScen][1:end-1]
            hourly_results[!, "U"] = u[iScen][1:end-1]
            hourly_results[!, "XX"] = w_xx[iScen][1:end-1]
            hourly_results[!, "YY"] = w_yy[iScen][1:end-1]
            hourly_results[!, "ZZ"] = w_zz[iScen][1:end-1]
            hourly_results[!, "UU"] = w_uu[iScen][1:end-1]
            hourly_results[!, "XY"] = w_xy[iScen][1:end-1]
            hourly_results[!, "XZ"] = w_xz[iScen][1:end-1]
            hourly_results[!, "ZY"] = w_zy[iScen][1:end-1]

            XLSX.writetable("$iScen scenario .xlsx", overwrite=true,                                       #$nameFile
            results_stages = (collect(DataFrames.eachcol(stage_results)),DataFrames.names(stage_results)),
            hourly_values = (collect(DataFrames.eachcol(hourly_results)),DataFrames.names(hourly_results)),
            )

    end

    cd(main)             # ritorno nella cartella di salvataggio dati


    return println("Saved data in xlsx")
end


function data_saving_without_deg_4(Results_No_Deg_4::ResultWithoutDeg_4, Results_ex_post::results_ex_post, NSteps, NStages)

    @unpack (revenues_per_stage_no_deg, soc_no_deg, soc_quad_no_deg, charge_no_deg, discharge_no_deg, x_no_deg, y_no_deg,z_no_deg, u_no_deg) = Results_No_Deg_4;
    @unpack (deg_singola, deg_stage, costo_stage, net_revenues, vector_prices, vector_stages_index, vector_downtime_stages) = Results_ex_post;

    cart = "Test 4 BIN 27.11"
    mkdir(cart)
    cd(cart)
    main = pwd()

    c = "Risultati senza degradazione"
    mkdir(c)
    cd(c)
    
    # SAVE TOTAL RESULTS IN A GENERAL EXCEL
    general_no_deg = DataFrame()

    tot_wem = zeros(NScen);
    tot_costs = zeros(NScen);
    tot_net = zeros(NScen);
    tot_deg = zeros(NScen);
    
    for iScen=1:NScen
        tot_wem[iScen] = sum(revenues_per_stage_no_deg[iScen,:])
        tot_costs[iScen] = sum(costo_stage[iScen,:])
        tot_net[iScen] = sum(net_revenues[iScen,:])
        tot_deg[iScen] = sum(deg_stage[iScen,:])
    end

    general_no_deg[!, "Scenario"] = 1:1:NScen
    general_no_deg[!, "Degradation MWWh"] = tot_deg[:]
    general_no_deg[!, "WEM revenues €"] = tot_wem[:]
    general_no_deg[!, "Revamping costs €"] = tot_costs[:]
    general_no_deg[!, "Net revenues €"] = tot_net[:]

    XLSX.writetable("general results.xlsx", overwrite=true,                                       #$nameFile
    Total = (collect(DataFrames.eachcol(general_no_deg)),DataFrames.names(general_no_deg)),
    )

    #SAVE SINGLE VALUES FOR EACH SCENARIO
    initial_cap = zeros(NScen, NStages);
    final_cap = zeros(NScen,NStages);

    for iScen=1:NScen

        for iStage=1:NStages
            final_cap[iScen, iStage] = min_SOH
            initial_cap[iScen, iStage] = final_cap[iScen, iStage] + deg_stage[iScen, iStage]
        end

        results_stage = DataFrame()
        results_stage[!, "Stage"] = 1:1:NStages
        results_stage[!, "Initial capacity MWh"] = initial_cap[iScen,:]
        results_stage[!, "Final capacity MWh"] = final_cap[iScen,:]
        results_stage[!, "Degradation MWh"] = deg_stage[iScen,:]
        results_stage[!, "WEM revenues €"] = revenues_per_stage_no_deg[iScen,:]
        results_stage[!, "Revamping costs €"] = costo_stage[iScen,:]
        results_stage[!, "Net revenues €"] = net_revenues[iScen,:]
        results_stage[!, "Battery cost€"] = Battery_price_purchase[1:NStages]

        hourly_no_deg = DataFrame()
        hourly_no_deg[!, "Steps"] = 1:1:NSteps
        hourly_no_deg[!, "Power prices €/MWh"] = Pp[:, iScen]
        hourly_no_deg[!, "SOC MW"] = soc_no_deg[iScen, 1:end-1]
        hourly_no_deg[!, "SOC_quad MW"] = soc_quad_no_deg[iScen, 1:end-1]
        hourly_no_deg[!, "Discharged MW"] = discharge_no_deg[iScen, :]
        hourly_no_deg[!, "Charge MWs"] = charge_no_deg[iScen, :]
        hourly_no_deg[!, "Degradation MWh"] = deg_singola[iScen, :]
        hourly_no_deg[!,"X"]= x_no_deg[iScen,1:end-1]
        hourly_no_deg[!,"Y"]= y_no_deg[iScen,1:end-1]
        hourly_no_deg[!,"Z"]= z_no_deg[iScen,1:end-1]
        hourly_no_deg[!,"U"]= u_no_deg[iScen,1:end-1]

        filtered_pp = DataFrame()
        filtered_pp[!,"Filtered energy prices €/Mwh"] = vector_prices[iScen]

        vector_indexes = DataFrame()
        vector_indexes[!,"Stage indexes"] = vector_stages_index[iScen,:]

        downtime_indexes = DataFrame()
        downtime_indexes[!,"Filtered downtime"] = vector_downtime_stages[iScen,:]

        XLSX.writetable("Scenario $iScen.xlsx", overwrite=true,                                       #$nameFile
        Total = (collect(DataFrames.eachcol(results_stage)),DataFrames.names(results_stage)),
        Hourly_results = (collect(DataFrames.eachcol(hourly_no_deg)),DataFrames.names(hourly_no_deg)),
        Filtered_prices = (collect(DataFrames.eachcol(filtered_pp)),DataFrames.names(filtered_pp)),
        New_stage_vector = (collect(DataFrames.eachcol(vector_indexes)),DataFrames.names(vector_indexes)),
        Downtime_vector = (collect(DataFrames.eachcol(downtime_indexes)),DataFrames.names(downtime_indexes)),
        )
    end

    return (main)

end


function data_saving_5(InputParameters::InputParam, ResultsOpt_5::Results_5, Results_ex_post::results_ex_post, Results_statistics::Data_analysed, Battery::BatteryParam,  NScen)

    @unpack (NYears, NMonths, NStages, Big, NHoursStep, bin) = InputParameters;       #NSteps,NHoursStage
    @unpack (charge, discharge, rev, cap, soc, soc_quad, deg, deg_stage, WEM_stage, cost_rev, net_revenues_per_stage, x, y, z,u,t, w_xx, w_yy, w_zz, w_uu,w_tt, w_xy, w_xz, w_zy, w_xu) = ResultsOpt_5;

    @unpack (vector_prices, vector_stages_index, NSteps_scenario) = Results_ex_post;
    @unpack (WEM_rev_stat_frame, cost_stat_frame, net_rev_stat_frame, deg_stat_frame, Initial_cap_stat_frame, Final_cap_stat_frame) = Results_statistics;
    @unpack (min_SOC, max_SOC, min_P, max_P, Eff_charge, Eff_discharge, max_SOH, min_SOH, Nfull) = Battery ; 

    hour=string(now())
    a=replace(hour,':'=> '-')

    nameF= "OCRA filtrato $NScen scenari disc=16, E_Max=$max_SOH "
    folder = "$nameF"
    mkdir(folder)
    cd(folder)

    a = pwd()
    cd(a)

    tot_WEM_rev = zeros(NScen)
    tot_rev_costs = zeros(NScen)
    tot_net_rev = zeros(NScen)
    tot_deg = zeros(NScen)
    tot_rev = zeros(NScen)
    initial_cap = zeros(NScen)

    for iScen=1:NScen
        tot_WEM_rev[iScen] = sum(WEM_stage[iScen,:])
        tot_rev_costs[iScen] = sum(cost_rev[iScen,:])
        tot_net_rev[iScen] = sum(net_revenues_per_stage[iScen,:])
        tot_deg[iScen] = sum(deg_stage[iScen,:])
        tot_rev[iScen] = sum(rev[iScen,:])
        initial_cap[iScen] = cap[iScen][2]
    end

    general = DataFrame()
    general[!,"Scenario"] = 1:1:NScen
    general[!,"Initial capacity MWh"] = initial_cap[:]
    general[!,"Revamping MWh"] = tot_rev[:]
    general[!,"Degradation MWh"] = tot_deg[:]
    general[!,"WEM revenues €"] = tot_WEM_rev[:]
    general[!,"Cost revamping"] = tot_rev_costs[:]
    general[!,"Net_Revenues €"] = tot_net_rev[:]

    box_whiskers_WEM = DataFrame()
    box_whiskers_costs = DataFrame()
    box_whiskers_net = DataFrame()

    for iScen=1:NScen
        box_whiskers_WEM[!,"Scen $iScen"] = WEM_stage[iScen,:]
        box_whiskers_costs[!,"Scen $iScen"] = cost_rev[iScen,:]
        box_whiskers_net[!,"Scen $iScen"] = net_revenues_per_stage[iScen,:]
    end

    XLSX.writetable("General results WITH degradation.xlsx", overwrite=true,                                       #$nameFile
    results_scenarios = (collect(DataFrames.eachcol(general)),DataFrames.names(general)),
    WEM_rev = (collect(DataFrames.eachcol(WEM_rev_stat_frame)),DataFrames.names(WEM_rev_stat_frame)),
    Revamping_costs = (collect(DataFrames.eachcol(cost_stat_frame)),DataFrames.names(cost_stat_frame)),
    Net_revenues = (collect(DataFrames.eachcol(net_rev_stat_frame)),DataFrames.names(net_rev_stat_frame)),
    Degradation = (collect(DataFrames.eachcol(deg_stat_frame)),DataFrames.names(deg_stat_frame)),
    Initial_capacity = (collect(DataFrames.eachcol(Initial_cap_stat_frame)),DataFrames.names(Initial_cap_stat_frame)),
    Final_capacity = (collect(DataFrames.eachcol(Final_cap_stat_frame)),DataFrames.names(Final_cap_stat_frame)),

    arbitrage_rev = (collect(DataFrames.eachcol(box_whiskers_WEM)),DataFrames.names(box_whiskers_WEM)),
    revamping_costs = (collect(DataFrames.eachcol(box_whiskers_costs)),DataFrames.names(box_whiskers_costs)),
    net_rev = (collect(DataFrames.eachcol(box_whiskers_net)),DataFrames.names(box_whiskers_net)),
    )

    initial_cap_stage = zeros(NScen, NStages)
    final_cap_stage = zeros(NScen, NStages)

    indice_uno= 0
    indice_due = 0

    # SALVO I SINGOLI SCENARI
    for iScen=1:NScen
            for iStage=1:NStages
                indice_uno = vector_stages_index[iScen,iStage]+2
                indice_due = vector_stages_index[iScen,iStage+1]+1
                initial_cap_stage[iScen, iStage] =cap[iScen][indice_uno]
                final_cap_stage[iScen, iStage] = cap[iScen][indice_due]
            end

            stage_results = DataFrame()

            stage_results[!, "Stages"] = 1:1:NStages
            stage_results[!, "Initial capacity MWh"] = initial_cap_stage[iScen,:]
            stage_results[!, "Final capacity MWh"] = final_cap_stage[iScen,:]
            stage_results[!, "Revamping MWh"] = rev[iScen, :]
            stage_results[!, "Degradation MWh"] = deg_stage[iScen, :]
            stage_results[!, "WEM revenues €"] = WEM_stage[iScen, :]
            stage_results[!, "Revamping costs €"] = cost_rev[iScen,:]
            stage_results[!, "Net revenues €"] = net_revenues_per_stage[iScen,:]

            hourly_results = DataFrame()

            hourly_results[!, "Steps"] = 1:1:NSteps_scenario[iScen]
            hourly_results[!, "Energy pries €/MWh"] = vector_prices[iScen][:]
            hourly_results[!, "Capacity MWh"] = cap[iScen][1:end-1]
            hourly_results[!, "SOC Mwh"] = soc[iScen][1:end-1]
            hourly_results[!, "SOC_quad Mwh"] = soc_quad[iScen][1:end-1]
            hourly_results[!, "Charge MWh"] = charge[iScen][:]
            hourly_results[!, "Discharge MWh"] = discharge[iScen][:]
            hourly_results[!, "X"] = x[iScen][1:end-1]
            hourly_results[!, "Y"] = y[iScen][1:end-1]
            hourly_results[!, "Z"] = z[iScen][1:end-1]
            hourly_results[!, "U"] = u[iScen][1:end-1]
            hourly_results[!, "T"] = t[iScen][1:end-1]
            hourly_results[!, "XX"] = w_xx[iScen][1:end-1]
            hourly_results[!, "YY"] = w_yy[iScen][1:end-1]
            hourly_results[!, "ZZ"] = w_zz[iScen][1:end-1]
            hourly_results[!, "UU"] = w_uu[iScen][1:end-1]
            hourly_results[!, "TT"] = w_tt[iScen][1:end-1]
            hourly_results[!, "XY"] = w_xy[iScen][1:end-1]
            hourly_results[!, "XZ"] = w_xz[iScen][1:end-1]
            hourly_results[!, "ZY"] = w_zy[iScen][1:end-1]

            XLSX.writetable("$iScen scenario .xlsx", overwrite=true,                                       #$nameFile
            results_stages = (collect(DataFrames.eachcol(stage_results)),DataFrames.names(stage_results)),
            hourly_values = (collect(DataFrames.eachcol(hourly_results)),DataFrames.names(hourly_results)),
            )

    end

    cd(main)             # ritorno nella cartella di salvataggio dati


    return println("Saved data in xlsx")
end


function data_saving_without_deg_5(Results_No_Deg_5::ResultWithoutDeg_5, Results_ex_post::results_ex_post, NSteps, NStages)

    @unpack (revenues_per_stage_no_deg, soc_no_deg, soc_quad_no_deg, charge_no_deg, discharge_no_deg, x_no_deg, y_no_deg,z_no_deg, u_no_deg, t_no_deg) = Results_No_Deg_5;
    @unpack (deg_singola, deg_stage, costo_stage, net_revenues, vector_prices, vector_stages_index, vector_downtime_stages) = Results_ex_post;

    cart = "Test 5 BIN 27.11"
    mkdir(cart)
    cd(cart)
    main = pwd()

    c = "Risultati senza degradazione"
    mkdir(c)
    cd(c)
    
    # SAVE TOTAL RESULTS IN A GENERAL EXCEL
    general_no_deg = DataFrame()

    tot_wem = zeros(NScen);
    tot_costs = zeros(NScen);
    tot_net = zeros(NScen);
    tot_deg = zeros(NScen);
    
    for iScen=1:NScen
        tot_wem[iScen] = sum(revenues_per_stage_no_deg[iScen,:])
        tot_costs[iScen] = sum(costo_stage[iScen,:])
        tot_net[iScen] = sum(net_revenues[iScen,:])
        tot_deg[iScen] = sum(deg_stage[iScen,:])
    end

    general_no_deg[!, "Scenario"] = 1:1:NScen
    general_no_deg[!, "Degradation MWWh"] = tot_deg[:]
    general_no_deg[!, "WEM revenues €"] = tot_wem[:]
    general_no_deg[!, "Revamping costs €"] = tot_costs[:]
    general_no_deg[!, "Net revenues €"] = tot_net[:]

    XLSX.writetable("general results.xlsx", overwrite=true,                                       #$nameFile
    Total = (collect(DataFrames.eachcol(general_no_deg)),DataFrames.names(general_no_deg)),
    )

    #SAVE SINGLE VALUES FOR EACH SCENARIO
    initial_cap = zeros(NScen, NStages);
    final_cap = zeros(NScen,NStages);

    for iScen=1:NScen

        for iStage=1:NStages
            final_cap[iScen, iStage] = min_SOH
            initial_cap[iScen, iStage] = final_cap[iScen, iStage] + deg_stage[iScen, iStage]
        end

        results_stage = DataFrame()
        results_stage[!, "Stage"] = 1:1:NStages
        results_stage[!, "Initial capacity MWh"] = initial_cap[iScen,:]
        results_stage[!, "Final capacity MWh"] = final_cap[iScen,:]
        results_stage[!, "Degradation MWh"] = deg_stage[iScen,:]
        results_stage[!, "WEM revenues €"] = revenues_per_stage_no_deg[iScen,:]
        results_stage[!, "Revamping costs €"] = costo_stage[iScen,:]
        results_stage[!, "Net revenues €"] = net_revenues[iScen,:]
        results_stage[!, "Battery cost€"] = Battery_price_purchase[1:NStages]

        hourly_no_deg = DataFrame()
        hourly_no_deg[!, "Steps"] = 1:1:NSteps
        hourly_no_deg[!, "Power prices €/MWh"] = Pp[:, iScen]
        hourly_no_deg[!, "SOC MW"] = soc_no_deg[iScen, 1:end-1]
        hourly_no_deg[!, "SOC_quad MW"] = soc_quad_no_deg[iScen, 1:end-1]
        hourly_no_deg[!, "Discharged MW"] = discharge_no_deg[iScen, :]
        hourly_no_deg[!, "Charge MWs"] = charge_no_deg[iScen, :]
        hourly_no_deg[!, "Degradation MWh"] = deg_singola[iScen, :]
        hourly_no_deg[!,"X"]= x_no_deg[iScen,1:end-1]
        hourly_no_deg[!,"Y"]= y_no_deg[iScen,1:end-1]
        hourly_no_deg[!,"Z"]= z_no_deg[iScen,1:end-1]
        hourly_no_deg[!,"U"]= u_no_deg[iScen,1:end-1]
        hourly_no_deg[!,"T"]= t_no_deg[iScen,1:end-1]

        filtered_pp = DataFrame()
        filtered_pp[!,"Filtered energy prices €/Mwh"] = vector_prices[iScen]

        vector_indexes = DataFrame()
        vector_indexes[!,"Stage indexes"] = vector_stages_index[iScen,:]

        downtime_indexes = DataFrame()
        downtime_indexes[!,"Filtered downtime"] = vector_downtime_stages[iScen,:]

        XLSX.writetable("Scenario $iScen.xlsx", overwrite=true,                                       #$nameFile
        Total = (collect(DataFrames.eachcol(results_stage)),DataFrames.names(results_stage)),
        Hourly_results = (collect(DataFrames.eachcol(hourly_no_deg)),DataFrames.names(hourly_no_deg)),
        Filtered_prices = (collect(DataFrames.eachcol(filtered_pp)),DataFrames.names(filtered_pp)),
        New_stage_vector = (collect(DataFrames.eachcol(vector_indexes)),DataFrames.names(vector_indexes)),
        Downtime_vector = (collect(DataFrames.eachcol(downtime_indexes)),DataFrames.names(downtime_indexes)),
        )
    end

    return (main)

end
