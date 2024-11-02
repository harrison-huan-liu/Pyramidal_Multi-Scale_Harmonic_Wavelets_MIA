function [Classification_Group, SUVR_name_set, Phi_set, All_SUVR]=load_suvr()
    %% SUVR
    demographic_header_table = ["AD", "CN", "EMCI", "LMCI", "SMC"];
    Classification_Group = ["AD_vs_CN", "EMCI_vs_CN", "LMCI_vs_EMCI", "LMCI_vs_CN"];
    SUVR_name_set = ["Amyloid", "Tau", "FDG"];
    Phi_set = ["multiPhi", "singlePhi", "Phi_ave", ""];
    All_SUVR = struct();

    wb_preprocess = parwaitbar(size(SUVR_name_set,2),'WaitMessage','Load SUVR data...','FinalMessage','Done!');
    for i = 1:size(SUVR_name_set,2)
        preprocess_function_str = append('Proprocess_original_',SUVR_name_set(i),'_data');
        preprocess_function = str2func(preprocess_function_str);
        SUVR_datafolder = append('../Data/signal_Data/', SUVR_name_set(i), '_SUVR.xlsx');
        [AD_SUVR,CN_SUVR,EMCI_SUVR,LMCI_SUVR,SMC_SUVR]=preprocess_function(SUVR_datafolder);

        All_SUVR.(SUVR_name_set(i)).AD = AD_SUVR;
        All_SUVR.(SUVR_name_set(i)).CN = CN_SUVR;
        All_SUVR.(SUVR_name_set(i)).EMCI = EMCI_SUVR;
        All_SUVR.(SUVR_name_set(i)).LMCI = LMCI_SUVR;
        All_SUVR.(SUVR_name_set(i)).SMC = SMC_SUVR;

        % Age_Gender_Demography_Statistic
        [Female_table,Male_table,demographic_table] = Age_Gender_Demography_Statistic(AD_SUVR,CN_SUVR,EMCI_SUVR,LMCI_SUVR,SMC_SUVR,demographic_header_table,SUVR_name_set(i));

        wb_preprocess.progress();
    end
end