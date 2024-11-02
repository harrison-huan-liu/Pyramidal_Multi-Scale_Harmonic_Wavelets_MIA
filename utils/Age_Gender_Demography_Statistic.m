function [Female_table,Male_table,demographic_table] = Age_Gender_Demography_Statistic(AD_SUVR,CN_SUVR,EMCI_SUVR,LMCI_SUVR,SMC_SUVR,demographic_header_table,SUVR_name)
    %% output: Female_table,Male_table,demographic_table,
    % AD_female_age,CN_female_age,EMCI_female_age,LMCI_female_age,SMC_female_age
    % AD_male_age,CN_male_age,EMCI_male_age,LMCI_male_age,SMC_male_age
    Age_female_value_table = zeros(1, size(demographic_header_table, 2));
    Age_male_value_table = zeros(1, size(demographic_header_table, 2));
    Num_female_value_table = zeros(1, size(demographic_header_table, 2));
    Num_male_value_table = zeros(1, size(demographic_header_table, 2));

    for i = 1:size(demographic_header_table, 2)
        SUVR = eval([demographic_header_table{i}, '_SUVR']);
        female_age = [];
        male_age = [];
        
        for j = 1:size(SUVR.SubjectID, 1)
            if strcmp(SUVR.Gender{j, 1}, 'Female')
                female_age(end + 1) = SUVR.Age{j, 1};
            else
                male_age(end + 1) = SUVR.Age{j, 1};
            end
        end
        
        Age_female_value_table(1, i) = mean(female_age);
        Age_male_value_table(1, i) = mean(male_age);
        Num_female_value_table(1, i) = numel(female_age);
        Num_male_value_table(1, i) = numel(male_age);
    end

    Age_female_value_table = num2str(Age_female_value_table, '%.3f');
    Age_male_value_table = num2str(Age_male_value_table, '%.3f');
    Num_female_value_table = num2str(Num_female_value_table, '%u');
    Num_male_value_table = num2str(Num_male_value_table, '%u');

    Female_table = [demographic_header_table; Age_female_value_table; Num_female_value_table];
    Male_table = [demographic_header_table; Age_male_value_table; Num_male_value_table];

    demographic_value_table = zeros(1, size(demographic_header_table, 2));
    for i = 1:size(demographic_header_table, 2)
        SUVR = eval([demographic_header_table{i}, '_SUVR']);
        demographic_value_table(1, i) = size(SUVR.(SUVR_name), 2);
    end

    demographic_value_table = num2str(demographic_value_table, '%u');
    demographic_table = [demographic_header_table; demographic_value_table];
end