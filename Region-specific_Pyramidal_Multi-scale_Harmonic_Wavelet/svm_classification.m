function [Beta, p_value, SVM_para, node_index_index] = svm_classification(Phi_set, SUVR_name_set, Classification_Group, hub_nodes, scale, wavelet_num, p, paran, All_SUVR)
    %% SVM classification
    performance = ["Accuracy", "Specificity", "Sensitivity", "Fscore"];
    value_static = ["mean", "std"];
    row_num = size(Phi_set, 2) * size(SUVR_name_set, 2);
    col_num = size(Classification_Group, 2);
    for i = 1:size(performance, 2)
        eval(sprintf('%s_fold = zeros(row_num * paran, col_num)', performance(i)));
        for j = 1:size(value_static, 2)
            eval(sprintf('%s_fold_%s = zeros(row_num, col_num)', performance(i), value_static(j)));
        end
    end
    Beta = cell(row_num, col_num);
    p_value = zeros(row_num * size(performance, 2), col_num);

    wb_svm = parwaitbar(size(Classification_Group, 2) * size(SUVR_name_set, 2) * size(Phi_set, 2), 'WaitMessage', 'SVM classifying...', 'FinalMessage', 'Done!');
    for i = 1:size(Classification_Group, 2)
        Group_name_split_index = strfind(Classification_Group(i), '_');
        a_index = Group_name_split_index(1) - 1;
        b_index = Group_name_split_index(2) + 1;
        Group_a_name = Classification_Group(i)(1:a_index);
        Group_b_name = Classification_Group(i)(b_index:end);
        for j = 1:size(SUVR_name_set, 2)
            for k = 1:size(Phi_set, 2)
                index_tmp = k + (j - 1) * size(Phi_set, 2);

                if isempty(Phi_set(k))
                    group_a_signdata = All_SUVR.(SUVR_name_set(j)).(Group_a_name);
                    group_b_signdata = All_SUVR.(SUVR_name_set(j)).(Group_b_name);
                else
                    group_a_signdata = transpose(Phi_set{k}) * All_SUVR.(SUVR_name_set(j)).(Group_a_name);
                    group_b_signdata = transpose(Phi_set{k}) * All_SUVR.(SUVR_name_set(j)).(Group_b_name);
                end
                row_start = (j - 1) * paran * size(Phi_set, 2) + (k - 1) * paran + 1;
                row_end = (j - 1) * paran * size(Phi_set, 2) + (k - 1) * paran + 10;
                [Accuracy_fold(row_start:row_end, i), Sensitivity_fold(row_start:row_end, i), Specificity_fold(row_start:row_end, i), Fscore_fold(row_start:row_end, i), Beta{index_tmp, i}] = SVM_signaldata(group_a_signdata, group_b_signdata);

                for l = 1:size(performance, 2)
                    array_fold = eval(sprintf('%s_fold(row_start:row_end, i)', performance(l)));
                    [mean_value, std_value] = Mean_Std_calculation(array_fold);
                    for m = 1:size(value_static, 2)
                        eval(sprintf('%s_fold_%s(index_tmp, i) = %s_value', performance(l), value_static(m), value_static(m)));
                    end
                    t_test_vector_a = eval(sprintf('%s_fold(row_start:row_end, i)', performance(l)));
                    row_start_fixed_phi = (j - 1) * paran * size(Phi_set, 2) + 1;
                    row_end_fixed_phi = (j - 1) * paran * size(Phi_set, 2) + 10;
                    t_test_vector_b = eval(sprintf('%s_fold(row_start_fixed_phi:row_end_fixed_phi, i)', performance(l)));
                    [p_value(index_tmp, i), ~, ~] = permutationTest(t_test_vector_a, t_test_vector_b, 10000);
                end

                wb_svm.progress();
            end
        end
    end

    SVM_para = cell(nhub, size(Classification_Group, 2));
    for j = 1:size(Classification_Group, 2)
        for i = 1:nhub
            SVM_para{i, j} = Beta{1, j}(1 + 18 * (i - 1):18 * i, 1:10);
            writematrix(SVM_para{i, j}.', '../results/Region-specific_Pyramidal_Multi-scale_Harmonic_Wavelet/SVM_weight.xls', 'WriteMode', 'append');
        end
    end

    save('../results/Region-specific_Pyramidal_Multi-scale_Harmonic_Wavelet/SVM_para.mat', 'SVM_para');

    for i = 1:size(performance, 2)
        for j = 1:size(value_static, 2)
            save_static_preformance_var_name = eval(sprintf('%s_fold_%s', performance(i), value_static(j)));
            save_static_preformance_filename = append('../results/Region-specific_Pyramidal_Multi-scale_Harmonic_Wavelet/', performance(i), '_', value_static(j), '.xls');
            writematrix(save_static_preformance_var_name, save_static_preformance_filename, 'WriteMode', 'append');
        end
    end

    select_scale = cell(nhub, size(Classification_Group, 2));
    select_m_scale = zeros(scale * nhub, size(Classification_Group));
    select_wavelets_index = zeros(scale * nhub, size(Classification_Group));

    for j = 1:size(Classification_Group, 2)
        for i = 1:nhub
            select_scale{i, j} = mean(abs(SVM_para{i, j}), 2);
            for k = 1:scale
                [select_m_scale((k - 1) * nhub + i, j), select_wavelets_index((k - 1) * nhub + i, j)] = max(select_scale{i, j}(wavelet_num * (k - 1) + 1:wavelet_num * k, 1));
            end
            writematrix(select_scale{i, j}.', '../results/Region-specific_Pyramidal_Multi-scale_Harmonic_Wavelet/select_scale.xls', 'WriteMode', 'append');
        end
    end

    for j = 1:size(Classification_Group, 2)
        for k = 1:scale
            [~, node_index_by_scale] = maxk(select_m_scale((k - 1) * nhub + 1:k * nhub, j), 6);
        end
        [~, node_index_ignore_scale] = maxk(select_m_scale, 18);
    end
    node_index_index = hub_nodes(node_index_by_scale);
    writematrix(node_index_index, '../results/Region-specific_Pyramidal_Multi-scale_Harmonic_Wavelet/node_index_index.xls', 'WriteMode', 'append');
end