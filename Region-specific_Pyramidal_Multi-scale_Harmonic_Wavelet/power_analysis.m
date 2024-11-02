function [energy, consider_group]=power_analysis(SUVR_name_set, multiPhi, nhub, scale, wavelet_num, All_SUVR)
    %% Power of Wavelets and Biomarkers
    consider_group = ["CN", "EMCI", "LMCI"];
    energy = struct();
    for i = 1:size(SUVR_name_set,2)
        consider_group_num = 0;
        All_signdata = [];
        for j = 1:size(consider_group,2)
            eval(sprintf('SUVR = All_SUVR.%s.%s;',SUVR_name_set(i),consider_group(j)));
            consider_group_num = consider_group_num + size(SUVR,2);
            All_signdata = [All_signdata, SUVR];
        end
        Energy_signdata = zeros(consider_group_num,1);
        Energy_signdata_square = zeros(consider_group_num,1);
        Energy_signdata_each = zeros(consider_group_num,nhub*scale*wavelet_num);

        for i = 1:consider_group_num
            for j = 1:nhub*scale*wavelet_num
                Energy_signdata(i,1) = Energy_signdata(i,1) + dot(All_signdata(1:148,i),multiPhi(1:148,j));
                Energy_signdata_square(i,1) = Energy_signdata(i,1)*Energy_signdata(i,1);
                Energy_signdata_each(i,j) = dot(All_signdata(1:148,i),multiPhi(1:148,j));
            end
        end
        colmin = min(Energy_signdata_each);
        colmax = max(Energy_signdata_each);
        Energy_signdata_norm = rescale(Energy_signdata_each,'InputMin',colmin,'InputMax',colmax);
        Energy_signdata_norm_num = sum(Energy_signdata_norm,2);
        Energy_signdata_norm_num_square = sum(Energy_signdata_norm.*Energy_signdata_norm,2);
        writematrix(Energy_signdata_norm_num_square,['../results/Region-specific_Pyramidal_Multi-scale_Harmonic_Wavelet/', SUVR_name_set(i), '_energy_norm_value.xls'],'WriteMode','append');

        energy.(SUVR_name_set(i)).Energy_signdata = Energy_signdata;
        energy.(SUVR_name_set(i)).Energy_signdata_square = Energy_signdata_square;
        energy.(SUVR_name_set(i)).Energy_signdata_each = Energy_signdata_each;
        energy.(SUVR_name_set(i)).consider_group_num = consider_group_num;
        energy.(SUVR_name_set(i)).Energy_signdata_norm = Energy_signdata_norm;
        energy.(SUVR_name_set(i)).Energy_signdata_norm_num = Energy_signdata_norm_num;
        energy.(SUVR_name_set(i)).Energy_signdata_norm_num_square = Energy_signdata_norm_num_square;
    end
end