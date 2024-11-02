function [Energy_signdata_norm_num,Energy_signdata_norm_num_square] = energy_cal(CL_num,LC_signdata,multiPhi,dir_name)
Energy_signdata = zeros(CL_num,1);
Energy_signdata_square = zeros(CL_num,1);
Energy_signdata_each = zeros(CL_num,9);
for i = 1:CL_num
    for j = 1:9
        Energy_signdata(i,1) = Energy_signdata(i,1) + dot(LC_signdata(1:148,i),multiPhi(1:148,j));
        Energy_signdata_square(i,1) = Energy_signdata(i,1)*Energy_signdata(i,1);
        Energy_signdata_each(i,j) = dot(LC_signdata(1:148,i),multiPhi(1:148,j));
    end
end

colmin = min(Energy_signdata_each);
colmax = max(Energy_signdata_each);
Energy_signdata_norm = rescale(Energy_signdata_each,'InputMin',colmin,'InputMax',colmax);

Energy_signdata_norm_num = sum(Energy_signdata_norm,2);
Energy_signdata_norm_num_square = sum(Energy_signdata_each.*Energy_signdata_each,2);

filename1 = ['results_',dir_name,'\Energy_signdata_norm_num.mat'];
save(filename1,'Energy_signdata_norm_num');

filename2 = ['results_',dir_name,'\Energy_signdata_norm_num_square.mat'];
save(filename2,'Energy_signdata_norm_num_square');
end