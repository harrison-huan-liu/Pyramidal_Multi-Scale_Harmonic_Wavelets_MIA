function [performance_fold,performance,performance_std,p_fold,SVM_weight_multi_max,SVM_weight_multi_max_ce,SVM_weight_multi_max_el] = Classification_svm(AD_data,LMCI_data,EMCI_data,CN_data,multiPhi,singlePhi,Phi_ave,dir_name)
%% accuracy,accuracy_std,specificity,specificity_std,sensitivity,sensitivity_std,fscore,fscore_std
paran = 100;
AC_num = size(AD_data,2)+size(CN_data,2);
LC_num = size(LMCI_data,2)+size(CN_data,2);
CE_num = size(CN_data,2)+size(EMCI_data,2);
EL_num = size(EMCI_data,2)+size(LMCI_data,2);
SVM_weight_multi = cell(paran,4);

accuracy_ac = zeros(paran,4);
accuracy_ce = zeros(paran,4);
accuracy_el = zeros(paran,4);
accuracy_lc = zeros(paran,4);

specificity_ac = zeros(paran,4);
specificity_ce = zeros(paran,4);
specificity_el = zeros(paran,4);
specificity_lc = zeros(paran,4);

sensitivity_ac = zeros(paran,4);
sensitivity_ce = zeros(paran,4);
sensitivity_el = zeros(paran,4);
sensitivity_lc = zeros(paran,4);

fscore_ac = zeros(paran,4);
fscore_ce = zeros(paran,4);
fscore_el = zeros(paran,4);
fscore_lc = zeros(paran,4);

accuracy = zeros(4,4);
accuracy_std = zeros(4,4);
specificity = zeros(4,4);
specificity_std = zeros(4,4);
sensitivity = zeros(4,4);
sensitivity_std = zeros(4,4);
fscore = zeros(4,4);
fscore_std = zeros(4,4);

% performance_fold = cell(16*paran,4);
% performance = zeros(16,4);
% performance_std = zeros(16,4);
p_fold = zeros(16,4);

for j = 1:paran
%     Classification_original
    [accuracy_ac(j,1),sensitivity_ac(j,1),specificity_ac(j,1),fscore_ac(j,1),SVM_weight_multi{j,1}] = Classification_original(multiPhi.'*CN_data,multiPhi.'*AD_data,AC_num);
    [accuracy_ce(j,1),sensitivity_ce(j,1),specificity_ce(j,1),fscore_ce(j,1),SVM_weight_multi{j,2}] = Classification_original(multiPhi.'*CN_data,multiPhi.'*EMCI_data,CE_num);
%     [accuracy_el(j,1),sensitivity_el(j,1),specificity_el(j,1),fscore_el(j,1),SVM_weight_multi{j,3}] = Classification_original(multiPhi.'*EMCI_data,multiPhi.'*LMCI_data,EL_num);
%     [accuracy_lc(j,1),sensitivity_lc(j,1),specificity_lc(j,1),fscore_lc(j,1),SVM_weight_multi{j,4}] = Classification_original(multiPhi.'*CN_data,multiPhi.'*LMCI_data,LC_num);
% 
%     [accuracy_ac(j,2),sensitivity_ac(j,2),specificity_ac(j,2),fscore_ac(j,2),~] = Classification_original(singlePhi.'*CN_data,singlePhi.'*AD_data,AC_num);
%     [accuracy_ce(j,2),sensitivity_ce(j,2),specificity_ce(j,2),fscore_ce(j,2),~] = Classification_original(singlePhi.'*CN_data,singlePhi.'*EMCI_data,CE_num);
%     [accuracy_el(j,2),sensitivity_el(j,2),specificity_el(j,2),fscore_el(j,2),~] = Classification_original(singlePhi.'*EMCI_data,singlePhi.'*LMCI_data,EL_num);
%     [accuracy_lc(j,2),sensitivity_lc(j,2),specificity_lc(j,2),fscore_lc(j,2),~] = Classification_original(singlePhi.'*CN_data,singlePhi.'*LMCI_data,LC_num);
% 
%     [accuracy_ac(j,3),sensitivity_ac(j,3),specificity_ac(j,3),fscore_ac(j,3),~] = Classification_original(Phi_ave.'*CN_data,Phi_ave.'*AD_data,AC_num);
%     [accuracy_ce(j,3),sensitivity_ce(j,3),specificity_ce(j,3),fscore_ce(j,3),~] = Classification_original(Phi_ave.'*CN_data,Phi_ave.'*EMCI_data,CE_num);
%     [accuracy_el(j,3),sensitivity_el(j,3),specificity_el(j,3),fscore_el(j,3),~] = Classification_original(Phi_ave.'*EMCI_data,Phi_ave.'*LMCI_data,EL_num);
%     [accuracy_lc(j,3),sensitivity_lc(j,3),specificity_lc(j,3),fscore_lc(j,3),~] = Classification_original(Phi_ave.'*CN_data,Phi_ave.'*LMCI_data,LC_num);
% 
%     [accuracy_ac(j,4),sensitivity_ac(j,4),specificity_ac(j,4),fscore_ac(j,4),~] = Classification_original(CN_data,AD_data,AC_num);
%     [accuracy_ce(j,4),sensitivity_ce(j,4),specificity_ce(j,4),fscore_ce(j,4),~] = Classification_original(CN_data,EMCI_data,CE_num);
%     [accuracy_el(j,4),sensitivity_el(j,4),specificity_el(j,4),fscore_el(j,4),~] = Classification_original(EMCI_data,LMCI_data,EL_num);
%     [accuracy_lc(j,4),sensitivity_lc(j,4),specificity_lc(j,4),fscore_lc(j,4),~] = Classification_original(CN_data,LMCI_data,LC_num);

    fprintf('Now the code is running in the %d loop; \n',j);
end

[paran_max_index,paran_max] = sort(accuracy_ac(:,1),'descend');
% paran_max_index = mean(paran_max_index);
SVM_weight_multi_max = zeros(720,10);
% for i = 1:10
%     SVM_weight_multi_max(:,i) = abs(mean(SVM_weight_multi{paran_max(i),4}(:,1:10),2));
% end

SVM_weight_multi_max_ce = zeros(720,10);
for i = 1:10
    SVM_weight_multi_max_ce(:,i) = abs(mean(SVM_weight_multi{paran_max(i),2}(:,1:10),2));
end
SVM_weight_multi_max_el = zeros(720,10);
% for i = 1:10
%     SVM_weight_multi_max_el(:,i) = abs(mean(SVM_weight_multi{paran_max(i),3}(:,1:10),2));
% end

accuracy(1,:) = mean(accuracy_ac);
accuracy(2,:) = mean(accuracy_lc);
accuracy(3,:) = mean(accuracy_ce);
accuracy(4,:) = mean(accuracy_el);
accuracy_std(1,:) = std(accuracy_ac);
accuracy_std(2,:) = std(accuracy_lc);
accuracy_std(3,:) = std(accuracy_ce);
accuracy_std(4,:) = std(accuracy_el);

specificity(1,:) = mean(specificity_ac);
specificity(2,:) = mean(specificity_lc);
specificity(3,:) = mean(specificity_ce);
specificity(4,:) = mean(specificity_el);
specificity_std(1,:) = std(specificity_ac);
specificity_std(2,:) = std(specificity_lc);
specificity_std(3,:) = std(specificity_ce);
specificity_std(4,:) = std(specificity_el);

sensitivity(1,:) = mean(sensitivity_ac);
sensitivity(2,:) = mean(sensitivity_lc);
sensitivity(3,:) = mean(sensitivity_ce);
sensitivity(4,:) = mean(sensitivity_el);
sensitivity_std(1,:) = std(sensitivity_ac);
sensitivity_std(2,:) = std(sensitivity_lc);
sensitivity_std(3,:) = std(sensitivity_ce);
sensitivity_std(4,:) = std(sensitivity_el);

fscore(1,:) = mean(fscore_ac);
fscore(2,:) = mean(fscore_lc);
fscore(3,:) = mean(fscore_ce);
fscore(4,:) = mean(fscore_el);
fscore_std(1,:) = std(fscore_ac);
fscore_std(2,:) = std(fscore_lc);
fscore_std(3,:) = std(fscore_ce);
fscore_std(4,:) = std(fscore_el);

performance_fold = [accuracy_ac;accuracy_lc;accuracy_ce;accuracy_el;specificity_ac;specificity_lc;specificity_ce;specificity_el;sensitivity_ac;sensitivity_lc;sensitivity_ce;sensitivity_el;fscore_ac;fscore_lc;fscore_ce;fscore_el];
performance = [accuracy;specificity;sensitivity;fscore];
performance_std = [accuracy_std;specificity_std;sensitivity_std;fscore_std];

for i = 1:16
    for j = 1:3
        p_fold(i,j) = permutationTest(performance_fold(1+paran*(i-1):paran*i,1),performance_fold(1+paran*(i-1):paran*i,j+1),10000);
    end
end

% filename = ['results_',dir_name,'\SVM_weight_multi.mat'];
% save(filename,'SVM_weight_multi');
end