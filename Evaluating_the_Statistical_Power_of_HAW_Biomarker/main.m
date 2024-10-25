% dbstop if all error
close all;clear;clc
%% load Brain Connectivity Network and Amyloid Signdata
p = 55;% p is the number of eignvectors %55
Graph = Preprocess_network_data('../Data/DataTS.csv','../Data/AD-Data/',p);
m = size(Graph,2);% m is the number of reasearch sample
n = size(Graph(1).L,1);% n is the number of node in each reasearch sample
[AD_amy,CN_amy,EMCI_amy,LMCI_amy,SMC_amy] = Proprocess_original_Amyloid_data('../Data/signal_Data/Amyloid_SUVR.xlsx');% Amyloid data
[AD_tau,CN_tau,EMCI_tau,LMCI_tau,SMC_tau] = Proprocess_original_Tau_data('../Data/signal_Data/Tau_SUVR.xlsx');% Tau data

%% COUNT THE DEMOGRAPHIC INFORMATION
[CN_female,CN_male]=count_demographic_information(CN_tau);
[EMCI_female,EMCI_male]=count_demographic_information(EMCI_tau);
[LMCI_female,LMCI_male]=count_demographic_information(LMCI_tau);
FEMALE=[CN_female,EMCI_female,LMCI_female];
MALE=[CN_male,EMCI_male,LMCI_male];
max(FEMALE)
min(FEMALE)
mean(FEMALE)

max(MALE)
min(MALE)
mean(MALE)

ALL=[FEMALE,MALE];
mean(ALL)

%%

LMCI_tau_data = LMCI_tau.Tau;
AD_tau_data = AD_tau.Tau;

% 合并LMCI和AD的数据
LMCI_tau_num = size(LMCI_tau_data,2);
AD_tau_num = size(AD_tau_data,2);
for i = 1:AD_tau_num
    LMCI_tau_data(:,LMCI_tau_num+i)=AD_tau_data(:,i);
end

% [AD_fdg,CN_fdg,EMCI_fdg,LMCI_fdg,SMC_fdg] = Proprocess_original_FDG_data('Data\signal_Data\FDG_SUVR.xlsx');% FDG data

%% Extraction of Brain Regions With High Centrality
% Transform the Adjacency Matrix into MST
[MST] = Transform_into_MST(m,n,Graph);
% identify underlying hubs
[hub_nodes,CommonNetwork,node_name_left,node_name_right] = Hub_identification(m,n,MST);
% writematrix(CommonNetwork,'CommonNetwork.txt','Delimiter','tab');

% aaaa = zeros(148,1);
% for i = 1:40
%     aaaa(hub_nodes(i),1) = 1;
% end

%% Construction of Region-specific Multi-layer Harmonic Wavelets
% Estimation of Common Harmonic Waves
[Graph_MST,LatentLaplacian,CommonHarmonics] = Estimate_common_harmonic_wave(p,m,Graph,MST,CommonNetwork);
% Construction of RM-mask
[region_mask,node_idx,u_vector] = Construct_RM_mask(n,hub_nodes,CommonNetwork);
% Identification of Region-specific Multi-layer Harmonic Wavelets
[Harmonicwavelets,multiPhi,singlePhi] = Identify_RM_harmonic_wavelets(n,hub_nodes,LatentLaplacian,CommonHarmonics,u_vector);

%% identify the hyperparameter
% % 153
% identiedhyperparameter = 400;
% paran_max_index = zeros(identiedhyperparameter,1);
% for i = 1:100
%     u1 = 0.1 * (mod(i-1,10)+1);
%     u2 = 0.1 * (floor((i-1)/10)+1);
%     [Harmonicwavelets,multiPhi,singlePhi] = Identify_RM_harmonic_wavelets(n,hub_nodes,LatentLaplacian,CommonHarmonics,u_vector,u1,u2);
%     [Amy_performance_fold,Amy_performance,Amy_performance_std,Amy_p_fold,Amy_SVM_weight_multi,paran_max_index(i,1)] = Classification_svm(AD_amy.Amyloid,LMCI_amy.Amyloid,EMCI_amy.Amyloid,CN_amy.Amyloid,multiPhi,singlePhi,CommonHarmonics,'amyloid');
% end
% 
% for i = 101:200
%     u1 = 1 * (mod(i-1,10)+1);
%     u2 = 1 * (floor((i-101)/10)+1);
%     [Harmonicwavelets,multiPhi,singlePhi] = Identify_RM_harmonic_wavelets(n,hub_nodes,LatentLaplacian,CommonHarmonics,u_vector,u1,u2);
%     [Amy_performance_fold,Amy_performance,Amy_performance_std,Amy_p_fold,Amy_SVM_weight_multi,paran_max_index(i,1)] = Classification_svm(AD_amy.Amyloid,LMCI_amy.Amyloid,EMCI_amy.Amyloid,CN_amy.Amyloid,multiPhi,singlePhi,CommonHarmonics,'amyloid');
% end
% 
% for i = 201:identiedhyperparameter
%     u1 = 10 * (mod(i-1,10)+1);
%     u2 = 10 * (floor((i-201)/10)+1);
%     [Harmonicwavelets,multiPhi,singlePhi] = Identify_RM_harmonic_wavelets(n,hub_nodes,LatentLaplacian,CommonHarmonics,u_vector,u1,u2);
%     [Amy_performance_fold,Amy_performance,Amy_performance_std,Amy_p_fold,Amy_SVM_weight_multi,paran_max_index(i,1)] = Classification_svm(AD_amy.Amyloid,LMCI_amy.Amyloid,EMCI_amy.Amyloid,CN_amy.Amyloid,multiPhi,singlePhi,CommonHarmonics,'amyloid');
% end
% 
% for i = 301:identiedhyperparameter
%     u1 = 10 * (mod(i-1,10)+1);
%     u2 = 10 * (floor((i-201)/10)+1);
%     [Harmonicwavelets,multiPhi,singlePhi] = Identify_RM_harmonic_wavelets(n,hub_nodes,LatentLaplacian,CommonHarmonics,u_vector,u1,u2);
%     [Amy_performance_fold,Amy_performance,Amy_performance_std,Amy_p_fold,Amy_SVM_weight_multi,paran_max_index(i,1)] = Classification_svm(AD_amy.Amyloid,LMCI_amy.Amyloid,EMCI_amy.Amyloid,CN_amy.Amyloid,multiPhi,singlePhi,CommonHarmonics,'amyloid');
% end
% save('paran_max_index.mat','paran_max_index')

%% Appliaction in Network Neuroscience
[Amy_performance_fold,Amy_performance,Amy_performance_std,Amy_p_fold,Amy_SVM_weight_multi,Amy_SVM_weight_multi_ce,Amy_SVM_weight_multi_el] = Classification_svm(AD_amy.Amyloid,LMCI_amy.Amyloid,EMCI_amy.Amyloid,CN_amy.Amyloid,multiPhi,singlePhi,CommonHarmonics,'amyloid');
[Tau_performance_fold,Tau_performance,Tau_performance_std,Tau_p_fold,Tau_SVM_weight_multi,Tau_SVM_weight_multi_ce,Tau_SVM_weight_multi_el] = Classification_svm(AD_tau.Tau,LMCI_tau_data,EMCI_tau.Tau,CN_tau.Tau,multiPhi,singlePhi,CommonHarmonics,'tau');
% [FDG_performance_fold,FDG_performance,FDG_performance_std,FDG_p_fold,FDG_SVM_weight_multi] = Classification_svm(AD_fdg.FDG,LMCI_fdg.FDG,EMCI_fdg.FDG,CN_fdg.FDG,multiPhi,singlePhi,CommonHarmonics);

% [Amy_select_scale,Amy_select_node_index,Amy_select_scale_index,Amy_select_node_index_index,Amy_select_node_index_index_name,Amy_node_scale,Amy_node_index_scale,Amy_multiPhi_scale] = select_node(Amy_SVM_weight_multi,hub_nodes,multiPhi,u_vector,CommonNetwork,'amyloid');
[Amy_select_scale_ce,Amy_select_node_index_ce,Amy_select_scale_index_ce,Amy_select_node_index_index_ce,Amy_select_node_index_index_name_ce,Amy_node_scale_ce,Amy_node_index_scale_ce,Amy_multiPhi_scale_ce] = select_node(Amy_SVM_weight_multi_ce,hub_nodes,multiPhi,u_vector,CommonNetwork,'amyloid');
% [Amy_select_scale_el,Amy_select_node_index_el,Amy_select_scale_index_el,Amy_select_node_index_index_el,Amy_select_node_index_index_name_el,Amy_node_scale_el,Amy_node_index_scale_el,Amy_multiPhi_scale_el] = select_node(Amy_SVM_weight_multi_el,hub_nodes,multiPhi,u_vector,CommonNetwork,'amyloid');

% [Tau_select_scale,Tau_select_node_index,Tau_select_scale_index,Tau_select_node_index_index,Tau_select_node_index_index_name,tau_node_scale,tau_node_index_scale,tau_multiPhi_scale] = select_node(Tau_SVM_weight_multi,hub_nodes,multiPhi,u_vector,CommonNetwork,'tau');
[Tau_select_scale_ce,Tau_select_node_index_ce,Tau_select_scale_index_ce,Tau_select_node_index_index_ce,Tau_select_node_index_index_name_ce,tau_node_scale_ce,tau_node_index_scale_ce,tau_multiPhi_scale_ce] = select_node(Tau_SVM_weight_multi_ce,hub_nodes,multiPhi,u_vector,CommonNetwork,'tau');
% [Tau_select_scale_el,Tau_select_node_index_el,Tau_select_scale_index_el,Tau_select_node_index_index_el,Tau_select_node_index_index_name_el,tau_node_scale_el,tau_node_index_scale_el,tau_multiPhi_scale_el] = select_node(Tau_SVM_weight_multi_el,hub_nodes,multiPhi,u_vector,CommonNetwork,'tau');

[Amy_Betweenness_node_select,Amy_Pagerank_node_select,Amy_Participation_coefficient_node_select,Amy_Degree_node_select,Amy_brain_node,Amy_Pselect_node_index_index_name_order,Amy_Bselect_node_index_index_name_order] = research_node(CommonNetwork,Amy_select_scale_index,Amy_select_node_index_index,Amy_select_node_index_index_name,'amyloid',MST,m,n);
[~,~,~,~,Amy_brain_node_ce,~,~] = research_node(CommonNetwork,Amy_select_scale_index_ce,Amy_select_node_index_index_ce,Amy_select_node_index_index_name_ce,'amyloid_ce',MST,m,n);
[~,~,~,~,Amy_brain_node_el,~,~] = research_node(CommonNetwork,Amy_select_scale_index_el,Amy_select_node_index_index_el,Amy_select_node_index_index_name_el,'amyloid_el',MST,m,n);

[Tau_Betweenness_node_select,Tau_Pagerank_node_select,Tau_Participation_coefficient_node_select,Tau_Degree_node_select,Tau_brain_node,Tau_Pselect_node_index_index_name_order,Tau_Bselect_node_index_index_name_order] = research_node(CommonNetwork,Tau_select_scale_index,Tau_select_node_index_index,Tau_select_node_index_index_name,'tau',MST,m,n);
[~,~,~,~,Tau_brain_node_ce,~,~] = research_node(CommonNetwork,Tau_select_scale_index_ce,Tau_select_node_index_index_ce,Tau_select_node_index_index_name_ce,'tau_ce',MST,m,n);
[~,~,~,~,Tau_brain_node_el,~,~] = research_node(CommonNetwork,Tau_select_scale_index_el,Tau_select_node_index_index_el,Tau_select_node_index_index_name_el,'tau_el',MST,m,n);

brain_node = zeros(148,1);
brain_node_num = 0;
brain_node_ce = zeros(148,1);
brain_node_ce_num = 0;
brain_node_el = zeros(148,1);
brain_node_el_num = 0;
for i = 1:148
    if Amy_brain_node(i)==1&Tau_brain_node(i)==1
        brain_node(i) = 1;
        brain_node_num = brain_node_num+1;
    end
end

for i = 1:148
    if Amy_brain_node_ce(i)==1&Tau_brain_node_ce(i)==1
        brain_node_ce(i) = 1;
        brain_node_ce_num = brain_node_ce_num+1;
    end
end

for i = 1:148
    if Amy_brain_node_el(i)==1&Tau_brain_node_el(i)==1
        brain_node_el(i) = 1;
        brain_node_el_num = brain_node_el_num+1;
    end
end

node_numindex = 1;
SURV_amy_cn = zeros(8,284);
SURV_amy_lmci = zeros(8,226);
SURV_tau_cn = zeros(8,157);
SURV_tau_lmci = zeros(8,110);
for i=1:148
    if brain_node(i) == 1
        SURV_amy_cn(node_numindex,:) = CN_amy.Amyloid(i,:);
        SURV_amy_lmci(node_numindex,:) = LMCI_amy.Amyloid(i,:);
        SURV_tau_cn(node_numindex,:) = CN_tau.Tau(i,:);
        SURV_tau_lmci(node_numindex,:) = LMCI_tau_data(i,:);
        node_numindex = node_numindex + 1;
    end
end

save('SURV_amy_cn.mat','SURV_amy_cn','-v6')
save('SURV_amy_lmci.mat','SURV_amy_lmci','-v6')
save('SURV_tau_cn.mat','SURV_tau_cn','-v6')
save('SURV_tau_lmci.mat','SURV_tau_lmci','-v6')

save('brain_node.txt','brain_node','-ascii');
save('brain_node_ce.txt','brain_node_ce','-ascii');
save('brain_node_el.txt','brain_node_el','-ascii');

% energy
CL_num_amy = size(CN_amy.Amyloid,2)+size(LMCI_amy.Amyloid,2);
CL_num_tau = size(CN_tau.Tau,2)+size(LMCI_tau_data,2);
LC_signdata_amy = [CN_amy.Amyloid,LMCI_amy.Amyloid];
LC_signdata_tau = [CN_tau.Tau,LMCI_tau_data];
[Energy_signdata_norm_num_amy,Energy_signdata_norm_num_square_amy] = energy_cal(CL_num_amy,LC_signdata_amy,Amy_multiPhi_scale,'amyloid');
[Energy_signdata_norm_num_tau,Energy_signdata_norm_num_square_tau] = energy_cal(CL_num_tau,LC_signdata_tau,tau_multiPhi_scale,'tau');
mean(Energy_signdata_norm_num_square_amy(1:284))
mean(Energy_signdata_norm_num_square_amy(285:510))
mean(Energy_signdata_norm_num_square_tau(1:157))
mean(Energy_signdata_norm_num_square_tau(158:267))
std(Energy_signdata_norm_num_square_amy(1:284))
std(Energy_signdata_norm_num_square_amy(285:510))
std(Energy_signdata_norm_num_square_tau(1:157))
std(Energy_signdata_norm_num_square_tau(158:267))

%% Draw picture
save('results20220703/MST.mat','MST')