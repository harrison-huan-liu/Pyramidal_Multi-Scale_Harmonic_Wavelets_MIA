%% set path for utils
addpath(fullfile(pwd,'..','utils'))

%% load brain connectivity network
filter_label = 0;
[Graph, CommonNetwork, LatentLaplacian, Phi_ave, SubjectNum, NodeNum] = Load_Network_Data_and_Calculate_Laplacian(filter_label);

load('../results/Generate_MST/Adjacency_MST.mat');

%% calculate L_MST % CommonNetwork = Adjacency_MST{SubjectNum+2};

%% load Mask data
load('../results/MST_mask/u_vector.mat');
load('../results/MST_mask/region_mask.mat');
load('../results/MST_mask/node_idx.mat');

%% initialize parameter
[nhub, paran, u1, u2, scale, wavelet_num, p, q, hub_nodes, gama, alpha] = parameter_initialize(u_vector, region_mask, node_idx, LatentLaplacian, Phi_ave);

%% validate inputs
assert(isnumeric(nhub) && nhub > 0, 'nhub must be a positive number');
assert(isnumeric(paran) && paran > 0, 'paran must be a positive number');
assert(isnumeric(u1) && u1 > 0, 'u1 must be a positive number');
assert(isnumeric(u2) && u2 > 0, 'u2 must be a positive number');
assert(isnumeric(scale) && scale > 0, 'scale must be a positive number');
assert(isnumeric(wavelet_num) && wavelet_num > 0, 'wavelet_num must be a positive number');
assert(isnumeric(p) && all(p(:) > 0), 'p must be a positive number');
assert(isnumeric(q) && all(q(:) > 0), 'q must be a positive number');
assert(isnumeric(hub_nodes) && all(hub_nodes(:) > 0), 'hub_nodes must be a positive number');
assert(isnumeric(gama) && all(gama(:) > 0), 'gama must be a positive number');
assert(isnumeric(alpha) && alpha > 0, 'alpha must be a positive number');

%% initialize phi_k in three scales
Theta_k = cell(nhub, scale);
phi_k = cell(nhub, scale + 1);

for i = 1:nhub
    for j = 1:scale
        Theta_k{i,j} = LatentLaplacian + ...
        u1 .* diag(ones(size(u_vector{i,j})) - u_vector{i,j}) + ...
        u2 .* (Phi_ave * Phi_ave.');
        [phi_k{i,j}, ~] = eig(Theta_k{i,j});
        phi_k{i,scale+1} = [phi_k{i,scale+1} phi_k{i,j}(:,1:p(i,j))];
    end
end

%% solve phi
[fixed_item, limited_item, mask_item, phi_iter, identify_wavelets_orthogonal, loss, function_value, iteration_num, multiPhi, singlePhi] = Multiscale_harmonic_wavelets_calculation(nhub, scale, alpha, LatentLaplacian, NodeNum, Phi_ave, phi_k, u_vector, u1, u2, gama);

%% load all suvr data
[Classification_Group, SUVR_name_set, Phi_set, All_SUVR] = load_suvr();

%% validate suvr data
assert(iscell(Classification_Group) && ~isempty(Classification_Group), 'Classification_Group must be a non-empty cell array');
assert(iscell(SUVR_name_set) && ~isempty(SUVR_name_set), 'SUVR_name_set must be a non-empty cell array');
assert(iscell(Phi_set) && ~isempty(Phi_set), 'Phi_set must be a non-empty cell array');
assert(iscell(All_SUVR) && ~isempty(All_SUVR), 'All_SUVR must be a non-empty cell array');

%% svm classification
[Beta, p_value, SVM_para, node_index_index] = svm_classification(Phi_set, SUVR_name_set, Classification_Group, hub_nodes, scale, wavelet_num, p, paran, All_SUVR);

%% post analysis of classification results
[results, brain_node] = post_analysis(SubjectNum, NodeNum, Graph, Adjacency_MST, hub_nodes, node_index_index);

%% power analysis of wavelets
[energy, consider_group] = power_analysis(SUVR_name_set, multiPhi, nhub, scale, wavelet_num, All_SUVR);

%% Visualization and Orthogonality validation of wavelets
[draw_wavelet, ortho] = wavelet_analysis(NodeNum, nhub, scale, multiPhi, Phi_ave, p);
