%% set path for utils
addpath(fullfile(pwd,'..','utils'))

%% load brain connectivity network
filter_label = 0;
[Graph,CommonNetwork,LatentLaplacian,Phi_ave,SubjectNum,NodeNum] = Load_Network_Data_and_Calculate_Laplacian(filter_label)

%% load MST_occurrence and MST tree
load ../results/Generate_MST/Adjacency_MST.mat
MST_common = Adjacency_MST{SubjectNum+1};
occurrence_matrix = Adjacency_MST{SubjectNum+2};
MST_occurrence = Adjacency_MST{SubjectNum+3};


%% HUB
hub_nodes = [3, 6, 10, 15, 16, 23, 24, 27, 30, 35, 44, 48, 51, 55, 56, 59, 65, 66, 68, 69, 73, 77, 80, 84, 89, 90, 93, 97, 98, 104, 118, 122, 123, 127, 128, 129, 140, 142, 144, 147];
hub_number = size(hub_nodes);

%% identify node_idx and region_mask
scale = 3;
node_idx = cell(hub_number,scale);
region_mask = cell(hub_number,scale);
for i = 1:hub_number
    [node_idx_scale,region_mask_scale] = Identify_n_hop_mask(hub_nodes(i), occurrence_matrix, scale)
    node_idx{i,:} = node_idx_scale;
    region_mask{i,:} = region_mask_scale;
end
save('../results/MST_mask/region_mask.mat','region_mask');
save('../results/MST_mask/node_idx.mat','node_idx');

u_vector = cell(hub_number,scale);
for j = 1:hub_number
    for k = 1:scale
        u_vector{j,k} = zeros(NodeNum,1);
        u_vector{j,k}(node_idx{j,k},1) = 1;
    end
end
save('../results/MST_mask/u_vector.mat','u_vector');
