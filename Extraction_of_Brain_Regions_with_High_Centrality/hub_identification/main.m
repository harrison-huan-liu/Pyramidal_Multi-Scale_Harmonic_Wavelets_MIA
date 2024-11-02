%% set path for utils
addpath(fullfile(pwd,'..','utils'))

%% load brain connectivity network
filter_label = 0;
[Graph,CommonNetwork,LatentLaplacian,Phi_ave,SubjectNum,NodeNum] = Load_Network_Data_and_Calculate_Laplacian(filter_label);

%% MST
load ../results/Generate_MST/Adjacency_MST.mat
occurrence_matrix = Adjacency_MST{SubjectNum+2};

%% hub index
hub_nodes = [3, 6, 10, 15, 16, 23, 24, 27, 30, 35, 44, 48, 51, 55, 56, 59, 65, 66, 68, 69, 73, 77, 80, 84, 89, 90, 93, 97, 98, 118, 122, 123, 127, 128, 129, 104, 140, 142, 144, 147];

metrics = {'Betweenness', 'Pagerank', 'Participation', 'Degree'};
results = struct();

for i = 1:length(metrics)
    metric = metrics{i};
    [Centrality_node, Centrality_subject_mean, Centrality_subject_mean_sort] = Hub_Centrality_Calculation(occurrence_matrix, SubjectNum, NodeNum, Graph, metric);
    results.(metric).node = Centrality_node;
    results.(metric).subject_mean_sort = Centrality_subject_mean_sort;
end

% Normalize and sort centrality measures
[Centrality_node_norm, Sum_norm, Sum_norm_sort, Sum_norm_sort_left, Sum_norm_sort_right] = Centrality_Normalization(results, metrics);

%% select hub nodes
for i = 1:length(metrics)
    metric = metrics{i};
    Centrality_node = results.(metric).node;
    Centrality_subject_mean_sort = results.(metric).subject_mean_sort;

    Centrality_mean = mean(Centrality_node);
    Centrality_sort = sort(Centrality_node, 'descend');
    Centrality_index = find(Centrality_node >= Centrality_sort(12)); % Betweenness_mean+450 % Pagerank_mean+0.0025 % Participation_coefficient_mean+0.23 % Degree_mean+6.2
    
    %% compare the overlap of node index on hub identification based on common network and subject mean
    hub_overlap = ismember(Centrality_subject_mean_sort, hub_nodes);
    hub_overlap_index = zeros(1, NodeNum);
    hub_overlap_index(hub_overlap) = 1;

    writematrix(hub_overlap_index, [metric, '_hub_overlap_index.xls']);
end
