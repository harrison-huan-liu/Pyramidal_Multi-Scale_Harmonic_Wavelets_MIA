%% set path for utils
addpath(fullfile(pwd,'..','utils'))

%% Preprocess Data
% Label SMC to LMCI or only save AD and CN data if filter_label is 1, no change if filter_label is 0
filter_label = 1
[Graph,CommonNetwork,LatentLaplacian,Phi_ave,SubjectNum,NodeNum] = Load_Network_Data_and_Calculate_Laplacian(filter_label)

Adjacency_MST = cell(SubjectNum+3,1);
%% 1. average the network first and then calculate the MST
[MST_Tree,~] = Calculate_MST_Tree(CommonNetwork,SubjectNum,NodeNum,Graph);
Adjacency_MST{SubjectNum+1} = MST_Tree;

%% 2. calculate the MST trees of each network first and then average all the MST trees
for i = 1:SubjectNum
    M_connectivity = tril(Graph(i).W);
    [MST_Tree,~] = Calculate_MST_Tree(M_connectivity,SubjectNum,NodeNum,Graph);
    Adjacency_MST{i} = MST_Tree;
end

Adjacency_MST{SubjectNum+2} = zeros(NodeNum,NodeNum);
for i = 1:SubjectNum
    Adjacency_MST{SubjectNum+2} = Adjacency_MST{SubjectNum+2}() + Adjacency_MST{i}();
end
occurrence_matrix = Adjacency_MST{SubjectNum+2};
occurrence_matrix = (occurrence_matrix+occurrence_matrix')/2;
[MST_occurrence,~] = Calculate_MST_Tree(occurrence_matrix,SubjectNum,NodeNum,Graph);

Adjacency_MST{SubjectNum+3} = MST_occurrence;
% save('../results/Generate_MST/Adjacency_MST.mat','Adjacency_MST')

%% Degree analysis of MST_occurrence and select those nodes with high degree
Degree = zeros(NodeNum,1);
for i = 1:NodeNum
    row = nnz(MST_occurrence(i,:));
    column = nnz(MST_occurrence(:,i));
    Degree(i,1) = (row+column)/2;
end

d_mean = mean(Degree);
d_index = find(Degree>=2);