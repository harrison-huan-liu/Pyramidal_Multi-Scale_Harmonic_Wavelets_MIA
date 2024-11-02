function [MST_Tree,CommonNetwork] = Calculate_MST_Tree(CommonNetwork,SubjectNum,NodeNum,Graph)
%% this function is used to calculate the MST tree of CommonNetwork, and the NodeNum is the node number of Input Graph, which is CommonNetwork.
MST_Tree = zeros(NodeNum,NodeNum);
com_x = 1;
while com_x <= NodeNum-1
    [maximum_column,index_column] = max(CommonNetwork,[],1);
    [maximun,index] = max(maximum_column);
    CommonNetwork(index_column(index),index) = 0;
    MST_Tree(index_column(index),index) = maximun;
    MST_Tree(index,index_column(index)) = maximun;
    z = MST_Tree;
    for k = 2:size(MST_Tree, 1) % ???
        z = z | (z*MST_Tree);
    end

    if any(diag(z))
        accz = 1;
    else
        accz = 0;
    end

    if accz == 0 % the reason of 0
        MST_Tree(index_column(index),index) = 0;
        MST_Tree(index,index_column(index)) = 0;
        com_x = com_x-1;
    end
    com_x = com_x+1;
end
MST_Tree = (MST_Tree+MST_Tree')/2;

%% there are two ways to calculate the MST tree: 1. average the network first and then calculate the MST; 2. calculate the MST trees of each network first and then average all the MST trees. In there, we save those two MST tree into adjMatrix1.txt and adjMatrix2.txt respectively, and compare those two matrices in adjMatrix3.txt.
% adjacency(G) return the sparse matrix of adjacency matrix

Test_two_ways_of_calculating_MST = 0;

if Test_two_ways_of_calculating_MST
    % First way:
    MST_Tree_a = minspantree(graph(CommonNetwork~=0));

    % Second way:
    occurrence_matrix = zeros(NodeNum,NodeNum);
    for i = 1:SubjectNum
        occurrence_matrix = occurrence_matrix+full(adjacency(minspantree(graph(Graph(i).W~=0))));
    end
    occurrence_matrix=occurrence_matrix/SubjectNum;
    MST_Tree_b=minspantree(graph(occurrence_matrix~=0));

    % adjMatrix_a = graph(MST_Tree_a.Edges(:,1),MST_Tree_a.Edges(:,2));
    % adjMatrix_b = graph(MST_Tree_b.Edges(:,1),MST_Tree_b.Edges(:,2));
    adjMatrix_a = full(adjacency(MST_Tree_a));
    adjMatrix_b = full(adjacency(MST_Tree_b));
    adjMatrix_overlap = adjMatrix_a & adjMatrix_b;
    adjMatrix_overlap = double(adjMatrix_overlap);
    save ../results/Generate_MST/adjMatrix_a.txt -ascii adjMatrix_a
    save ../results/Generate_MST/adjMatrix_b.txt -ascii adjMatrix_b
    save ../results/Generate_MST/adjMatrix_overlap.txt -ascii adjMatrix_overlap

    save('adjMatrix_a.mat','adjMatrix_a','-v6')
    save('adjMatrix_b.mat','adjMatrix_b','-v6')
end

end