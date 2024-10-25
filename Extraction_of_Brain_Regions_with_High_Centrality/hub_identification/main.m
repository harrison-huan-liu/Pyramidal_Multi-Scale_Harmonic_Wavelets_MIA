% close all;clear;clc
%% load brain connectivity network
p=55;% p is the number of eignvectors %55
Graph=Preprocess_network_data('../Data/DataTS.csv','../Data/AD-Data/',p);
m = size(Graph,2);
n = size(Graph(1).L,1);
MST = cell(m+2,1);

%% MST
for i = 1:m+2
    MST{i} = zeros(148,148);
end

for i = 1:m
    M_connectivity = tril(Graph(i).W);
    j = 1;
    while j <= n-1
        [maximum_column,index_column] = max(M_connectivity,[],1);
        [maximun,index] = max(maximum_column);
%         temp = M_connectivity(index_column(index),index);
        M_connectivity(index_column(index),index) = 0;
        MST{i}(index_column(index),index) = maximun;
        MST{i}(index,index_column(index)) = maximun;
        t = MST{i};
        for k = 2:size(MST{i}, 1)
            t = t | (t*MST{i});
        end
 
        if any(diag(t))
            acc = 1;
        else
            acc = 0;
        end
        
        if acc == 0 %the reason of 0
            MST{i}(index_column(index),index) = 0;
            MST{i}(index,index_column(index)) = 0;
            j = j-1;
        end
        j = j+1;
    end
%     MST{i} = MST{i} + MST{i}';
end

for i = 1:m
    MST{m+1} = MST{m+1}() + MST{i}();
end



% CommonNetwork = zeros(n);
Betweenness_hub = zeros(m,n);
Pagerank_hub = zeros(m,n);
Participation_hub = zeros(m,n);
Degree_hub = zeros(m,n);
% for i = 1:m
%     CommonNetwork=CommonNetwork+Graph(i).W;
% end

% for j = 1:m
%     Betweenness_hub(j,:) = betweenness_wei(Graph(j).W)';
%     Pagerank_hub(j,:) = pagerank_centrality(Graph(j).W,0.85,Graph(j).W)';
%     [community_vector, community_statistic] = community_louvain(Graph(j).W);
%     Participation_hub(j,:) = participation_coef(Graph(j).W,community_vector,0)';
%     Degree_hub(j,:) = degrees_und(Graph(j).W);
% end

% CommonNetwork=CommonNetwork/m;
% CommonNetwork(CommonNetwork<0.02)=0;
% CommonNetwork=(CommonNetwork+CommonNetwork')/2;
Betweenness_node = betweenness_wei(MST{m+1});
Pagerank_node = pagerank_centrality(MST{m+1},0.85,MST{m+1})';
[community_vector, community_statistic] = community_louvain(MST{m+1});
Participation_node = participation_coef(MST{m+1},community_vector,0)';
Degree_node = degrees_und(MST{m+1});
Efficiency_node = efficiency_wei(MST{m+1},2);
Clustering_coef_node = clustering_coef_wd(MST{m+1});

Betweenness_node_norm = normalize(Betweenness_node,'range');
Pagerank_node_norm = normalize(Pagerank_node,'range');
Participation_node_norm = normalize(Participation_node,'range');
Degree_node_norm = normalize(Degree_node,'range');
Sum_norm = Betweenness_node_norm' + Pagerank_node_norm + Participation_node_norm + Degree_node_norm;
Sum_norm_left = Sum_norm(1:74);
Sum_norm_right = Sum_norm(75:148);
[Sum_norm_rank,Sum_norm_sort] = sort(Sum_norm,'descend');
[~,Sum_norm_sort_left] = sort(Sum_norm_left,'descend');
[~,Sum_norm_sort_right] = sort(Sum_norm_right,'descend');


Betweenness_node_norm = Betweenness_node_norm(Sum_norm_sort);
Pagerank_node_norm = Pagerank_node_norm(Sum_norm_sort);
Participation_node_norm = Participation_node_norm(Sum_norm_sort);
Degree_node_norm = Degree_node_norm(Sum_norm_sort);

writematrix(Betweenness_node_norm','Betweenness_node_norm.xls')
writematrix(Pagerank_node_norm,'Pagerank_node_norm.xls')
writematrix(Participation_node_norm,'Participation_node_norm.xls')
writematrix(Degree_node_norm,'Degree_node_norm.xls')
writematrix(Sum_norm_sort,'Sum_norm_sort.xls')

Betweenness_hub_mean = mean(Betweenness_hub);
Pagerank_hub_mean = mean(Pagerank_hub);
Participation_hub_mean = mean(Participation_hub);
Degree_hub_mean = mean(Degree_hub);

[Betweenness_hub_rank,Betweenness_hub_sort] = sort(Betweenness_hub_mean,'descend');
[Pagerank_hub_rank,Pagerank_hub_sort] = sort(Pagerank_hub_mean,'descend');
[Participation_hub_rank,Participation_hub_sort] = sort(Participation_hub_mean,'descend');
[Degree_hub_rank,Degree_hub_sort] = sort(Degree_hub_mean,'descend');

% for k = 1:n
% end

Betweenness_hub = Betweenness_hub(:,Betweenness_hub_sort);
Pagerank_hub = Pagerank_hub(:,Pagerank_hub_sort);
Participation_hub = Participation_hub(:,Participation_hub_sort);
Degree_hub = Degree_hub(:,Degree_hub_sort);

writematrix(Betweenness_hub,'Betweenness_hub.xls')
writematrix(Pagerank_hub,'Pagerank_hub.xls')
writematrix(Participation_hub,'Participation_hub.xls')
writematrix(Degree_hub,'Degree_hub.xls')

hub_nodes = [3, 6, 10, 15, 16, 23, 24, 27, 30, 35, 44, 48, 51, 55, 56, 59, 65, 66, 68, 69, 73, 77, 80, 84, 89, 90, 93, 97, 98, 118, 122, 123, 127, 128, 129, 104, 140, 142, 144, 147];
Betweenness_hub_is = ismember(Betweenness_hub_sort,hub_nodes);
Pagerank_hub_is = ismember(Pagerank_hub_sort,hub_nodes);
Participation_hub_is = ismember(Participation_hub_sort,hub_nodes);
Degree_hub_is = ismember(Degree_hub_sort,hub_nodes);

Betweenness_hub_int = zeros(1,148);
Pagerank_hub_int = zeros(1,148);
Participation_hub_int = zeros(1,148);
Degree_hub_int = zeros(1,148);

Betweenness_hub_int(Betweenness_hub_is) = 1;
Pagerank_hub_int(Pagerank_hub_is) = 1;
Participation_hub_int(Participation_hub_is) = 1;
Degree_hub_int(Degree_hub_is) = 1;

writematrix(Betweenness_hub_int,'Betweenness_hub_int.xls')
writematrix(Pagerank_hub_int,'Pagerank_hub_int.xls')
writematrix(Participation_hub_int,'Participation_hub_int.xls')
writematrix(Degree_hub_int,'Degree_hub_int.xls')

%% betweenness
% Betweenness_node = zeros(m,n);
% for i=1:m
%     Betweenness_node(i,:) = betweenness_wei(Graph(i).W);
% end
% 
% %% pagerank
% Pagerank_node = zeros(m,n);
% for i=1:m
%     Pagerank_node(i,:) = pagerank_centrality(Graph(i).W,0.85,Graph(i).W);
% end
% 
% %% participation coefficient
% Participation_coefficient_node = zeros(m,n);
% for i = 1:m
%     [community_vector, community_statistic] = community_louvain(Graph(i).W);
%     Participation_coefficient_node(i,:) = participation_coef(Graph(i).W,community_vector,0);
% end
% 
% %% degree
% Degree_node = zeros(m,n);
% for i = 1:m
%     Graph(i).W(Graph(i).W~=0) = 1;
%     Degree_node(i,:) = degrees_und(Graph(i).W);
% end

%% select hub nodes
% Betweenness = mean(Betweenness_node);
% Pagerank = mean(Pagerank_node);
% Participation_coefficient = mean(Participation_coefficient_node);
% Degree = mean(Degree_node);

Betweenness_mean = mean(Betweenness_node);
Pagerank_mean = mean(Pagerank_node);
Participation_coefficient_mean = mean(Participation_coefficient_node);
Degree_mean = mean(Degree_node);

Betweenness_sort = sort(Betweenness_node,'descend');
Pagerank_sort = sort(Pagerank_node,'descend');
Participation_coefficient_sort = sort(Participation_coefficient_node,'descend');
Degree_sort = sort(Degree_node,'descend');

Betweenness_index = find(Betweenness_node>=Betweenness_sort(12));% Betweenness_mean+450
Pagerank_index = find(Pagerank_node>=Pagerank_sort(12));% Pagerank_mean+0.0025
Participation_coefficient_index = find(Participation_coefficient_node>=Participation_coefficient_sort(12));% Participation_coefficient_mean+0.23
Degree_index = find(Degree_node>=Degree_sort(12));% Degree_mean+6.2

hub_nodes = [3, 6, 10, 15, 16, 23, 24, 27, 30, 35, 44, 48, 51, 55, 56, 59, 65, 66, 68, 69, 73, 77, 80, 84, 89, 90, 93, 97, 98, 118, 122, 123, 127, 128, 129, 104, 140, 142, 144, 147];
