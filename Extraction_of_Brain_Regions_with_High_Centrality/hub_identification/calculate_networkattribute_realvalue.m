function [Betweenness_node,Pagerank_node,Participation_node,Degree_node,Efficiency_local_node,Clustering_coef_node]=calculate_networkattribute_realvalue(Graph)
m = size(Graph,2);%m是研究对象的数目
n = size(Graph(1).L,1);%每个研究对象的结点数目

Amyloid_node_index=[24,98,142,88,122,52,92,136,83,84,12,26,68,8,10];
Tau_node_index=[24,98,142,136,128,92,111,45,10,56,68,122,8,90,3];

load('results_amyloid\select_node_index_index.mat')
Betweenness_node_norm_amy=Betweenness_node(select_node_index_index)/max(Betweenness_node);
Pagerank_node_norm_amy=Pagerank_node(select_node_index_index)/max(Pagerank_node);
Participation_node_norm_amy=Participation_node(select_node_index_index)/max(Participation_node);
Degree_node_norm_amy=Degree_node(select_node_index_index)/max(Degree_node);
load('results_tau\select_node_index_index.mat')

Betweenness_node_norm_amy=Betweenness_node(Amyloid_node_index)/max(Betweenness_node);
Pagerank_node_norm_amy=Pagerank_node(Amyloid_node_index)/max(Pagerank_node);
Participation_node_norm_amy=Participation_node(Amyloid_node_index)/max(Participation_node);
Degree_node_norm_amy=Degree_node(Amyloid_node_index)/max(Degree_node);
Betweenness_node_norm_tau=Betweenness_node(Tau_node_index)/max(Betweenness_node);
Pagerank_node_norm_tau=Pagerank_node(Tau_node_index)/max(Pagerank_node);
Participation_node_norm_tau=Participation_node(Tau_node_index)/max(Participation_node);
Degree_node_norm_tau=Degree_node(Tau_node_index)/max(Degree_node);

xlswrite('results_revised\Betweenness_node_norm_amy.xls',Betweenness_node_norm_amy')
xlswrite('results_revised\Pagerank_node_norm_amy.xls',Pagerank_node_norm_amy)
xlswrite('results_revised\Participation_node_norm_amy.xls',Participation_node_norm_amy)
xlswrite('results_revised\Degree_node_norm_amy.xls',Degree_node_norm_amy)
xlswrite('results_revised\Betweenness_node_norm_tau.xls',Betweenness_node_norm_tau')
xlswrite('results_revised\Pagerank_node_norm_tau.xls',Pagerank_node_norm_tau)
xlswrite('results_revised\Participation_node_norm_tau.xls',Participation_node_norm_tau)
xlswrite('results_revised\Degree_node_norm_tau.xls',Degree_node_norm_tau)

%% the network attributes of commonnetwork
CommonNetwork = zeros(n);
for i = 1:m
    CommonNetwork=CommonNetwork+Graph(i).W;
end

CommonNetwork=CommonNetwork/m;
CommonNetwork(CommonNetwork<0.02)=0;
CommonNetwork=(CommonNetwork+CommonNetwork')/2;

% centrality_metrics:
% degree,PageRank,participation_coefficient,betweenness
% network_attribute:
% betweeness,local_efficiency,global_efficiency(scalar),degree,PageRank,participation_coefficient,cluster_coefficient(optional)
% all_attribute:
Betweenness_node = betweenness_wei(CommonNetwork);
Pagerank_node = pagerank_centrality(CommonNetwork,0.85,CommonNetwork)';
[community_vector, community_statistic] = community_louvain(CommonNetwork);
Participation_node = participation_coef(CommonNetwork,community_vector,0)';
Degree_node = degrees_und(CommonNetwork);
% Efficiency_glocal_node = efficiency_wei(CommonNetwork);
Efficiency_local_node = efficiency_wei(CommonNetwork,2);
Clustering_coef_node = clustering_coef_wd(CommonNetwork);

Betweenness_node_max=max(Betweenness_node);
Pagerank_node_max=max(Pagerank_node);
Participation_node_max=max(Participation_node);
Degree_node_max=max(Degree_node);
Efficiency_local_node_max=max(Efficiency_local_node);

%% the network attributes of MST
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

Betweenness_hub = zeros(m,n);
Pagerank_hub = zeros(m,n);
Participation_hub = zeros(m,n);
Degree_hub = zeros(m,n);

for j = 1:m
    Betweenness_hub(j,:) = betweenness_wei(Graph(j).W)';
    Pagerank_hub(j,:) = pagerank_centrality(Graph(j).W,0.85,Graph(j).W)';
    [community_vector, community_statistic] = community_louvain(Graph(j).W);
    Participation_hub(j,:) = participation_coef(Graph(j).W,community_vector,0)';
    Degree_hub(j,:) = degrees_und(Graph(j).W);
end

Betweenness_hub_mean = mean(Betweenness_hub);
Pagerank_hub_mean = mean(Pagerank_hub);
Participation_hub_mean = mean(Participation_hub);
Degree_hub_mean = mean(Degree_hub);

Betweenness_node = betweenness_wei(MST{m+1});
Pagerank_node = pagerank_centrality(MST{m+1},0.85,MST{m+1})';
[community_vector, community_statistic] = community_louvain(MST{m+1});
Participation_node = participation_coef(MST{m+1},community_vector,0)';
Degree_node = degrees_und(MST{m+1});
Efficiency_local_node = efficiency_wei(MST{m+1},2);
Clustering_coef_node = clustering_coef_wd(MST{m+1});

Betweenness_node_max=max(Betweenness_node);
Pagerank_node_max=max(Pagerank_node);
Participation_node_max=max(Participation_node);
Degree_node_max=max(Degree_node);
Efficiency_local_node_max=max(Efficiency_local_node);
end
