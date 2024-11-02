function [hub_nodes,CommonNetwork,node_name_left,node_name_right] = Hub_identification(m,n,MST)
CommonNetwork = zeros(n);
for i = 1:m
    CommonNetwork=CommonNetwork+MST{i};
end
CommonNetwork=CommonNetwork/m;
% CommonNetwork(CommonNetwork<0.02)=0;
CommonNetwork=(CommonNetwork+CommonNetwork')/2;

%% Betweenness,PageRank,Participation coefficient,Degree
Betweenness_node = betweenness_wei(CommonNetwork);
Betweenness_colmin = min(Betweenness_node);
Betweenness_colmax = max(Betweenness_node);
Betweenness = rescale(Betweenness_node,'InputMin',Betweenness_colmin,'InputMax',Betweenness_colmax);

Pagerank_node = pagerank_centrality(CommonNetwork,0.85,CommonNetwork);
Pagerank_colmin = min(Pagerank_node);
Pagerank_colmax = max(Pagerank_node);
Pagerank = rescale(Pagerank_node,'InputMin',Pagerank_colmin,'InputMax',Pagerank_colmax);

[community_vector, community_statistic] = community_louvain(CommonNetwork);
Participation_coefficient_node = participation_coef(CommonNetwork,community_vector,0);
Participation_coefficient_colmin = min(Participation_coefficient_node);
Participation_coefficient_colmax = max(Participation_coefficient_node);
Participation_coefficient = rescale(Participation_coefficient_node,'InputMin',Participation_coefficient_colmin,'InputMax',Participation_coefficient_colmax);

% CommonNetwork(CommonNetwork~=0) = 1;
Degree_node = degrees_und(CommonNetwork);
Degree_colmin = min(Degree_node);
Degree_colmax = max(Degree_node);
Degree = rescale(Degree_node,'InputMin',Degree_colmin,'InputMax',Degree_colmax);

%% select hub nodes
Rank_select = Betweenness+Pagerank+Participation_coefficient+Degree';
Rank_select_sort = sort(Rank_select,'descend');
hub_nodes = find(Rank_select>=Rank_select_sort(40));

%% save results
[~,Rank_select_rank] = sort(Rank_select,'descend');
Betweenness_node_norm = Betweenness(Rank_select_rank);
Pagerank_node_norm = Pagerank(Rank_select_rank);
Participation_node_norm = Participation_coefficient(Rank_select_rank);
Degree_node_norm = Degree(Rank_select_rank);

% writematrix(Betweenness_node_norm','Betweenness_node_norm.xls')
% writematrix(Pagerank_node_norm','Pagerank_node_norm.xls')
% writematrix(Participation_node_norm','Participation_node_norm.xls')
% writematrix(Degree_node_norm,'Degree_node_norm.xls')
% writematrix(Rank_select_rank,'Rank_select_rank.xls')

%% draw picture
node_name_left = cell(40,1);
node_name_right = cell(40,1);
destriux_148 = readcell('destriux_148.xlsx');
lloop = 0;
rloop = 0;
node_l = cell(40,1);
node_r = cell(40,1);
for i = 1:148
    if Rank_select_rank(i)<=74&&lloop<40
        lloop = lloop + 1;
        node_name_left{lloop,1} = destriux_148{Rank_select_rank(i),6};
        node_l{lloop,1} = destriux_148{Rank_select_rank(i),1}/1.1;
        node_l{lloop,2} = destriux_148{Rank_select_rank(i),2}/1.1;
        node_l{lloop,3} = destriux_148{Rank_select_rank(i),3}/1.1;
        node_l{lloop,4} = 2;
        node_l{lloop,5} = Rank_select_sort(i);
        node_l{lloop,6} = destriux_148{Rank_select_rank(i),6};
    elseif Rank_select_rank(i)>74&&rloop<40
        rloop = rloop + 1;
        node_name_right{rloop,1} = destriux_148{Rank_select_rank(i),6};
        node_r{rloop,1} = destriux_148{Rank_select_rank(i),1}/1.1;
        node_r{rloop,2} = destriux_148{Rank_select_rank(i),2}/1.1;
        node_r{rloop,3} = destriux_148{Rank_select_rank(i),3}/1.1;
        node_r{rloop,4} = 2;
        node_r{rloop,5} = Rank_select_sort(i);
        node_r{rloop,6} = destriux_148{Rank_select_rank(i),6};
    end
end

% writecell(node_l,'node_l.txt','Delimiter','tab')
% writecell(node_r,'node_r.txt','Delimiter','tab')
end