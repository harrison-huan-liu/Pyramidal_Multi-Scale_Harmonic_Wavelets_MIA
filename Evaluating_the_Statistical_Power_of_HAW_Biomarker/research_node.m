function [Betweenness_node_select,Pagerank_node_select,Participation_coefficient_node_select,Degree_node_select,brain_node,Pselect_node_index_index_name_order,Bselect_node_index_index_name_order] = research_node(CommonNetwork,select_scale_index,select_node_index_index,select_node_index_index_name,dir_name,MST,m,n)
CommonNetwork = zeros(n);
for i = 1:m
    CommonNetwork=CommonNetwork+MST{i};
end
CommonNetwork=CommonNetwork/m;
% CommonNetwork(CommonNetwork<0.02)=0;
CommonNetwork=(CommonNetwork+CommonNetwork')/2;

select_node_index_index_order = zeros(15,1);
select_node_index_index_name_order = cell(15,1);
j=0;
k=0;
l=0;
for i = 1:15
    if select_scale_index(1,i) == 1
        j=j+1;
        select_node_index_index_order(j,1) = select_node_index_index(i,1);
        select_node_index_index_name_order{j,1} = select_node_index_index_name{i,1};
    elseif select_scale_index(1,i) == 2
        k=k+1;
        select_node_index_index_order(k+5,1) = select_node_index_index(i,1);
        select_node_index_index_name_order{k+5,1} = select_node_index_index_name{i,1};
    else
        l=l+1;
        select_node_index_index_order(l+10,1) = select_node_index_index(i,1);
        select_node_index_index_name_order{l+10,1} = select_node_index_index_name{i,1};
    end
end
Betweenness_node = betweenness_wei(CommonNetwork);
Pagerank_node = efficiency_wei(CommonNetwork,2);
Participation_coefficient_node = clustering_coef_wd(CommonNetwork);
Degree_node = degrees_und(CommonNetwork);
Betweenness_node_select = Betweenness_node(select_node_index_index_order);
Pagerank_node_select = Pagerank_node(select_node_index_index_order);
Participation_coefficient_node_select = Participation_coefficient_node(select_node_index_index_order);
Degree_node_select = Degree_node(select_node_index_index_order);

Pagerank_node_select_order = zeros(15,1);
Betweenness_node_select_order = zeros(15,1);
for i = 1:3
    [Pagerank_node_select(1+5*(i-1):5+5*(i-1)),Pagerank_node_select_order(1+5*(i-1):5+5*(i-1))] = sort(Pagerank_node_select(1+5*(i-1):5+5*(i-1)),'descend');
    [Betweenness_node_select(1+5*(i-1):5+5*(i-1)),Betweenness_node_select_order(1+5*(i-1):5+5*(i-1))] = sort(Betweenness_node_select(1+5*(i-1):5+5*(i-1)));
    Pagerank_node_select_order(1+5*(i-1):5+5*(i-1)) = Pagerank_node_select_order(1+5*(i-1):5+5*(i-1)) + 5*(i-1);
    Betweenness_node_select_order(1+5*(i-1):5+5*(i-1)) = Betweenness_node_select_order(1+5*(i-1):5+5*(i-1)) + 5*(i-1);
end
Pselect_node_index_index_order = select_node_index_index_order(Pagerank_node_select_order);
Pselect_node_index_index_name_order = select_node_index_index_name_order(Pagerank_node_select_order);
Bselect_node_index_index_order = select_node_index_index_order(Betweenness_node_select_order);
Bselect_node_index_index_name_order = select_node_index_index_name_order(Betweenness_node_select_order);

brain_node = zeros(148,1);
for i = 1:15
    brain_node(select_node_index_index(i),1) = 1;
end
filename_node = ['results_',dir_name,'\brain_significant_node.txt'];
save(filename_node,'brain_node','-ascii');

% filename_betweeness = ['results_',dir_name,'\betweenness_wei0620.xls'];
% writematrix(Betweenness_node_select.',filename_betweeness,'WriteMode','append')
% 
% filename_degree = ['results_',dir_name,'\degrees_und0620.xls'];
% writematrix(Degree_node_select,filename_degree,'WriteMode','append')
% 
% filename_pagerank = ['results_',dir_name,'\efficiency_wei0620.xls'];
% writematrix(Pagerank_node_select.',filename_pagerank,'WriteMode','append')
% 
% filename_participation = ['results_',dir_name,'\clustering_coef_wd0620.xls'];
% writematrix(Participation_coefficient_node_select.',filename_participation,'WriteMode','append')
% 
% filename1 = ['results_',dir_name,'\Pselect_node_index_index_name_order.mat'];
% save(filename1,'Pselect_node_index_index_name_order');
% 
% filename2 = ['results_',dir_name,'\Bselect_node_index_index_name_order.mat'];
% save(filename2,'Bselect_node_index_index_name_order');
end