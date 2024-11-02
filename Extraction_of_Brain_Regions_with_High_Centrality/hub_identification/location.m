% close all;clear;clc
%% load brain connectivity network
p=55;% p is the number of eignvectors %55
Graph=Preprocess_network_data('../Data/DataTS.csv','../Data/AD-Data/',p);

%% txt
destriux_148 = readcell('destriux_148.xlsx');
hub_nodes = [3, 6, 10, 15, 16, 23, 24, 27, 30, 35, 44, 48, 51, 55, 56, 59, 65, 66, 68, 69, 73, 77, 80, 84, 89, 90, 93, 97, 98, 104, 118, 122, 123, 127, 128, 129, 140, 142, 144, 147, 155];
% x = cell(148,1);
% y = cell(148,1);
% z = cell(148,1);
% color = cell(148,1);
% degree = cell(148,1);
node = destriux_148;
hub_nodes = [Sum_norm_sort_left(1:40),Sum_norm_sort_right(1:40)];

j = 1;
% temp = degrees_und(Graph(1).W);
temp = Sum_norm;
for i = 1:148
    if i == hub_nodes(j)
        node{j,1} = destriux_148{hub_nodes(j),1}/1.1;
        node{j,2} = destriux_148{hub_nodes(j),2}/1.1;
        node{j,3} = destriux_148{hub_nodes(j),3}/1.1;
%         color{i} = 2;
%         node{i,4} = color{i};
%         temp = degrees_und(Graph(1).W);
%         degree{i} = temp(hub_nodes(j));
%         node{i,5} = degree{i};
        node{j,4} = 1.5 + 0.5*rand;
        node{j,5} = temp(hub_nodes(j));
        node{j,6} = destriux_148{hub_nodes(j),6};
        j = j + 1;
    else
%         node{i,1} = destriux_148{i,1};
%         node{i,2} = destriux_148{i,2};
%         node{i,3} = destriux_148{i,3};
% %         color{i} = 1;
% %         node{i,4} = color{i};
% %         temp = degrees_und(Graph(1).W);
% %         degree{i} = temp(hub_nodes(i));
% %         node{i,5} = degree{i};
%         node{i,5} = temp(i);
    end
end

writecell(node,'hub40_randi_520.txt','Delimiter','tab')
% load MST
% writematrix(MST{139},'edge40.txt','Delimiter','tab')

%% degree 148 nodes
p=55;% p is the number of eignvectors %55
Graph=Preprocess_network_data('Data\DataTS.csv','Data\AD-Data\',p);
CommonNetwork = zeros(n);
for i = 1:138
    CommonNetwork=CommonNetwork+Graph(i).W;
end
CommonNetwork=CommonNetwork/m;
CommonNetwork(CommonNetwork<0.002)=0;
CommonNetwork=(CommonNetwork+CommonNetwork')/2;
temp = degrees_und(CommonNetwork);
destriux_148 = readcell('destriux_148.xlsx');
node = destriux_148;
for i = 1:148
    node{i,1} = destriux_148{i,1}/1.1;
    node{i,2} = destriux_148{i,2}/1.1;
    node{i,3} = destriux_148{i,3}/1.1;
    node{i,4} = temp(i)/148;
    node{i,5} = temp(i)/148;
    node{i,6} = destriux_148{i,6};
end

writecell(node,'node148_degree.txt','Delimiter','tab')

load b725
load b727
node_a = cell(9,1);
node_n = zeros(1,9);
node_index = cell(9,1);
for i = 1:9
    node_index{i,1} = find(abs(b725(:,i))>0);
    node_index{i,1} = union(node_index{i,1},b727(i));
    node_a{i,1} = destriux_148(node_index{i,1},:);
    node_n(1,i) = length(node_index{i,1});
    for j = 1:node_n(1,i)
        node_a{i,1}{j,6} = b725(node_index{i,1}(j),i);
        if b727(1,i) == (node_index{i,1}(j))
            node_a{i,1}{j,4} = 2;
            node_a{i,1}{j,5} = 3;
        end
    end
end
for i = 1:9
    filename = ['node_',num2str(i),'.txt'];
    writecell(node_a{i,1},filename,'Delimiter','tab');
end


CommonNetwork = zeros(n);
for i = 1:m
    CommonNetwork=CommonNetwork+Graph(i).W;
end
CommonNetwork=CommonNetwork/m;
CommonNetwork(CommonNetwork<0.002)=0;
CommonNetwork=(CommonNetwork+CommonNetwork')/2;
edge_a = cell(9,1);
for i = 1:9
    edge_a{i,1} = zeros(node_n(1,i),node_n(1,i));
    for j = 1:node_n(1,i)
        for k = 1:node_n(1,i)
            edge_a{i,1}(j,k) = CommonNetwork(node_index{i,1}(j),node_index{i,1}(k));
        end
    end
end
for i = 1:9
    filename = ['edge_',num2str(i),'.txt'];
    writematrix(edge_a{i,1},filename,'Delimiter','tab');
end

sign_data = cell(3,1);
for i = 1:3
    sign_data{i,1} = destriux_148;
    for j = 1:148
        sign_data{i,1}{j,4}=AD_amy.Amyloid(j,i);
        sign_data{i,1}{j,5}=AD_amy.Amyloid(j,i);
    end
end
for i = 1:3
    filename = ['sign_data_',num2str(i),'.txt'];
    writecell(sign_data{i,1},filename,'Delimiter','tab');
end


%% framework picture
b728 = zeros(148,9);
for i = 1:9
    b728(:,i) = multiPhi(:,6*(i-1)+1);
end
b729 = [3,3,3,6,6,6,10,10,10];

% for i = 1:148
%     for j = 1:9
%         if abs(b728(i,j))<0.1
%             b728(i,j) = 0;
%         end
%     end
% end

for i = 1:2
    for j = 1
        if abs(b728(i,j))<0.1
            b728(i,j) = 0;
        end
    end
end

for i = 4:148
    for j = 1
        if abs(b728(i,j))<0.1
            b728(i,j) = 0;
        end
    end
end

for i = 1:5
    for j = 2
        if abs(b728(i,j))<0.1
            b728(i,j) = 0;
        end
    end
end

for i = 7:148
    for j = 2
        if abs(b728(i,j))<0.1
            b728(i,j) = 0;
        end
    end
end

for i = 1:9
    for j = 3
        if abs(b728(i,j))<0.1
            b728(i,j) = 0;
        end
    end
end

for i = 11:148
    for j = 3
        if abs(b728(i,j))<0.1
            b728(i,j) = 0;
        end
    end
end

%% out

fnode_a = cell(9,1);
fnode_n = zeros(1,9);
fnode_index = cell(9,1);
for i = 1:9
    fnode_n(1,i) = length(find(abs(b728(:,i))>0));
end
for i = 1:9
    fnode_index{i,1} = find(abs(b728(:,i))>0);
    fnode_a{i,1} = destriux_148(fnode_index{i,1},:);
    for j = 1:fnode_n(1,i)
        fnode_a{i,1}{j,6} = b728(fnode_index{i,1}(j),i);
        if b729(1,i) == (fnode_index{i,1}(j))
            fnode_a{i,1}{j,4} = 2;
        end
    end
end
for i = 1:9
    filename = ['fnode_',num2str(i),'.txt'];
    writecell(fnode_a{i,1},filename,'Delimiter','tab');
end


CommonNetwork = zeros(n);
for i = 1:m
    CommonNetwork=CommonNetwork+Graph(i).W;
end
CommonNetwork=CommonNetwork/m;
CommonNetwork(CommonNetwork<0.002)=0;
CommonNetwork=(CommonNetwork+CommonNetwork')/2;
fedge_a = cell(9,1);
for i = 1:9
    fedge_a{i,1} = zeros(fnode_n(1,i),fnode_n(1,i));
    for j = 1:fnode_n(1,i)
        for k = 1:fnode_n(1,i)
            fedge_a{i,1}(j,k) = CommonNetwork(fnode_index{i,1}(j),fnode_index{i,1}(k));
        end
    end
end
for i = 1:9
    filename = ['fedge_',num2str(i),'.txt'];
    writematrix(fedge_a{i,1},filename,'Delimiter','tab');
end

save('b728.mat','b728');
save('b729.mat','b729');