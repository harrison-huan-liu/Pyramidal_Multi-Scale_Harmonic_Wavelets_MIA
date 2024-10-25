close all;clear;clc
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
occurrence_matrix = MST{m+1};

x = 1;
while x <= n-1
    [maximum_column,index_column] = max(occurrence_matrix,[],1);
    [maximun,index] = max(maximum_column);
    occurrence_matrix(index_column(index),index) = 0;
    MST{m+2}(index_column(index),index) = maximun;
    MST{m+2}(index,index_column(index)) = maximun;
    y = MST{m+2};
    for k = 2:size(MST{m+2}, 1)
        y = y | (y*MST{m+2});
    end

    if any(diag(y))
        accu = 1;
    else
        accu = 0;
    end

    if accu == 0 %the reason of 0
        MST{m+2}(index_column(index),index) = 0;
        MST{m+2}(index,index_column(index)) = 0;
        x = x-1;
    end
    x = x+1;
end
MST_occurrence = MST{m+2};

%% HUB
hub_nodes = [3, 6, 10, 15, 16, 23, 24, 27, 30, 35, 44, 48, 51, 55, 56, 59, 65, 66, 68, 69, 73, 77, 80, 84, 89, 90, 93, 97, 98, 104, 118, 122, 123, 127, 128, 129, 140, 142, 144, 147];
node_idx = cell(40,1);
region_mask = cell(40,1);
region_mask_scale_1 = cell(40,1);
region_mask_scale_3 = cell(40,1);
node_idx_scale_1 = cell(40,1);
node_idx_scale_3 = cell(40,1);
node_idx_scale_4 = cell(40,1);
node_idx_scale_5 = cell(40,1);
for i = 1:40
    node_idx{i} = hub_nodes(i);
    node_idx_scale_1_un = find(MST{m+1}(:,hub_nodes(i)));
    node_idx_scale_1{i} = union(node_idx{i},node_idx_scale_1_un);
    for j = 1:148
        if MST{m+1}(j,hub_nodes(i)) ~= 0
            node_idx1 = j;
            node_idx2 = find(MST{m+1}(:,j));
            node_idx3 = union(node_idx1,node_idx2);
            node_idx{i} = union(node_idx{i},node_idx3);
        end
    end
    
    for k = 1:size(node_idx{i})
        node_idx4 = node_idx{i}(k);
        node_idx5 = find(MST{m+1}(:,node_idx{i}(k)));
        node_idx6 = union(node_idx4,node_idx5);
        node_idx_scale_3{i} = union(node_idx_scale_3{i},node_idx6);
    end
    
    for l = 1:size(node_idx_scale_3{i})
        node_idx7 = node_idx_scale_3{i}(l);
        node_idx8 = find(MST{m+1}(:,node_idx_scale_3{i}(l)));
        node_idx9 = union(node_idx7,node_idx8);
        node_idx_scale_4{i} = union(node_idx_scale_4{i},node_idx9);
    end
    
    for z = 1:size(node_idx_scale_4{i})
        node_idx10 = node_idx_scale_4{i}(z);
        node_idx11 = find(MST{m+1}(:,node_idx_scale_4{i}(z)));
        node_idx12 = union(node_idx10,node_idx11);
        node_idx_scale_5{i} = union(node_idx_scale_5{i},node_idx12);
    end
    region_mask{i} = zeros(148,148);
    region_mask_scale_1{i} = zeros(148,148);
    region_mask_scale_3{i} = zeros(148,148);
    for k = 1:size(node_idx_scale_3{i})
        region_mask{i}(:,node_idx_scale_3{i}) = MST{m+1}(:,node_idx_scale_3{i});
        region_mask{i}(node_idx_scale_3{i},:) = MST{m+1}(node_idx_scale_3{i},:);
    end
    region_mask{i}(find(region_mask{i}(:,:))) = 1;
    for k = 1:size(node_idx_scale_1{i})
        region_mask_scale_1{i}(:,node_idx_scale_1{i}) = MST{m+1}(:,node_idx_scale_1{i});
        region_mask_scale_1{i}(node_idx_scale_1{i},:) = MST{m+1}(node_idx_scale_1{i},:);
    end
    region_mask_scale_1{i}(find(region_mask_scale_1{i}(:,:))) = 1;
    for k = 1:size(node_idx_scale_5{i})
        region_mask_scale_3{i}(:,node_idx_scale_5{i}) = MST{m+1}(:,node_idx_scale_5{i});
        region_mask_scale_3{i}(node_idx_scale_5{i},:) = MST{m+1}(node_idx_scale_5{i},:);
    end
    region_mask_scale_3{i}(find(region_mask_scale_3{i}(:,:))) = 1;
end
save('region_mask.mat','region_mask');
save('region_mask_scale_1.mat','region_mask_scale_1');
save('region_mask_scale_3.mat','region_mask_scale_3');
save('node_idx.mat','node_idx');
save('node_idx_scale_1.mat','node_idx_scale_1');
save('node_idx_scale_3.mat','node_idx_scale_3');
save('MST.mat','MST');

u_vector = cell(40,1);
u_vector_scale_1 = cell(40,1);
u_vector_scale_3 = cell(40,1);
for i = 1:40
    u_vector{i} = zeros(148,1);
    u_vector_scale_1{i} = zeros(148,1);
    u_vector_scale_3{i} = zeros(148,1);
    for j = 1:size(node_idx_scale_3{i})
        u_vector{i}(node_idx_scale_3{i}(j),1) = 1;
    end
    for j = 1:size(node_idx_scale_1{i})
        u_vector_scale_1{i}(node_idx_scale_1{i}(j),1) = 1;
    end
    for j = 1:size(node_idx_scale_5{i})
        u_vector_scale_3{i}(node_idx_scale_5{i}(j),1) = 1;
    end
end
save('u_vector.mat','u_vector');
save('u_vector_scale_1.mat','u_vector_scale_1');
save('u_vector_scale_3.mat','u_vector_scale_3');
