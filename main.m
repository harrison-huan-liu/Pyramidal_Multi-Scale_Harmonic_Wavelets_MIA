close all;clear;clc
p = 60;% p is the number of eignvectors %60
Graph = Preprocess_data_np_noTransM_01('Data/DataTS-processed.xlsx',p);%Graph包括关于脑连接信息的邻接矩阵，度矩阵，拉普拉斯矩阵及其对应的特征向量矩阵
m = size(Graph,2);%m是研究对象的数目
n = size(Graph(1).L,1);%每个研究对象的结点数目
MST = cell(m+2,1);

SubjectNum=size(Graph,2);
NodeNum=size(Graph(1).W,1);
CommonNetwork=zeros(NodeNum);
for i=1:SubjectNum
    CommonNetwork=CommonNetwork+Graph(i).W;
end
CommonNetwork=CommonNetwork/SubjectNum;
CommonNetwork(CommonNetwork<0.02)=0;
CommonNetwork=(CommonNetwork+CommonNetwork')/2;

MST{97}=zeros(148,148);
com_x = 1;
while com_x <= n-1
    [maximum_column,index_column] = max(CommonNetwork,[],1);
    [maximun,index] = max(maximum_column);
    CommonNetwork(index_column(index),index) = 0;
    MST{97}(index_column(index),index) = maximun;
    MST{97}(index,index_column(index)) = maximun;
    z = MST{97};
    for k = 2:size(MST{97}, 1) % ???
        z = z | (z*MST{97});
    end

    if any(diag(z))
        accz = 1;
    else
        accz = 0;
    end

    if accz == 0 %the reason of 0
        MST{97}(index_column(index),index) = 0;
        MST{97}(index,index_column(index)) = 0;
        com_x = com_x-1;
    end
    com_x = com_x+1;
end

a=minspantree(graph(CommonNetwork~=0));
occurrence_matrix=zeros(148,148);
for i=1:m
    occurrence_matrix=occurrence_matrix+full(adjacency(minspantree(graph(Graph(i).W~=0))));
end
occurrence_matrix=occurrence_matrix/m;
b=minspantree(graph(occurrence_matrix~=0));
% G1 = graph(a.Edges(:,1),a.Edges(:,2));
adjMatrix1 = full(adjacency(a)); % adjacency(G)返回邻接矩阵的稀疏矩阵
% G2 = graph(b.Edges(:,1),b.Edges(:,2));
adjMatrix2 = full(adjacency(b)); % adjacency(G)返回邻接矩阵的稀疏矩阵
save('adjMatrix1.mat','adjMatrix1','-v6')
save('adjMatrix2.mat','adjMatrix2','-v6')
adjMatrix3 = adjMatrix1 & adjMatrix2;
adjMatrix3 = double(adjMatrix3);
save adjMatrix1.txt -ascii adjMatrix1
save adjMatrix2.txt -ascii adjMatrix2
save adjMatrix3.txt -ascii adjMatrix3

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
% save('MST.mat','MST')

for i = 1:m
    MST{95} = MST{95}() + MST{i}();
end
occurrence_matrix = MST{95};
occurrence_matrix=(occurrence_matrix+occurrence_matrix')/2;

x = 1;
while x <= n-1
    [maximum_column,index_column] = max(occurrence_matrix,[],1);
    [maximun,index] = max(maximum_column);
    occurrence_matrix(index_column(index),index) = 0;
    MST{96}(index_column(index),index) = maximun;
    MST{96}(index,index_column(index)) = maximun;
    y = MST{96};
    for k = 2:size(MST{96}, 1)
        y = y | (y*MST{96});
    end

    if any(diag(y))
        accu = 1;
    else
        accu = 0;
    end

    if accu == 0 %the reason of 0
        MST{96}(index_column(index),index) = 0;
        MST{96}(index,index_column(index)) = 0;
        x = x-1;
    end
    x = x+1;
end
MST_occurrence = MST{96};

Degree = zeros(148,1);
for i = 1:n
    row = nnz(MST_occurrence(i,:));
    column = nnz(MST_occurrence(:,i));
    Degree(i,1) = (row+column)/2;
end

d_mean = mean(Degree);
d_index = find(Degree>=2);