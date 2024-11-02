
Graph = Preprocess_data_np_noTransM_01('../Data/DataTS-processed.xlsx',p);

MST = cell(m+2,1);


% G1 = graph(a.Edges(:,1),a.Edges(:,2));
adjMatrix1 = full(adjacency(a)); 
% G2 = graph(b.Edges(:,1),b.Edges(:,2));
adjMatrix2 = full(adjacency(b)); 
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
