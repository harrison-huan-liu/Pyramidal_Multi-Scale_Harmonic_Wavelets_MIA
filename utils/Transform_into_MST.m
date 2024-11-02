function [MST] = Transform_into_MST(m,n,Graph)
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
% save('results\MST.mat','MST');
end