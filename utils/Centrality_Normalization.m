function [Centrality_node_norm, Sum_norm, Sum_norm_sort, Sum_norm_sort_left, Sum_norm_sort_right] = Centrality_Normalization(results, metrics)
    % Initialize Sum_norm
    Sum_norm = zeros(size(results.(metrics{1}).node));
    
    % Normalize and reshape centrality measures
    for i = 1:length(metrics)
        metric = metrics{i};
        Centrality_node = results.(metric).node;
        Centrality_node_norm = normalize(Centrality_node, 'range')(:);

        % Sum normalized centrality measures
        Sum_norm = Sum_norm + Centrality_node_norm;
    end
    
    % Sort whole brain centrality
    [~, Sum_norm_sort] = sort(Sum_norm, 'descend');
    
    % Sort left and right brain centrality
    [~, Sum_norm_sort_left] = sort(Sum_norm(1:74), 'descend');
    [~, Sum_norm_sort_right] = sort(Sum_norm(75:end), 'descend');
    
    % Write sorted centrality to file
    writematrix(Sum_norm_sort, '../results/Hub_Identification/Sum_norm_sort.xls');
end