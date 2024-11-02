function [Centrality_node, Centrality_subject_mean, Centrality_subject_mean_sort] = Hub_Centrality_Calculation(occurrence_matrix, SubjectNum, NodeNum, Graph, metric)
    %% hub identification based on centrality
    switch metric
        % betweenness
        case 'Betweenness'
            Centrality_node = betweenness_wei(occurrence_matrix);
        case 'Pagerank'
            Centrality_node = pagerank_centrality(occurrence_matrix, 0.85, occurrence_matrix)';
        case 'Participation'
            [community_vector, community_statistic] = community_louvain(occurrence_matrix);
            Centrality_node = participation_coef(occurrence_matrix, community_vector, 0)';
        % degree
        case 'Degree'
            Centrality_node = degrees_und(occurrence_matrix);
        % local efficiency
        case 'Efficiency'
            Centrality_node = efficiency_wei(occurrence_matrix,2);
        % clustering coefficient
        case 'Clustering_coef'
            Centrality_node = clustering_coef_wd(occurrence_matrix);
        otherwise
            error('Unknown metric');
    end

    %% hub identification based on subject mean
    Centrality_subject = zeros(SubjectNum, NodeNum);

    for i = 1:SubjectNum
        switch metric
            case 'Betweenness'
                Centrality_subject(i, :) = betweenness_wei(Graph(i).W)';
            case 'Pagerank'
                Centrality_subject(i, :) = pagerank_centrality(Graph(i).W, 0.85, Graph(i).W)';
            case 'Participation'
                [community_vector, community_statistic] = community_louvain(Graph(i).W);
                Centrality_subject(i, :) = participation_coef(Graph(i).W, community_vector, 0)';
            case 'Degree'
                % Graph(i).W(Graph(i).W~=0) = 1; % only for degree, don't need to set threshold for other centrality measures
                Centrality_subject(i, :) = degrees_und(Graph(i).W);
            case 'Efficiency'
                Centrality_subject(i, :) = efficiency_wei(Graph(i).W,2);
            case 'Clustering_coef'
                Centrality_subject(i, :) = clustering_coef_wd(Graph(i).W);
            otherwise
                error('Unknown metric');
        end
    end

    Centrality_subject_mean = mean(Centrality_subject);
    [Centrality_subject_mean_rank, Centrality_subject_mean_sort] = sort(Centrality_subject_mean, 'descend');
end