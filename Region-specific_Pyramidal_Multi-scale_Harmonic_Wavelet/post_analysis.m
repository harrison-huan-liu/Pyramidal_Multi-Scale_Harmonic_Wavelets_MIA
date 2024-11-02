function [results, brain_node]=post_analysis(SubjectNum, NodeNum, Graph, Adjacency_MST, hub_nodes, node_index_index)
    %% Post Analysis
    metrics = {'Betweenness', 'Efficiency', 'Clustering_coef', 'Degree'};
    results = struct();

    CommonNetwork = Adjacency_MST{SubjectNum+2};
    for i = 1:length(metrics)
        metric = metrics{i};
        [Centrality_node, Centrality_subject_mean, Centrality_subject_mean_sort] = Hub_Centrality_Calculation(CommonNetwork, SubjectNum, NodeNum, Graph, metric);
        results.(metric).node = Centrality_node;
        results.(metric).subject_mean = Centrality_subject_mean;
        results.(metric).subject_mean_sort = Centrality_subject_mean_sort;

        hub_centrality_top = Centrality_subject_mean(node_index_index);
        writematrix(hub_centrality_top, ['../results/Region-specific_Pyramidal_Multi-scale_Harmonic_Wavelet/', metric, '_hub_top.xls'],'WriteMode','append');
    end

    for i = 1:length(metrics)
        metric = metrics{i};
        Centrality_node = results.(metric).node;
        Centrality_node_hub_top = Centrality_node(:)(node_index_index);
        results.(metric).Centrality_node_hub_top = Centrality_node_hub_top;
        writematrix(Centrality_node_hub_top, ['../results/Region-specific_Pyramidal_Multi-scale_Harmonic_Wavelet/', metric, '_commonnetwork_hub_top.xls'],'WriteMode','append');
    end

    brain_node = zeros(NodeNum,1);
    brain_node(node_index_index) = 1;
    save('../results/Region-specific_Pyramidal_Multi-scale_Harmonic_Wavelet/brain_significant_node.txt','brain_node','-ascii');
end