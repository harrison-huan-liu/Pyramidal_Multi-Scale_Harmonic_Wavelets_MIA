function [nhub, paran, u1, u2, scale, wavelet_num, p, q, hub_nodes, gama, alpha]=parameter_initialize(u_vector, region_mask, node_idx, LatentLaplacian, Phi_ave)
    %% identified paraemter
    nhub = size(u_vector,1);
    paran = 10; % ten fold
    u1 = 0.8; % scale=1: 3; scale=2: 2; scale=3: 1.
    u2 = 0.7;
    scale = 3;

    %% select the number of wavelets
    wavelet_num = 6;
    p = wavelet_num * ones(nhub,scale);
    q = sum(p,2);

    hub_nodes = [3, 6, 10, 15, 16, 23, 24, 27, 30, 35, 44, 48, 51, 55, 56, 59, 65, 66, 68, 69, 73, 77, 80, 84, 89, 90, 93, 97, 98, 118, 122, 123, 127, 128, 129, 104, 140, 142, 144, 147];

    %% initialize the p and q based on the degree of region_mask
    degree = cell(nhub,scale);
    for i = 1:nhub
        for j = 1:scale
            degree{i,j} = degrees_und(region_mask{i,j});
            if j>1
                node_number_diff_scale = size(node_idx{i,j-1},1);
                node_idx_diff_scale = node_idx{i,j-1};
            else
                node_number_diff_scale = 1;
                node_idx_diff_scale = hub_nodes(i);
            end
            for k = 1:node_number_diff_scale
                p(i,j) = p(i,j) + degree{i,j}(node_idx_diff_scale(k));
            end
            p(i,j) = p(i,j)/node_number_diff_scale;
        end
    end

    %% identify gama and alpha
    gama = cell(nhub,3);
    for i = 1:nhub
        moving_index = 0;
        for j = 1:scale
            gama{i,j} = zeros(q(i),q(i));
            for k = 1:p(i,j)
                gama{i,j}(moving_index+k,moving_index+k) = 1;
            end
            moving_index = moving_index + p(i,j);
        end
    end

    alpha = 0;
    for i = 1:nhub
        fixed_value = LatentLaplacian + u2.*(Phi_ave*Phi_ave.');
        for j = 1:scale
            fixed_value = fixed_value + ...
            u1.*diag(ones(size(u_vector{i,j})) - u_vector{i,j});
        end
        [~,eigenvalue] = eig(fixed_value);
        if max(max(eigenvalue))>alpha
            alpha = max(max(eigenvalue));
        end
    end
end