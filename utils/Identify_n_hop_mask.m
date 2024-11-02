function [node_idx_scale,region_mask_scale] = Identify_n_hop_mask(hub_node_idx, occurrence_matrix, scale)
    %% node_idx_scale
    node_idx_temp = hub_node_idx;
    node_idx_scale = cell{1,scale};
    for i = 1:scale
        if i > 1
            node_idx_temp = node_idx_scale{1,i-1}

            for j = 1:size(node_idx_temp)
                node_idx_a = node_idx_temp(j);
                node_idx_b = find(occurrence_matrix(:,node_idx_temp(j)));
                node_idx_c = union(node_idx_a,node_idx_b);
                node_idx_scale{1,i} = union(node_idx_scale{1,i},node_idx_c);
            end
        end
        node_idx_connect_node = find(occurrence_matrix(:,hub_node_idx));
        node_idx_scale{1,i} = union(node_idx_temp,node_idx_connect_node);
    end

    %% region_mask_scale
    region_mask_scale = cell{1,scale};
    for k = 1:scale
        region_mask_scale{1,k} = zeros(size(occurrence_matrix,1),size(occurrence_matrix,1));

        region_mask_scale{1,k}(:,node_idx_scale{1,k}) = occurrence_matrix(:,node_idx_scale{1,k});
        region_mask_scale{1,k}(node_idx_scale{1,k},:) = occurrence_matrix(node_idx_scale{1,k},:);
        region_mask_scale{1,k}(find(region_mask_scale{1,k}(:,:))) = 1;
    end

end