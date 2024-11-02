function [region_mask,node_idx_modules,u_vector] = Construct_RM_mask(n,hub_nodes,CommonNetwork)
%% HUB
select_modules = 3;
nhub = size(hub_nodes,1);
region_mask = cell(nhub,select_modules);
node_idx = cell(nhub,select_modules);
% node_idx = cell(nhub,2*select_modules-1); % scale: 1,3,5
node_idx_modules = cell(nhub,select_modules);
for i = 1:nhub
    node_idx_un = find(CommonNetwork(:,hub_nodes(i)));
    node_idx{i,1} = union(hub_nodes(i),node_idx_un);
    for k = 1:select_modules-1 % 2*select_modules-2 % scale: 1,3,5
        for j = 1:size(node_idx{i,k})
            node_idx1 = node_idx{i,k}(j);
            node_idx2 = find(CommonNetwork(:,node_idx{i,k}(j)));
            node_idx3 = union(node_idx1,node_idx2);
            node_idx{i,k+1} = union(node_idx{i,k+1},node_idx3);
        end
    end
    
    for l = 1:select_modules
        node_idx_modules{i,l} = node_idx{i,l}; % 2*l-1 % scale: 1,3,5
        region_mask{i,l} = zeros(n,n);
        region_mask{i,l}(:,node_idx_modules{i,l}) = CommonNetwork(:,node_idx_modules{i,l});
        region_mask{i,l}(node_idx_modules{i,l},:) = CommonNetwork(node_idx_modules{i,l},:);
        region_mask{i,l}(find(region_mask{i,l}(:,:))) = 1;
    end
end
% save('results\region_mask.mat','region_mask');
% save('results\node_idx.mat','node_idx');
% save('results\node_idx_modules.mat','node_idx_modules');

u_vector = cell(nhub,select_modules);
for i = 1:nhub
    for l = 1:select_modules
        u_vector{i,l} = zeros(n,1);
        for j = 1:size(node_idx_modules{i,l})
            u_vector{i,l}(node_idx_modules{i,l}(j),1) = 1;
        end
    end
end
% save('results\u_vector.mat','u_vector');
end