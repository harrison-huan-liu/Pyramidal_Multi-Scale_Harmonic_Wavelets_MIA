function [select_scale,select_node_index,select_scale_index,select_node_index_index,select_node_index_index_name,node_scale,node_index_scale,multiPhi_scale] = select_node(SVM_para,hub_nodes,multiPhi,u_vector,CommonNetwork,dir_name)
nhub = size(hub_nodes,1);
destriux_148 = readcell('destriux_148.xlsx');

Tau_SVM_weight_multi = cell(1,nhub);
for i = 1:nhub
    Tau_SVM_weight_multi{1,i} = SVM_para(1+18*(i-1):18*i,1:10);
end
filename_svmweight = ['results_',dir_name,'\SVM_weight_lc0620.xls'];
% for writexlsx = 1:nhub
%     writematrix(Tau_SVM_weight_multi{1,writexlsx}.',filename_svmweight,'WriteMode','append')
% end

% select wavelets
SVM_para_ave = mean(SVM_para,2);
[SVM_para_order,SVM_para_index] = sort(SVM_para_ave);
SVM_para_scale1 = zeros(nhub*6,1);
SVM_para_scale2 = zeros(nhub*6,1);
SVM_para_scale3 = zeros(nhub*6,1);
for i = 1:nhub
    SVM_para_scale1(1+6*(i-1):6*i,1) = SVM_para_ave(1+18*(i-1):6+18*(i-1),1);
    SVM_para_scale2(1+6*(i-1):6*i,1) = SVM_para_ave(7+18*(i-1):12+18*(i-1),1);
    SVM_para_scale3(1+6*(i-1):6*i,1) = SVM_para_ave(13+18*(i-1):18+18*(i-1),1);
end

%% new_ex
% SVM_para_scale1_ave = zeros(nhub,1);
% SVM_para_scale2_ave = zeros(nhub,1);
% SVM_para_scale3_ave = zeros(nhub,1);
% for i = 1:nhub
%     SVM_para_scale1_ave(i) = mean(SVM_para_scale1(1+6*(i-1):6*i,1));
%     SVM_para_scale2_ave(i) = mean(SVM_para_scale2(1+6*(i-1):6*i,1));
%     SVM_para_scale3_ave(i) = mean(SVM_para_scale3(1+6*(i-1):6*i,1));
% end

%% end of the new_ex

[SVM_para_scale1_order,SVM_para_scale1_index] = sort(SVM_para_scale1,'descend');
[SVM_para_scale2_order,SVM_para_scale2_index] = sort(SVM_para_scale2,'descend');
[SVM_para_scale3_order,SVM_para_scale3_index] = sort(SVM_para_scale3,'descend');
multiPhi_scale1 = zeros(148,10);
multiPhi_scale2 = zeros(148,10);
multiPhi_scale3 = zeros(148,10);
node_scale1 = cell(1,10);
node_scale2 = cell(1,10);
node_scale3 = cell(1,10);
node_index_scale1 = zeros(1,10);
node_index_scale2 = zeros(1,10);
node_index_scale3 = zeros(1,10);
for i = 1:10
    node_index_scale1(1,i) = hub_nodes(floor((SVM_para_scale1_index(i)-1)/6) + 1);
    node_scale1{1,i} = destriux_148{node_index_scale1(1,i),6};
    node_index_scale2(1,i) = hub_nodes(floor((SVM_para_scale2_index(i)-1)/6) + 1);
    node_scale2{1,i} = destriux_148{node_index_scale2(1,i),6};
    node_index_scale3(1,i) = hub_nodes(floor((SVM_para_scale3_index(i)-1)/6) + 1);
    node_scale3{1,i} = destriux_148{node_index_scale3(1,i),6};
    multiPhi_scale1(:,i) = multiPhi(:,1+mod(SVM_para_scale1_index(i)-1,6)+floor((SVM_para_scale1_index(i)-1)/6)*18).*u_vector{floor((SVM_para_scale1_index(i)-1)/6) + 1,1};
    multiPhi_scale2(:,i) = multiPhi(:,7+mod(SVM_para_scale2_index(i)-1,6)+floor((SVM_para_scale2_index(i)-1)/6)*18).*u_vector{floor((SVM_para_scale2_index(i)-1)/6) + 1,2};
    multiPhi_scale3(:,i) = multiPhi(:,13+mod(SVM_para_scale3_index(i)-1,6)+floor((SVM_para_scale3_index(i)-1)/6)*18).*u_vector{floor((SVM_para_scale3_index(i)-1)/6) + 1,3};
end
node_scale = [node_scale1,node_scale2,node_scale3];
node_index_scale = [node_index_scale1,node_index_scale2,node_index_scale3];
multiPhi_scale = [multiPhi_scale1,multiPhi_scale2,multiPhi_scale3];
fnode = cell(30,1);
fnode_index = cell(30,1);
fedge = cell(30,1);
fnode_n = zeros(1,30);
for i = 1:30
    fnode_n(1,i) = length(find(abs(multiPhi_scale(:,i))>0));
    fnode_index{i,1} = find(abs(multiPhi_scale(:,i))>0);
    fnode{i,1} = destriux_148(fnode_index{i,1},:);
    for j = 1:fnode_n(1,i)
        fnode{i,1}{j,1} = fnode{i,1}{j,1}/1.1;
        fnode{i,1}{j,2} = fnode{i,1}{j,2}/1.1;
        fnode{i,1}{j,3} = fnode{i,1}{j,3}/1.1;
        fnode{i,1}{j,6} = multiPhi_scale(fnode_index{i,1}(j),i);
        if node_index_scale(1,i) == fnode_index{i,1}(j)
            fnode{i,1}{j,4} = 2;
            fnode{i,1}{j,5} = 4;
        end
    end
    
    filename_fnode = ['results_',dir_name,'\fnode_',num2str(i),'.txt'];
    writecell(fnode{i,1},filename_fnode,'Delimiter','tab');
    
    fedge{i,1} = zeros(fnode_n(1,i),fnode_n(1,i));
    for j = 1:fnode_n(1,i)
        for k = 1:fnode_n(1,i)
            fedge{i,1}(j,k) = CommonNetwork(fnode_index{i,1}(j),fnode_index{i,1}(k));
        end
    end
    
    filename_fedge = ['results_',dir_name,'\fedge_',num2str(i),'.txt'];
    writematrix(fedge{i,1},filename_fedge,'Delimiter','tab');
end

select_scale = cell(1,nhub);
select_scale1 = cell(1,nhub);
select_scale2 = cell(1,nhub);
select_scale3 = cell(1,nhub);
select_m_scale1 = zeros(1,nhub);
select_m_scale2 = zeros(1,nhub);
select_m_scale3 = zeros(1,nhub);
select_wavelets1 = zeros(1,nhub);
select_wavelets2 = zeros(1,nhub);
select_wavelets3 = zeros(1,nhub);
for i = 1:nhub
    select_scale{1,i} = mean(SVM_para(1+18*(i-1):18*i,1:10),2);
    select_scale1{1,i} = select_scale{1,i}(1:6,1);
    select_scale2{1,i} = select_scale{1,i}(7:12,1);
    select_scale3{1,i} = select_scale{1,i}(13:18,1);
%     select_m_scale1(1,i) = mean(select_scale1{1,i});
%     select_m_scale2(1,i) = mean(select_scale2{1,i});
%     select_m_scale3(1,i) = mean(select_scale3{1,i});
    [select_m_scale1(1,i),select_wavelets1(1,i)] = max(select_scale1{1,i});
    [select_m_scale2(1,i),select_wavelets2(1,i)] = max(select_scale2{1,i});
    [select_m_scale3(1,i),select_wavelets3(1,i)] = max(select_scale3{1,i});
end
% for writexlsx = 1:nhub
%     writematrix(select_scale{1,writexlsx}.','results\select_scale.xls','WriteMode','append')
% end
select_m_scale = [select_m_scale1,select_m_scale2,select_m_scale3];
[weight_order,select_index] = sort(select_m_scale,'descend');
node_index = mod(select_index-1,40)+1;
scale_index = floor((select_index-1)./40)+1;
select_node_index_num = 1;
select_node_index_repeat = 0;
select_node_index = [];
select_scale_index = [];
scale1 = 0;
scale2 = 0;
scale3 = 0;
while select_node_index_num<16
    if ismember(node_index(select_node_index_repeat+select_node_index_num),select_node_index)
        select_node_index_repeat = select_node_index_repeat + 1;
    else
        switch scale_index(select_node_index_repeat+select_node_index_num)
            case 1
                scale1 = scale1 + 1;
            case 2
                scale2 = scale2 + 1;
            case 3
                scale3 = scale3 + 1;
        end
        if scale1<6&&scale2<6&&scale3<6
            select_node_index = [select_node_index,node_index(select_node_index_repeat+select_node_index_num)];
            select_scale_index = [select_scale_index,scale_index(select_node_index_repeat+select_node_index_num)];
            select_node_index_num = select_node_index_num + 1;
        else
            select_node_index_repeat = select_node_index_repeat + 1;
            if scale1==6
                scale1 = scale1 - 1;
            elseif scale2==6
                scale2 = scale2 - 1;
            else
                scale3 = scale3 - 1;
            end
        end
    end
end
select_node_index_index = hub_nodes(select_node_index);

select_node_index_index_name = cell(16,1);
for i = 1:15
    select_node_index_index_name{i,1} = destriux_148{select_node_index_index(i),6};
end

multiPhi_scale_select = zeros(148,9);
% multiPhi_scale_order = [1,3,6,11,14,15,21,23,24]; %amy
multiPhi_scale_order = [1,3,5,11,12,15,22,25,21]; %tau
for i = 1:9
    multiPhi_scale_select(:,i) = multiPhi_scale(:,multiPhi_scale_order(i));
end
multiPhi_scale_select_ce = [multiPhi_scale1(:,1:3),multiPhi_scale2(:,1:3),multiPhi_scale3(:,1:3)];
filename5 = ['results_',dir_name,'\multiPhi_scale_select_ce.mat'];
% fprintf('Scale 1, Amyloid: %d; \n',SVM_para_scale1_index);
% fprintf('Scale 2, Amyloid: %d; \n',SVM_para_scale2_index);
% fprintf('Scale 3, Amyloid: %d; \n',SVM_para_scale3_index);
filename6 = ['results_',dir_name,'\SVM_para_scale1_index_ce.mat'];
filename7 = ['results_',dir_name,'\SVM_para_scale2_index_ce.mat'];
filename8 = ['results_',dir_name,'\SVM_para_scale3_index_ce.mat'];
save(filename5,'multiPhi_scale_select_ce');
save(filename6,'SVM_para_scale1_index');
save(filename7,'SVM_para_scale2_index');
save(filename8,'SVM_para_scale3_index');

% filename1 = ['results_',dir_name,'\select_node_index_index.mat'];
% save(filename1,'select_node_index_index');
% filename2 = ['results_',dir_name,'\select_node_index_index_name.mat'];
% save(filename2,'select_node_index_index_name');
% filename3 = ['results_',dir_name,'\node_scale.mat'];
% save(filename3,'node_scale');
% filename4 = ['results_',dir_name,'\node_index_scale.mat'];
% save(filename4,'node_index_scale');
% filename5 = ['results_',dir_name,'\multiPhi_scale.mat'];
% save(filename5,'multiPhi_scale');
end