close all;clear;clc
%% synthetic data (int)
node_num = 148;
net_num = 100;
repeat_num = 100;
node_num_max = 20; % 60
scale_max = 3;
sample_num = 200;
consider_node = 5:5:node_num_max;
sigma_error = 0.005:0.005:0.05;
MST = cell(repeat_num,net_num+1);
brain_network = cell(repeat_num,net_num+1);
node_degree = cell(repeat_num,1);
wavelet_num = cell(repeat_num,1);
Phi_ave = cell(repeat_num,1);
Phi_ave_wavelet = cell(repeat_num,1);

%% set wavelet
AD_wavelet = cell(repeat_num,1);
CN_wavelet = cell(repeat_num,1);
node_max = cell(repeat_num,1);
node_max_index = cell(repeat_num,1);
node_min = cell(repeat_num,1);
node_min_index = cell(repeat_num,1);
rand_max = cell(repeat_num,1);
rand_min = cell(repeat_num,1);
node_max_index_new = cell(repeat_num,1);
node_min_index_new = cell(repeat_num,1);
scale_num = cell(repeat_num,1);
set_wavelet = cell(repeat_num,node_num_max);
node_idx_dim_ad_after = cell(repeat_num,1);
node_idx_dim_cn_after = cell(repeat_num,1);

%% create signal data
CN_signdata = cell(repeat_num,scale_max*size(consider_node,2)*size(sigma_error,2));
AD_signdata = cell(repeat_num,scale_max*size(consider_node,2)*size(sigma_error,2));
for j = 1:scale_max*size(consider_node,2)*size(sigma_error,2)
    for i = 1:repeat_num
        CN_signdata{i,j} = zeros(node_num,sample_num);
        AD_signdata{i,j} = zeros(node_num,sample_num);
    end
end
weight_vector_cn = cell(repeat_num,sample_num);
weight_vector_ad = cell(repeat_num,sample_num);
sample_cn = cell(repeat_num,sample_num);
sample_ad = cell(repeat_num,sample_num);

multiwavelet = cell(repeat_num,1);
singlewavelet = cell(repeat_num,1);
node_idx_dim_ad_after_single = cell(repeat_num,1);
node_idx_dim_cn_after_single = cell(repeat_num,1);
AD_wavelet_single = cell(repeat_num,1);
CN_wavelet_single = cell(repeat_num,1);

multiwavelet_AD = cell(repeat_num,1);
multiwavelet_CN = cell(repeat_num,1);

ave_right_num = zeros(repeat_num,scale_max*size(consider_node,2));
SVM_model = cell(repeat_num,scale_max*size(consider_node,2));
ave_right_num_multi = zeros(repeat_num,scale_max*size(consider_node,2));
SVM_model_multi = cell(repeat_num,scale_max*size(consider_node,2));
ave_right_num_single = zeros(repeat_num,scale_max*size(consider_node,2));
SVM_model_single = cell(repeat_num,scale_max*size(consider_node,2));
ave_right_num_global = zeros(repeat_num,scale_max*size(consider_node,2));
SVM_model_global = cell(repeat_num,scale_max*size(consider_node,2));

%% repeat experiment
for repeat = 1:repeat_num
    random_network = cell(net_num,1);
    half_network = cell(net_num,1);
    
    %% calculate MST
    for brainnet_order = 1:net_num
        random_network{brainnet_order} = 0.02 .* abs(randn(node_num));
        half_network{brainnet_order} = tril(random_network{brainnet_order},-1);
        brain_network{repeat,brainnet_order} = half_network{brainnet_order} + half_network{brainnet_order}.';
        i = 1;
        MST{repeat,brainnet_order} = zeros(node_num);
        while i <= node_num
            [maximum_column,index_column] = max(brain_network{repeat,brainnet_order},[],1);
            [maximun,index] = max(maximum_column);
            brain_network{repeat,brainnet_order}(index_column(index),index) = 0;
            MST{repeat,brainnet_order}(index_column(index),index) = maximun;
            MST{repeat,brainnet_order}(index,index_column(index)) = maximun;
            t = MST{repeat,brainnet_order};
            for k = 2:size(MST{repeat,brainnet_order}, 1)
                t = t | (t*MST{repeat,brainnet_order});
            end

            if any(diag(t))
                acc = 1;
            else
                acc = 0;
            end

            if acc == 0 % the reason of 0
                MST{repeat,brainnet_order}(index_column(index),index) = 0;
                MST{repeat,brainnet_order}(index,index_column(index)) = 0;
                i = i-1;
            end
            i = i+1;
        end
    end

    brain_network{repeat,net_num+1} = zeros(node_num);
    MST{repeat,net_num+1} = zeros(node_num);
    for brainnet_order = 1:net_num
        MST{repeat,net_num+1} = MST{repeat,net_num+1} + MST{repeat,brainnet_order};
        brain_network{repeat,net_num+1} = brain_network{repeat,net_num+1} + brain_network{repeat,brainnet_order};
    end
    MST{repeat,net_num+1} = MST{repeat,net_num+1}./net_num;
    MST{repeat,net_num+1}(MST{repeat,net_num+1}<0.001) = 0;
    brain_network{repeat,net_num+1} = brain_network{repeat,net_num+1}./net_num;

    node_degree{repeat,1} = degrees_und(MST{repeat,net_num+1});
    wavelet_num{repeat,1} = round(mean(node_degree{repeat,1},2));
    
    %% calculate harmonic
    temp_D = diag(sum(brain_network{repeat,net_num+1},2));
    LatentLaplacian = temp_D - brain_network{repeat,net_num+1};
    [Phi_temp,~] = eig(LatentLaplacian);
    Phi_ave{repeat,1} = Phi_temp(:,1:60);
    Phi_ave_wavelet{repeat,1} = Phi_ave{repeat,1}(:,1:wavelet_num{repeat,1});

    %% 根据度选择与AD和CN相关的节点和小波
    [node_max{repeat,1},node_max_index{repeat,1}] = maxk(node_degree{repeat,1},node_num_max);
    [node_min{repeat,1},node_min_index{repeat,1}] = mink(node_degree{repeat,1},node_num_max);
    rand_max{repeat,1} = randperm(size(node_max_index{repeat,1},2));
    rand_min{repeat,1} = randperm(size(node_min_index{repeat,1},2));
    node_max_index_new{repeat,1} = node_max_index{repeat,1}(rand_max{repeat,1});
    node_min_index_new{repeat,1} = node_min_index{repeat,1}(rand_min{repeat,1});

    %% calculate wavelet(set ground truth)
    AD_wavelet{repeat,1} = zeros(node_num,node_num_max*scale_max*round(wavelet_num{repeat,1}/scale_max));
    CN_wavelet{repeat,1} = zeros(node_num,node_num_max*scale_max*round(wavelet_num{repeat,1}/scale_max));
    for node_order = 1:node_num_max
        node_idx_dim_ad_before = node_min_index_new{repeat,1}(node_order);
        node_idx_dim_cn_before = node_max_index_new{repeat,1}(node_order);
        for scale = 1:scale_max
            scale_num{repeat,1}(scale,1) = round(wavelet_num{repeat,1}/scale_max);
            set_wavelet{repeat,node_order}(scale,:) = randperm(wavelet_num{repeat,1},2*scale_num{repeat,1}(scale));
            for cal_scale = 1:scale
                for i = 1:size(node_idx_dim_ad_before)
                    node_idx_dim_ad_add = find(MST{repeat,net_num+1}(:,node_idx_dim_ad_before(i)));
                    node_idx_dim_ad_after{repeat} = union(node_idx_dim_ad_before,node_idx_dim_ad_add);
                end
                for i = 1:size(node_idx_dim_cn_before)
                    node_idx_dim_cn_add = find(MST{repeat,net_num+1}(:,node_idx_dim_cn_before(i)));
                    node_idx_dim_cn_after{repeat} = union(node_idx_dim_cn_before,node_idx_dim_cn_add);
                end
                node_idx_dim_ad_before = node_idx_dim_ad_after{repeat};
                node_idx_dim_cn_before = node_idx_dim_cn_after{repeat};
            end
            start_num = scale_max*round(wavelet_num{repeat,1}/scale_max)*(node_order-1)+round(wavelet_num{repeat,1}/scale_max)*(scale-1)+1;
            end_num = scale_max*round(wavelet_num{repeat,1}/scale_max)*(node_order-1)+round(wavelet_num{repeat,1}/scale_max)*scale;
            for i = 1:size(node_idx_dim_ad_after{repeat})
                AD_wavelet{repeat,1}(node_idx_dim_ad_after{repeat}(i),start_num:end_num) = Phi_ave_wavelet{repeat,1}(node_idx_dim_ad_after{repeat}(i),set_wavelet{repeat,node_order}(scale,1:scale_num{repeat,1}(scale)));
            end
            for i = 1:size(node_idx_dim_cn_after{repeat})
                CN_wavelet{repeat,1}(node_idx_dim_cn_after{repeat}(i),start_num:end_num) = Phi_ave_wavelet{repeat,1}(node_idx_dim_cn_after{repeat}(i),set_wavelet{repeat,node_order}(scale,scale_num{repeat,1}(scale)+1:2*scale_num{repeat,1}(scale)));
            end
        end
    end

    %% calculate signal data
    for sample = 1:sample_num
        weight_vector_cn{repeat,sample} = abs(randn(1,node_num_max*scale_max*round(wavelet_num{repeat,1}/scale_max)));
        weight_vector_ad{repeat,sample} = abs(randn(1,node_num_max*scale_max*round(wavelet_num{repeat,1}/scale_max)));
        sample_cn{repeat,sample} = weight_vector_cn{repeat,sample};
        sample_ad{repeat,sample} = weight_vector_ad{repeat,sample};
        
        %% consider the different scale
        for scale = 1:scale_max
            for sigma_error_count = 1:size(sigma_error,2)
                for conside_node_count = 1:size(consider_node,2)
                    loop_node = scale_max*consider_node(conside_node_count)*round(wavelet_num{repeat,1}/scale_max);
                    wavelet_order = scale_max*conside_node_count-scale+1+(scale_max*size(consider_node,2))*(sigma_error_count-1);
                    for i = 1:loop_node
                        CN_signdata{repeat,wavelet_order}(:,sample) = CN_signdata{repeat,wavelet_order}(:,sample) + weight_vector_cn{repeat,sample}(1,i) .* CN_wavelet{repeat,1}(:,i);
                        AD_signdata{repeat,wavelet_order}(:,sample) = AD_signdata{repeat,wavelet_order}(:,sample) + weight_vector_ad{repeat,sample}(1,i) .* AD_wavelet{repeat,1}(:,i);
                    end
                    CN_signdata{repeat,wavelet_order}(:,sample) = CN_signdata{repeat,wavelet_order}(:,sample)/sum(weight_vector_cn{repeat,sample}(1,1:loop_node)) + sigma_error(sigma_error_count) .* randn(node_num,1);
                    AD_signdata{repeat,wavelet_order}(:,sample) = AD_signdata{repeat,wavelet_order}(:,sample)/sum(weight_vector_cn{repeat,sample}(1,1:loop_node)) + sigma_error(sigma_error_count) .* randn(node_num,1);
                end
            end
            for i = 1:node_num_max
                start_wavelet_num = (scale_max-scale)*round(wavelet_num{repeat,1}/scale_max)+(i-1)*scale_max*round(wavelet_num{repeat,1}/scale_max)+1;
                end_wavelet_num = (scale_max-scale)*round(wavelet_num{repeat,1}/scale_max)+(i-1)*scale_max*round(wavelet_num{repeat,1}/scale_max)+round(wavelet_num{repeat,1}/scale_max);
                weight_vector_cn{repeat,sample}(start_wavelet_num:end_wavelet_num)=0;
                weight_vector_ad{repeat,sample}(start_wavelet_num:end_wavelet_num)=0;
            end
        end
    end

    %% SVM
% for repeat = 1:63
    for i = 1:scale_max*size(consider_node,2)*size(sigma_error,2)
        [right_num,SVM_model{repeat,i}] = SVM_signaldata(CN_signdata{repeat,i},AD_signdata{repeat,i},sample_num);
        ave_right_num(repeat,i) = mean(right_num);
    end
    
    multiwavelet_AD{repeat,1} = zeros(node_num,node_num_max*scale_max*wavelet_num{repeat,1});
    multiwavelet_CN{repeat,1} = zeros(node_num,node_num_max*scale_max*wavelet_num{repeat,1});
    for node_order = 1:node_num_max
        node_idx_dim_ad_before = node_min_index_new{repeat,1}(node_order);
        node_idx_dim_cn_before = node_max_index_new{repeat,1}(node_order);
        for scale = 1:scale_max
%             scale_num{repeat,1}(scale,1) = round(wavelet_num{repeat,1}/scale_max);
%             set_wavelet{repeat,node_order}(scale,:) = randperm(wavelet_num{repeat,1},2*scale_num{repeat,1}(scale));
            for cal_scale = 1:scale
                for i = 1:size(node_idx_dim_ad_before)
                    node_idx_dim_ad_add = find(MST{repeat,net_num+1}(:,node_idx_dim_ad_before(i)));
                    node_idx_dim_ad_after{repeat} = union(node_idx_dim_ad_before,node_idx_dim_ad_add);
                end
                for i = 1:size(node_idx_dim_cn_before)
                    node_idx_dim_cn_add = find(MST{repeat,net_num+1}(:,node_idx_dim_cn_before(i)));
                    node_idx_dim_cn_after{repeat} = union(node_idx_dim_cn_before,node_idx_dim_cn_add);
                end
                node_idx_dim_ad_before = node_idx_dim_ad_after{repeat};
                node_idx_dim_cn_before = node_idx_dim_cn_after{repeat};
            end
            start_num = scale_max*wavelet_num{repeat,1}*(node_order-1)+wavelet_num{repeat,1}*(scale-1)+1;
            end_num = scale_max*wavelet_num{repeat,1}*(node_order-1)+wavelet_num{repeat,1}*scale;
            for i = 1:size(node_idx_dim_ad_after{repeat})
                multiwavelet_AD{repeat,1}(node_idx_dim_ad_after{repeat}(i),start_num:end_num) = Phi_ave_wavelet{repeat,1}(node_idx_dim_ad_after{repeat}(i),1:wavelet_num{repeat,1});
            end
            for i = 1:size(node_idx_dim_cn_after{repeat})
                multiwavelet_CN{repeat,1}(node_idx_dim_cn_after{repeat}(i),start_num:end_num) = Phi_ave_wavelet{repeat,1}(node_idx_dim_cn_after{repeat}(i),1:wavelet_num{repeat,1});
            end
        end
    end
    multiwavelet{repeat,1} = [multiwavelet_AD{repeat,1} multiwavelet_CN{repeat,1}];
    for i = 1:scale_max*size(consider_node,2)*size(sigma_error,2)
        [right_num_multi,SVM_model_multi{repeat,i}] = SVM_signaldata(multiwavelet{repeat,1}.'*CN_signdata{repeat,i},multiwavelet{repeat,1}.'*AD_signdata{repeat,i},sample_num);
        ave_right_num_multi(repeat,i) = mean(right_num_multi);
    end
    
    AD_wavelet_single{repeat,1} = zeros(node_num,node_num_max*wavelet_num{repeat,1});
    CN_wavelet_single{repeat,1} = zeros(node_num,node_num_max*wavelet_num{repeat,1});
    for node_order = 1:node_num_max
        node_idx_dim_ad_before_single = node_min_index_new{repeat,1}(node_order);
        node_idx_dim_cn_before_single = node_max_index_new{repeat,1}(node_order);
        for i = 1:size(node_idx_dim_ad_before_single)
            node_idx_dim_ad_add_single = find(MST{repeat,net_num+1}(:,node_idx_dim_ad_before_single(i)));
            node_idx_dim_ad_after_single{repeat} = union(node_idx_dim_ad_before_single,node_idx_dim_ad_add_single);
        end
        for i = 1:size(node_idx_dim_cn_before_single)
            node_idx_dim_cn_add_single = find(MST{repeat,net_num+1}(:,node_idx_dim_cn_before_single(i)));
            node_idx_dim_cn_after_single{repeat} = union(node_idx_dim_cn_before_single,node_idx_dim_cn_add_single);
        end
        start_num_single = wavelet_num{repeat,1}*(node_order-1)+1;
        end_num_single = wavelet_num{repeat,1}*(node_order-1)+wavelet_num{repeat,1};
        for i = 1:size(node_idx_dim_ad_after_single{repeat})
            AD_wavelet_single{repeat,1}(node_idx_dim_ad_after_single{repeat}(i),start_num_single:end_num_single) = Phi_ave_wavelet{repeat,1}(node_idx_dim_ad_after_single{repeat}(i),:);
        end
        for i = 1:size(node_idx_dim_cn_after_single{repeat})
            CN_wavelet_single{repeat,1}(node_idx_dim_cn_after_single{repeat}(i),start_num_single:end_num_single) = Phi_ave_wavelet{repeat,1}(node_idx_dim_cn_after_single{repeat}(i),:);
        end
    end
    singlewavelet{repeat,1} = [AD_wavelet_single{repeat,1} CN_wavelet_single{repeat,1}];
    for i = 1:scale_max*size(consider_node,2)*size(sigma_error,2)
        [right_num_single,SVM_model_single{repeat,i}] = SVM_signaldata(singlewavelet{repeat,1}.'*CN_signdata{repeat,i},singlewavelet{repeat,1}.'*AD_signdata{repeat,i},sample_num);
        ave_right_num_single(repeat,i) = mean(right_num_single);
    end
    
    for i = 1:scale_max*size(consider_node,2)*size(sigma_error,2)
        [right_num_global,SVM_model_global{repeat,i}] = SVM_signaldata(Phi_ave{repeat,1}.'*CN_signdata{repeat,i},Phi_ave{repeat,1}.'*AD_signdata{repeat,i},sample_num);
        ave_right_num_global(repeat,i) = mean(right_num_global);
    end
    
    fprintf('Now the code is running in the %d loop; \n',repeat);
end
final_acc = mean(ave_right_num(1:repeat,:));
final_acc_var = std(ave_right_num(1:repeat,:));
final_acc_multi = mean(ave_right_num_multi(1:repeat,:));
final_acc_multi_var = std(ave_right_num_multi(1:repeat,:));
final_acc_single = mean(ave_right_num_single(1:repeat,:));
final_acc_single_var = std(ave_right_num_single(1:repeat,:));
final_acc_global = mean(ave_right_num_global(1:repeat,:));
final_acc_global_var = std(ave_right_num_global(1:repeat,:));

% ptest_index = zeros(10,1);
% for i = 1:10
%     ptest_index(i) = 36*(i-1)+27;
% end
% [p_acc_ms,~,~]=permutationTest(final_acc_multi(1,ptest_index(i)),final_acc_single(1,ptest_index(i)),10000);
% [p_acc_mg,~,~]=permutationTest(final_acc_multi(1,ptest_index(i)),final_acc_global(1,ptest_index(i)),10000);
% [p_acc_mo,~,~]=permutationTest(final_acc_multi(1,ptest_index(i)),final_acc(1,ptest_index(i)),10000);
% 
% p_acc_ms = zeros(10,1);
% p_acc_mg = zeros(10,1);
% p_acc_mo = zeros(10,1);
% for i = 1:10
%     [p_acc_ms(i),~,~]=permutationTest(ave_right_num_multi(1:repeat,ptest_index(i)),ave_right_num_single(1:repeat,ptest_index(i)),10000);
%     [p_acc_mg(i),~,~]=permutationTest(ave_right_num_multi(1:repeat,ptest_index(i)),ave_right_num_global(1:repeat,ptest_index(i)),10000);
%     [p_acc_mo(i),~,~]=permutationTest(ave_right_num_multi(1:repeat,ptest_index(i)),ave_right_num(1:repeat,ptest_index(i)),10000);
% end

save('results\brain_network.mat','brain_network');
save('results\MST.mat','MST');
save('results\node_degree.mat','node_degree');
save('results\Phi_ave.mat','Phi_ave');
save('results\Phi_ave_wavelet.mat','Phi_ave_wavelet');
save('results\wavelet_num.mat','wavelet_num');
save('results\AD_wavelet.mat','AD_wavelet');
save('results\CN_wavelet.mat','CN_wavelet');
save('results\node_max.mat','node_max');
save('results\node_max_index.mat','node_max_index');
save('results\node_min.mat','node_min');
save('results\node_min_index.mat','node_min_index');
save('results\rand_max.mat','rand_max');
save('results\rand_min.mat','rand_min');
save('results\node_max_index_new.mat','node_max_index_new');
save('results\node_min_index_new.mat','node_min_index_new');
save('results\scale_num.mat','scale_num');
save('results\set_wavelet.mat','set_wavelet');
save('results\CN_signdata.mat','CN_signdata');
save('results\AD_signdata.mat','AD_signdata');
save('results\sample_cn.mat','sample_cn');
save('results\sample_ad.mat','sample_ad');
save('results\node_idx_dim_ad_after.mat','node_idx_dim_ad_after');
save('results\node_idx_dim_cn_after.mat','node_idx_dim_cn_after');
save('results\SVM_model.mat','SVM_model');
save('results\ave_right_num.mat','ave_right_num');
save('results\final_acc.mat','final_acc');
save('results\SVM_model_multi.mat','SVM_model_multi');
save('results\ave_right_num_multi.mat','ave_right_num_multi');
save('results\final_acc_multi.mat','final_acc_multi');
save('results\SVM_model_single.mat','SVM_model_single');
save('results\ave_right_num_single.mat','ave_right_num_single');
save('results\final_acc_single.mat','final_acc_single');
save('results\SVM_model_global.mat','SVM_model_global');
save('results\ave_right_num_global.mat','ave_right_num_global');
save('results\final_acc_global.mat','final_acc_global');
save('results\multiwavelet.mat','multiwavelet');
save('results\singlewavelet.mat','singlewavelet');
save('results\node_idx_dim_ad_after_single.mat','node_idx_dim_ad_after_single');
save('results\node_idx_dim_cn_after_single.mat','node_idx_dim_cn_after_single');
save('results\AD_wavelet_single.mat','AD_wavelet_single');
save('results\CN_wavelet_single.mat','CN_wavelet_single');

writematrix(final_acc(1:scale_max*size(consider_node,2)),'accuracy_original.xls','Sheet',1)
for writexlsx = 1:size(sigma_error,2)-1
    start_num_writexlsx = scale_max*size(consider_node,2)*writexlsx+1;
    end_num_writexlsx = scale_max*size(consider_node,2)*(writexlsx+1);
    writematrix(final_acc(start_num_writexlsx:end_num_writexlsx),'accuracy_original.xls','WriteMode','append')
end

writematrix(final_acc_global(1:scale_max*size(consider_node,2)),'accuracy_global.xls','Sheet',1)
for writexlsx = 1:size(sigma_error,2)-1
    start_num_writexlsx = scale_max*size(consider_node,2)*writexlsx+1;
    end_num_writexlsx = scale_max*size(consider_node,2)*(writexlsx+1);
    writematrix(final_acc_global(start_num_writexlsx:end_num_writexlsx),'accuracy_global.xls','WriteMode','append')
end

writematrix(final_acc_multi(1:scale_max*size(consider_node,2)),'accuracy_multi.xls','Sheet',1)
for writexlsx = 1:size(sigma_error,2)-1
    start_num_writexlsx = scale_max*size(consider_node,2)*writexlsx+1;
    end_num_writexlsx = scale_max*size(consider_node,2)*(writexlsx+1);
    writematrix(final_acc_multi(start_num_writexlsx:end_num_writexlsx),'accuracy_multi.xls','WriteMode','append')
end

writematrix(final_acc_single(1:scale_max*size(consider_node,2)),'accuracy_single.xls','Sheet',1)
for writexlsx = 1:size(sigma_error,2)-1
    start_num_writexlsx = scale_max*size(consider_node,2)*writexlsx+1;
    end_num_writexlsx = scale_max*size(consider_node,2)*(writexlsx+1);
    writematrix(final_acc_single(start_num_writexlsx:end_num_writexlsx),'accuracy_single.xls','WriteMode','append')
end

writematrix(final_acc_var(1:scale_max*size(consider_node,2)),'accuracy_original_var.xls','Sheet',1)
for writexlsx = 1:size(sigma_error,2)-1
    start_num_writexlsx = scale_max*size(consider_node,2)*writexlsx+1;
    end_num_writexlsx = scale_max*size(consider_node,2)*(writexlsx+1);
    writematrix(final_acc_var(start_num_writexlsx:end_num_writexlsx),'accuracy_original_var.xls','WriteMode','append')
end

writematrix(final_acc_global_var(1:scale_max*size(consider_node,2)),'accuracy_global_var.xls','Sheet',1)
for writexlsx = 1:size(sigma_error,2)-1
    start_num_writexlsx = scale_max*size(consider_node,2)*writexlsx+1;
    end_num_writexlsx = scale_max*size(consider_node,2)*(writexlsx+1);
    writematrix(final_acc_global_var(start_num_writexlsx:end_num_writexlsx),'accuracy_global_var.xls','WriteMode','append')
end

writematrix(final_acc_multi_var(1:scale_max*size(consider_node,2)),'accuracy_multi_var.xls','Sheet',1)
for writexlsx = 1:size(sigma_error,2)-1
    start_num_writexlsx = scale_max*size(consider_node,2)*writexlsx+1;
    end_num_writexlsx = scale_max*size(consider_node,2)*(writexlsx+1);
    writematrix(final_acc_multi_var(start_num_writexlsx:end_num_writexlsx),'accuracy_multi_var.xls','WriteMode','append')
end

writematrix(final_acc_single_var(1:scale_max*size(consider_node,2)),'accuracy_single_var.xls','Sheet',1)
for writexlsx = 1:size(sigma_error,2)-1
    start_num_writexlsx = scale_max*size(consider_node,2)*writexlsx+1;
    end_num_writexlsx = scale_max*size(consider_node,2)*(writexlsx+1);
    writematrix(final_acc_single_var(start_num_writexlsx:end_num_writexlsx),'accuracy_single_var.xls','WriteMode','append')
end

%% node acuracy
% load('results\noise_[0.005_0.005_0.05]_node_[5_5_40]_scale_3\AD_wavelet.mat')
% load('results\noise_[0.005_0.005_0.05]_node_[5_5_40]_scale_3\CN_wavelet.mat')
% load('results\noise_[0.005_0.005_0.05]_node_[5_5_40]_scale_3\AD_signdata.mat')
% load('results\noise_[0.005_0.005_0.05]_node_[5_5_40]_scale_3\CN_signdata.mat')
% load('results\noise_[0.005_0.005_0.05]_node_[5_5_40]_scale_3\MST.mat')
% load('results\noise_[0.005_0.005_0.05]_node_[5_5_40]_scale_3\wavelet_num.mat')
% load('results\noise_[0.005_0.005_0.05]_node_[5_5_40]_scale_3\set_wavelet.mat')
% load('results\noise_[0.005_0.005_0.05]_node_[5_5_40]_scale_3\Phi_ave_wavelet.mat')
% load('results\noise_[0.005_0.005_0.05]_node_[5_5_40]_scale_3\node_idx_dim_ad_after.mat')
% load('results\noise_[0.005_0.005_0.05]_node_[5_5_40]_scale_3\node_idx_dim_cn_after.mat')
% load('results\noise_[0.005_0.005_0.05]_node_[5_5_40]_scale_3\node_min_index_new.mat')
% load('results\noise_[0.005_0.005_0.05]_node_[5_5_40]_scale_3\node_max_index_new.mat')
% load('results\noise_[0.005_0.005_0.05]_node_[5_5_40]_scale_3\SVM_model.mat')
% load('results\noise_[0.005_0.005_0.05]_node_[5_5_40]_scale_3\SVM_model_multi.mat')
% load('results\noise_[0.005_0.005_0.05]_node_[5_5_40]_scale_3\SVM_model_single.mat')
% load('results\noise_[0.005_0.005_0.05]_node_[5_5_40]_scale_3\SVM_model_global.mat')
% load('results\noise_[0.005_0.005_0.05]_node_[5_5_40]_scale_3\node_idx_dim_ad_after.mat');
% load('results\noise_[0.005_0.005_0.05]_node_[5_5_40]_scale_3\node_idx_dim_cn_after.mat');

node_acc_right = cell(repeat_num,scale_max*size(consider_node,2)*size(sigma_error,2));
node_acc = zeros(repeat_num,scale_max*size(consider_node,2)*size(sigma_error,2));
SVM_node_index = cell(repeat_num,scale_max*size(consider_node,2)*size(sigma_error,2));

node_acc_right_multi = cell(repeat_num,scale_max*size(consider_node,2)*size(sigma_error,2));
node_acc_multi = zeros(repeat_num,scale_max*size(consider_node,2)*size(sigma_error,2));
SVM_node_index_multi = cell(repeat_num,scale_max*size(consider_node,2)*size(sigma_error,2));

node_acc_right_single = cell(repeat_num,scale_max*size(consider_node,2)*size(sigma_error,2));
node_acc_single = zeros(repeat_num,scale_max*size(consider_node,2)*size(sigma_error,2));
SVM_node_index_single = cell(repeat_num,scale_max*size(consider_node,2)*size(sigma_error,2));

node_right_index_new = cell(repeat_num,1);
for  i = 1:repeat_num
    for j = 1:scale_max*size(consider_node,2)*size(sigma_error,2)
        node_correspond_num = floor(mod(j-1,12)/3)+1;
        node_right_index_new{i,1} = union(node_max_index_new{i,1}(1,1:consider_node(1,node_correspond_num)),node_min_index_new{i,1}(1,1:consider_node(1,node_correspond_num)));
        SVM_model{i,j} = normalize(abs(SVM_model{i,j}),'range');
        weight_rank = sort(abs(SVM_model{i,mod(j-1,12)+1}),'descend');
        SVM_node_index{i,j} = find(abs(SVM_model{i,j})>weight_rank(2*consider_node(1,node_correspond_num)));
        node_acc_right{i,j} = ismember(node_right_index_new{i,1}',SVM_node_index{i,j});
        node_acc(i,j) = sum(node_acc_right{i,j})/(2*consider_node(1,node_correspond_num));
        
        SVM_model_multi{i,j} = normalize(abs(SVM_model_multi{i,j}),'range');
        weight_rank_multi = sort(abs(SVM_model_multi{i,mod(j-1,12)+1}),'descend');
        SVM_node_index_multi{i,j} = find(abs(SVM_model_multi{i,j})>weight_rank_multi(size(SVM_model_multi{i,j},1)/node_num_max*consider_node(1,node_correspond_num))); % 
        SVM_node_index_multi{i,j} = floor((SVM_node_index_multi{i,j}-1)/(size(SVM_model_multi{i,j},1)/(2*node_num_max)))+1;
        node_acc_right_multi{i,j} = ismember(node_right_index_new{i,1}',SVM_node_index_multi{i,j});
        node_acc_multi(i,j) = sum(node_acc_right_multi{i,j})/(2*consider_node(1,node_correspond_num));
        
        SVM_model_single{i,j} = normalize(abs(SVM_model_single{i,j}),'range');
        weight_rank_single = sort(abs(SVM_model_single{i,mod(j-1,12)+1}),'descend');
        SVM_node_index_single{i,j} = find(abs(SVM_model_single{i,j})>weight_rank_single(size(SVM_model_single{i,j},1)/node_num_max*consider_node(1,node_correspond_num)));
        SVM_node_index_single{i,j} = floor((SVM_node_index_single{i,j}-1)/(size(SVM_model_single{i,j},1)/(2*node_num_max)))+1;
        node_acc_right_single{i,j} = ismember(node_right_index_new{i,1}',SVM_node_index_single{i,j});
        node_acc_single(i,j) = sum(node_acc_right_single{i,j})/(2*consider_node(1,node_correspond_num));
    end
end

node_acc_ave = mean(node_acc);
node_acc_ave_var = std(node_acc);
writematrix(node_acc_ave(1:scale_max*size(consider_node,2)),'node_accuracy_original.xls','Sheet',1)
for writexlsx = 1:size(sigma_error,2)-1
    start_num_writexlsx = scale_max*size(consider_node,2)*writexlsx+1;
    end_num_writexlsx = scale_max*size(consider_node,2)*(writexlsx+1);
    writematrix(node_acc_ave(start_num_writexlsx:end_num_writexlsx),'node_accuracy_original.xls','WriteMode','append')
end
writematrix(node_acc_ave_var(1:scale_max*size(consider_node,2)),'node_accuracy_original_var.xls','Sheet',1)
for writexlsx = 1:size(sigma_error,2)-1
    start_num_writexlsx = scale_max*size(consider_node,2)*writexlsx+1;
    end_num_writexlsx = scale_max*size(consider_node,2)*(writexlsx+1);
    writematrix(node_acc_ave_var(start_num_writexlsx:end_num_writexlsx),'node_accuracy_original_var.xls','WriteMode','append')
end

node_acc_ave_multi = mean(node_acc_multi);
node_acc_ave_multi_var = std(node_acc_multi);
writematrix(node_acc_ave_multi(1:scale_max*size(consider_node,2)),'node_accuracy_multi.xls','Sheet',1)
for writexlsx = 1:size(sigma_error,2)-1
    start_num_writexlsx = scale_max*size(consider_node,2)*writexlsx+1;
    end_num_writexlsx = scale_max*size(consider_node,2)*(writexlsx+1);
    writematrix(node_acc_ave_multi(start_num_writexlsx:end_num_writexlsx),'node_accuracy_multi.xls','WriteMode','append')
end
writematrix(node_acc_ave_multi_var(1:scale_max*size(consider_node,2)),'node_accuracy_multi_var.xls','Sheet',1)
for writexlsx = 1:size(sigma_error,2)-1
    start_num_writexlsx = scale_max*size(consider_node,2)*writexlsx+1;
    end_num_writexlsx = scale_max*size(consider_node,2)*(writexlsx+1);
    writematrix(node_acc_ave_multi_var(start_num_writexlsx:end_num_writexlsx),'node_accuracy_multi_var.xls','WriteMode','append')
end

node_acc_ave_single = mean(node_acc_single);
node_acc_ave_single_var = std(node_acc_single);
writematrix(node_acc_ave_single(1:scale_max*size(consider_node,2)),'node_accuracy_single.xls','Sheet',1)
for writexlsx = 1:size(sigma_error,2)-1
    start_num_writexlsx = scale_max*size(consider_node,2)*writexlsx+1;
    end_num_writexlsx = scale_max*size(consider_node,2)*(writexlsx+1);
    writematrix(node_acc_ave_single(start_num_writexlsx:end_num_writexlsx),'node_accuracy_single.xls','WriteMode','append')
end
writematrix(node_acc_ave_single_var(1:scale_max*size(consider_node,2)),'node_accuracy_single_var.xls','Sheet',1)
for writexlsx = 1:size(sigma_error,2)-1
    start_num_writexlsx = scale_max*size(consider_node,2)*writexlsx+1;
    end_num_writexlsx = scale_max*size(consider_node,2)*(writexlsx+1);
    writematrix(node_acc_ave_single_var(start_num_writexlsx:end_num_writexlsx),'node_accuracy_single_var.xls','WriteMode','append')
end

%% wavelet accuracy
wavelet_acc_right = cell(repeat_num,scale_max*size(consider_node,2)*size(sigma_error,2));
wavelet_acc = zeros(repeat_num,scale_max*size(consider_node,2)*size(sigma_error,2));
SVM_wavelet_index = cell(repeat_num,scale_max*size(consider_node,2)*size(sigma_error,2));

wavelet_acc_right_multi = cell(repeat_num,scale_max*size(consider_node,2)*size(sigma_error,2));
wavelet_acc_multi = zeros(repeat_num,scale_max*size(consider_node,2)*size(sigma_error,2));
SVM_wavelet_index_multi = cell(repeat_num,scale_max*size(consider_node,2)*size(sigma_error,2));

wavelet_acc_right_single = cell(repeat_num,scale_max*size(consider_node,2)*size(sigma_error,2));
wavelet_acc_single = zeros(repeat_num,scale_max*size(consider_node,2)*size(sigma_error,2));
SVM_wavelet_index_single = cell(repeat_num,scale_max*size(consider_node,2)*size(sigma_error,2));

wavelet_right_index_new = cell(repeat_num,2);
right_wavelet{repeat,1} = [AD_wavelet{repeat,1} CN_wavelet{repeat,1}];
ddd_num = scale_max*round(wavelet_num{i,1}/scale_max)*node_num_max;
for  i = 1:repeat_num
    wavelet_right_index_new{i,1} = zeros(1,2*ddd_num);
    wavelet_right_index_new{i,2} = zeros(1,2*ddd_num/scale_max);
    for no_index = 1:scale_max*size(consider_node,2)*size(sigma_error,2)
        wavelet_correspond_num = floor(mod(no_index-1,12)/3)+1;
        scale_correspond_num = mod(mod(no_index-1,12),3)+1;
%         wavelet_right_index_new{i,1} = union(wavelet_max_index_new{i,1}(1,1:consider_node(1,wavelet_correspond_num)),wavelet_min_index_new{i,1}(1,1:consider_node(1,wavelet_correspond_num)));
        for k = 1:3
            for j = 1:node_num_max
                aaa_num = scale_max*round(wavelet_num{i,1}/scale_max)*(j-1)+round(wavelet_num{i,1}/scale_max)*(k-1)+1;
                bbb_num = scale_max*round(wavelet_num{i,1}/scale_max)*(j-1)+round(wavelet_num{i,1}/scale_max)*k;
                ccc_num = scale_max*wavelet_num{i,1}*(j-1)+wavelet_num{i,1}*(k-1);
                wavelet_right_index_new{i,1}(1,aaa_num:bbb_num) = set_wavelet{i,j}(k,1:round(wavelet_num{i,1}/scale_max))+ccc_num;
                wavelet_right_index_new{i,1}(1,aaa_num+ddd_num:bbb_num+ddd_num) = set_wavelet{i,j}(k,round(wavelet_num{i,1}/scale_max)+1:2*round(wavelet_num{i,1}/scale_max))+ccc_num+wavelet_num{i,1}*scale_max*node_num_max;
            end
        end
        
        for j = 1:node_num_max
            aaa_num = round(wavelet_num{i,1}/scale_max)*(j-1)+1;
            bbb_num = round(wavelet_num{i,1}/scale_max)*j;
            ccc_num = wavelet_num{i,1}*(j-1);
            wavelet_right_index_new{i,2}(1,aaa_num:bbb_num) = set_wavelet{i,j}(1,1:round(wavelet_num{i,1}/scale_max))+ccc_num;
            wavelet_right_index_new{i,2}(1,aaa_num+ddd_num:bbb_num+ddd_num) = set_wavelet{i,j}(1,round(wavelet_num{i,1}/scale_max)+1:2*round(wavelet_num{i,1}/scale_max))+ccc_num+wavelet_num{i,1}*node_num_max;
        end
        
%         for ad_cn = 1:scale_max*round(wavelet_num{i,1}/scale_max)*node_num_max
%             wavelet_right_index_new{i,1}(1,ad_cn+scale_max*round(wavelet_num{i,1}/scale_max)*node_num_max) = wavelet_right_index_new{i,1}(1,ad_cn)+wavelet_num{i,1}*scale_max*node_num_max;
%         end

        weight_rank_multi = sort(abs(SVM_model_multi{i,mod(no_index-1,12)+1}),'descend');
        SVM_wavelet_index_multi{i,no_index} = find(abs(SVM_model_multi{i,no_index})>weight_rank_multi(2*2*consider_node(1,wavelet_correspond_num)*scale_correspond_num*round(wavelet_num{i,1}/scale_max)));
        wavelet_acc_right_multi{i,no_index} = ismember(wavelet_right_index_new{i,1}',SVM_wavelet_index_multi{i,no_index});
        wavelet_acc_multi(i,no_index) = sum(wavelet_acc_right_multi{i,no_index})/(2*consider_node(1,wavelet_correspond_num)*scale_correspond_num*round(wavelet_num{i,1}/scale_max));

        weight_rank_single = sort(abs(SVM_model_single{i,mod(no_index-1,12)+1}),'descend');
        SVM_wavelet_index_single{i,no_index} = find(abs(SVM_model_single{i,no_index})>weight_rank_single(2*2*consider_node(1,wavelet_correspond_num)*round(wavelet_num{i,1}/scale_max)));
        wavelet_acc_right_single{i,no_index} = ismember(wavelet_right_index_new{i,2}',SVM_wavelet_index_single{i,no_index});
        wavelet_acc_single(i,no_index) = sum(wavelet_acc_right_single{i,no_index})/(2*consider_node(1,wavelet_correspond_num)*scale_correspond_num*round(wavelet_num{i,1}/scale_max));
    end
end

wavelet_acc_ave_multi = mean(wavelet_acc_multi);
wavelet_acc_ave_multi_var = std(wavelet_acc_multi);
writematrix(wavelet_acc_ave_multi(1:scale_max*size(consider_node,2)),'wavelet_accuracy_multi.xls','Sheet',1)
for writexlsx = 1:size(sigma_error,2)-1
    start_num_writexlsx = scale_max*size(consider_node,2)*writexlsx+1;
    end_num_writexlsx = scale_max*size(consider_node,2)*(writexlsx+1);
    writematrix(wavelet_acc_ave_multi(start_num_writexlsx:end_num_writexlsx),'wavelet_accuracy_multi.xls','WriteMode','append')
end
writematrix(wavelet_acc_ave_multi_var(1:scale_max*size(consider_node,2)),'wavelet_accuracy_multi_var.xls','Sheet',1)
for writexlsx = 1:size(sigma_error,2)-1
    start_num_writexlsx = scale_max*size(consider_node,2)*writexlsx+1;
    end_num_writexlsx = scale_max*size(consider_node,2)*(writexlsx+1);
    writematrix(wavelet_acc_ave_multi_var(start_num_writexlsx:end_num_writexlsx),'wavelet_accuracy_multi_var.xls','WriteMode','append')
end

wavelet_acc_ave_single = mean(wavelet_acc_single);
wavelet_acc_ave_single_var = std(wavelet_acc_single);
writematrix(wavelet_acc_ave_single(1:scale_max*size(consider_node,2)),'wavelet_accuracy_single.xls','Sheet',1)
for writexlsx = 1:size(sigma_error,2)-1
    start_num_writexlsx = scale_max*size(consider_node,2)*writexlsx+1;
    end_num_writexlsx = scale_max*size(consider_node,2)*(writexlsx+1);
    writematrix(wavelet_acc_ave_single(start_num_writexlsx:end_num_writexlsx),'wavelet_accuracy_single.xls','WriteMode','append')
end
writematrix(wavelet_acc_ave_single_var(1:scale_max*size(consider_node,2)),'wavelet_accuracy_single_var.xls','Sheet',1)
for writexlsx = 1:size(sigma_error,2)-1
    start_num_writexlsx = scale_max*size(consider_node,2)*writexlsx+1;
    end_num_writexlsx = scale_max*size(consider_node,2)*(writexlsx+1);
    writematrix(wavelet_acc_ave_single_var(start_num_writexlsx:end_num_writexlsx),'wavelet_accuracy_single_var.xls','WriteMode','append')
end


%% nosense code
% Para_node = cell(repeat_num,scale_max*size(consider_node,2)*size(sigma_error,2));
% select_binary_order = cell(repeat_num,scale_max*size(consider_node,2)*size(sigma_error,2));
% select_binary = cell(repeat_num,scale_max*size(consider_node,2)*size(sigma_error,2));
% node_right = cell(repeat_num,scale_max*size(consider_node,2)*size(sigma_error,2));
% count_right_node = cell(repeat_num,scale_max*size(consider_node,2)*size(sigma_error,2));
% accuracy_node = zeros(repeat_num,scale_max*size(consider_node,2)*size(sigma_error,2));
% for i = 1:scale_max*size(consider_node,2)*size(sigma_error,2)
%     for j = 1:repeat_num
%         Para_node{j,i} = SVM_model{j,i};
%         select_binary_order{j,i} = find(abs(Para_node{j,i}) > 1.2);
%         select_binary{j,i} = abs(Para_node{j,i}) > 1.2;
%         node_right{j,i} = [node_idx_dim_ad_after{j,1}() node_idx_dim_cn_after{j,1}()];
%         count_right_node{j,i} = ismember(select_binary_order{j,i},node_right{j,i});
%         accuracy_node(j,i) = sum(count_right_node{j,i},2);
%     end
% end
% 
% save('results\Para_node.mat','Para_node');
% save('results\select_binary_order.mat','select_binary_order');
% save('results\select_binary.mat','select_binary');
