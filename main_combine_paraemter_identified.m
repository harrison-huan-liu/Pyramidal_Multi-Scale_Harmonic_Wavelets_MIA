close all;clear;clc
%% load brain connectivity network
load MST
p = 55;% p is the number of eignvectors %55
Graph = Preprocess_network_data('Data\DataTS.csv','Data\AD-Data\',p);
m = size(Graph,2);% m是研究对象的数目
n = size(Graph(1).L,1);% 每个研究对象的结点数目
% CommonNetwork = zeros(n);

%% calculate L_MST
% for i = 1:m
%     CommonNetwork=CommonNetwork+MST{i};
% end
% CommonNetwork=CommonNetwork/m;
% CommonNetwork(CommonNetwork<0.02)=0;
% CommonNetwork=(CommonNetwork+CommonNetwork')/2;
CommonNetwork = MST{m+2};
temp_D=diag(sum(CommonNetwork,2));
LatentLaplacian=temp_D-CommonNetwork;
[Phi_temp,~]=eig(LatentLaplacian);
Phi_ave=Phi_temp(:,1:p);

%% load Mask data
load u_vector
load u_vector_scale_1
load u_vector_scale_3
load region_mask
load region_mask_scale_1
load region_mask_scale_3
load node_idx
load node_idx_scale_1
load node_idx_scale_3

%% identify hyperparaemter
% paran = 1000;
nhub = size(u_vector_scale_1,1);
% u1 = cell(nhub,1);
% u2 = cell(nhub,1);
% for j = 1:paran
%     for i = 1:nhub
%         u1{i,j} = 0.1*(mod(j-1,10)+1);
%         u2{i,j} = 0.1*(floor((j-1)/10)+1);
%     end
% end

%% identified paraemter
paran = 10;
u1 = 0.8;
u2 = 0.7;

%% initialize phi_k1, phi_k2, phi_k3
Theta_k = cell(nhub,3);
phi_k = cell(nhub,4);

% for j = 1:paran
    for i = 1:nhub
%         Theta_k{i,1} = LatentLaplacian + u1{i,j}.*diag(ones(size(u_vector_scale_1{i}))-u_vector_scale_1{i}) + u2{i,j}.*(Phi_ave * Phi_ave.');
        Theta_k{i,1} = LatentLaplacian + u1.*diag(ones(size(u_vector_scale_1{i}))-u_vector_scale_1{i}) + u2.*(Phi_ave * Phi_ave.');
        [phi_k{i,1},~] = eig(Theta_k{i,1});
    %     phi_k{i,1} = diag(2.*ones(n,1));
%         Theta_k{i,2} = LatentLaplacian + u1{i,j}.*diag(ones(size(u_vector{i}))-u_vector{i}) + u2{i,j}.*(Phi_ave * Phi_ave.');
        Theta_k{i,2} = LatentLaplacian + u1.*diag(ones(size(u_vector{i}))-u_vector{i}) + u2.*(Phi_ave * Phi_ave.');
        [phi_k{i,2},~] = eig(Theta_k{i,2});
    %     phi_k{i,2} = phi_k{i,1};
%         Theta_k{i,3} = LatentLaplacian + u1{i,j}.*diag(ones(size(u_vector_scale_3{i}))-u_vector_scale_3{i}) + u2{i,j}.*(Phi_ave * Phi_ave.');
        Theta_k{i,3} = LatentLaplacian + u1.*diag(ones(size(u_vector_scale_3{i}))-u_vector_scale_3{i}) + u2.*(Phi_ave * Phi_ave.');
        [phi_k{i,3},~] = eig(Theta_k{i,3});
    %     phi_k{i,3} = phi_k{i,1};
    end
% end

%% select the number of wavelets
p1 = zeros(nhub,1);
p2 = zeros(nhub,1);
p3 = zeros(nhub,1);

hub_nodes = [3, 6, 10, 15, 16, 23, 24, 27, 30, 35, 44, 48, 51, 55, 56, 59, 65, 66, 68, 69, 73, 77, 80, 84, 89, 90, 93, 97, 98, 118, 122, 123, 127, 128, 129, 104, 140, 142, 144, 147];

for i = 1:nhub
%     for m = 1:size(node_idx_scale_1{i},1)
%         degree1 = degrees_und(region_mask_scale_1{i});
%         p1(i) = p1(i) + degree1(node_idx_scale_1{i}(m));
%     end
%     p1(i) = p1(i) / size(node_idx_scale_1{i},1);
%     for j = 1:size(node_idx{i},1)
%         degree2 = degrees_und(region_mask{i});
%         p2(i) = p2(i) + degree2(node_idx{i}(j));
%     end
%     p2(i) = p2(i) / size(node_idx{i},1);
%     for k = 1:size(node_idx_scale_3{i},1)
%         degree3 = degrees_und(region_mask_scale_3{i});
%         p3(i) = p3(i) + degree3(node_idx_scale_3{i}(k));
%     end
%     p3(i) = p3(i) / size(node_idx_scale_3{i},1);
    
    degree1 = degrees_und(region_mask_scale_1{i});
    p1(i) = degree1(hub_nodes(i));
    if p1(i) <= 2
        p1(i) = 3;
    end
    for j = 1:size(node_idx{i},1)
        degree2 = degrees_und(region_mask_scale_1{i});
        p2(i) = p2(i) + degree2(node_idx{i}(j));
    end
    p2(i) = 0.5 * (p2(i) - 2 * p1(i));
    for k = 1:size(node_idx_scale_3{i},1)
        degree3 = degrees_und(region_mask{i});
        p3(i) = p3(i) + degree3(node_idx_scale_3{i}(k));
    end
    p3(i) = 0.5 * (p3(i) - 2 * p2(i) - 2 * p1(i));
end

q1 = round(mean(p1));
q2 = round(mean(p2));
q3 = round(mean(p3));
for i = 1:nhub
    p1(i) = q1;
    p2(i) = q1;
    p3(i) = q1;
end
q = p1 + p2 + p3;

%% identify gama and alpha
gama1 = cell(nhub,1);
gama2 = cell(nhub,1);
gama3 = cell(nhub,1);
for i = 1:nhub
    gama1{i} = zeros(q(i),q(i));
    gama2{i} = zeros(q(i),q(i));
    gama3{i} = zeros(q(i),q(i));
    for j = 1:p1(i)
        gama1{i}(j,j) = 1;
    end
    for j = 1:p2(i)
        gama2{i}(p1(i)+j,p1(i)+j) = 1;
    end
    for j = 1:p3(i)
        gama3{i}(p1(i)+p2(i)+j,p1(i)+p2(i)+j) = 1;
    end
end
alpha = 0;
% for j = 1:paran
    for i = 1:nhub
%         [~,eigenvalue]=eig(LatentLaplacian+u2{1,j}.*(Phi_ave*Phi_ave.')+u1{1,j}.*diag(ones(size(u_vector_scale_1{i}))-u_vector_scale_1{i})+u1{1,j}.*diag(ones(size(u_vector{i}))-u_vector{i})+u1{1,j}.*diag(ones(size(u_vector_scale_3{i}))-u_vector_scale_3{i}));
        [~,eigenvalue]=eig(LatentLaplacian+u2.*(Phi_ave*Phi_ave.')+u1.*diag(ones(size(u_vector_scale_1{i}))-u_vector_scale_1{i})+u1.*diag(ones(size(u_vector{i}))-u_vector{i})+u1.*diag(ones(size(u_vector_scale_3{i}))-u_vector_scale_3{i}));
        if max(max(eigenvalue))>alpha
            alpha = max(max(eigenvalue));
        end
    end
% end

% while tem1>=0.01||tem2>=0.01||tem3>=0.01
%% solve phi_k
e = cell(2,2);
f = cell(2,2);
item1 = cell(2,2);
item2 = cell(2,2);
item3 = cell(2,2);
item4 = cell(2,2);
item5 = cell(2,2);
iterationnum = zeros(nhub,1);
    
    %% consider the first and fifth item
%     for i = 1:nhub
%         sim1 = alpha.* ones(n) - LatentLaplacian; % 
%         sim5 = Phi_ave * Phi_ave.';
%         phi_k{i,4} = [phi_k{i,1}(:,1:p1(i)) phi_k{i,2}(:,1:p2(i)) phi_k{i,3}(:,1:p3(i))];
%         looptime = 1;
%         e{looptime,i} = 1000;
%         f{looptime,i} = trace(phi_k{i,4}.'*(sim1 - u2{i,j}.*sim5)*phi_k{i,4});
%         item1{looptime,i} = trace(phi_k{i,4}.'*(sim1)*phi_k{i,4});
%         item5{looptime,i} = trace(phi_k{i,4}.'*(sim5)*phi_k{i,4});
% %         u1{i,j}=0.1*abs(item1{looptime,i})/(abs(item3{looptime,i}));
% %         u2{i,j}=100*abs(item1{looptime,i})/(abs(item5{looptime,i}));
%         while abs(e{looptime,i})>=0.01
% %             f = trace(phi_k{i,4}.'*(sim1)*phi_k{i,4});
%             M_gra = 2 .* (sim1 - u2{i,j}.*sim5) * phi_k{i,4};
%             [U,S,V] = svd(M_gra);
%             if q(i,1)>n
%                 V(:,149:q(i,1)) = [];
%                 phi_k{i,4} = U(:,1:n)*V.';
%             else
%                 phi_k{i,4} = U(:,1:q(i,1))*V.';
%             end
%             looptime = looptime + 1;
%             f{looptime,i} = trace(phi_k{i,4}.'*(sim1 - u2{i,j}.*sim5)*phi_k{i,4});
%             item1{looptime,i} = trace(phi_k{i,4}.'*sim1*phi_k{i,4});
%             item5{looptime,i} = trace(phi_k{i,4}.'*sim5*phi_k{i,4});
%             e{looptime,i} = f{looptime,i} - f{looptime-1,i};
%         end
%         iterationnum(i,1) = looptime;
%     end
    
     %% minimize
%     for i = 1:nhub
%         sim1 = LatentLaplacian; % 
%         sim2 = diag(ones(size(u_vector_scale_1{i}))-u_vector_scale_1{i});
%         sim3 = diag(ones(size(u_vector{i}))-u_vector{i});
%         sim4 = diag(ones(size(u_vector_scale_3{i}))-u_vector_scale_3{i});
%         sim5 = Phi_ave * Phi_ave.';
%         phi_k{i,4} = [phi_k{i,1}(:,1:p1(i)) phi_k{i,2}(:,1:p2(i)) phi_k{i,3}(:,1:p3(i))];
%         looptime = 1;
%         e{looptime,i} = 1000;
%         f{looptime,i} = trace(phi_k{i,4}.'*(sim1+u2{i,j}.*sim5)*phi_k{i,4} + u1{i,j}.*gama1{i}*phi_k{i,4}.'*sim2*phi_k{i,4} + u1{i,j}.*gama2{i}*phi_k{i,4}.'*sim3*phi_k{i,4} + u1{i,j}.*gama3{i}*phi_k{i,4}.'*sim4*phi_k{i,4});
%         item1{looptime,i} = trace(phi_k{i,4}.'*(sim1)*phi_k{i,4});
%         item2{looptime,i} = trace(gama1{i}*phi_k{i,4}.'*sim2*phi_k{i,4});
%         item3{looptime,i} = trace(gama2{i}*phi_k{i,4}.'*sim3*phi_k{i,4});
%         item4{looptime,i} = trace(gama3{i}*phi_k{i,4}.'*sim4*phi_k{i,4});
%         item5{looptime,i} = trace(phi_k{i,4}.'*(sim5)*phi_k{i,4});
% %         u1{i,j}=0.1*abs(item1{looptime,i})/(abs(item3{looptime,i}));
% %         u2{i,j}=100*abs(item1{looptime,i})/(abs(item5{looptime,i}));
%         while abs(e{looptime,i})>=0.01
% %             f = trace(phi_k{i,4}.'*(sim1)*phi_k{i,4}+gama1{i}*phi_k{i,4}.'*u1{i,j}.*sim2*phi_k{i,4}+gama2{i}*phi_k{i,4}.'*u1{i,j}.*sim3*phi_k{i,4}+gama3{i}*phi_k{i,4}.'*u1{i,j}.*sim4*phi_k{i,4});
%             M_gra = 2 .* (sim1 + u2{i,j}.*sim5) * phi_k{i,4} + 2.*(u1{i,j}.*sim2*phi_k{i,4}*gama1{i}+u1{i,j}.*sim3*phi_k{i,4}*gama2{i}+u1{i,j}.*sim4*phi_k{i,4}*gama3{i});
%             [U,S,V] = svd(M_gra);
%             if q(i,1)>n
%                 V(:,149:q(i,1)) = [];
%                 phi_k{i,4} = U(:,1:n)*V.';
%             else
%                 phi_k{i,4} = U(:,1:q(i,1))*V.';
%             end
%             looptime = looptime + 1;
%             f{looptime,i} = trace(phi_k{i,4}.'*(sim1+u2{i,j}.*sim5)*phi_k{i,4} + u1{i,j}.*gama1{i}*phi_k{i,4}.'*sim2*phi_k{i,4} + u1{i,j}.*gama2{i}*phi_k{i,4}.'*sim3*phi_k{i,4} + u1{i,j}.*gama3{i}*phi_k{i,4}.'*sim4*phi_k{i,4});
%             item1{looptime,i} = trace(phi_k{i,4}.'*sim1*phi_k{i,4});
%             item2{looptime,i} = trace(gama1{i}*phi_k{i,4}.'*sim2*phi_k{i,4});
%             item3{looptime,i} = trace(gama2{i}*phi_k{i,4}.'*sim3*phi_k{i,4});
%             item4{looptime,i} = trace(gama3{i}*phi_k{i,4}.'*sim4*phi_k{i,4});
%             item5{looptime,i} = trace(phi_k{i,4}.'*sim5*phi_k{i,4});
%             e{looptime,i} = f{looptime,i} - f{looptime-1,i};
%         end
%         iterationnum(i,1) = looptime;
%     end
    
    
%% maximize
inti = cell(paran,nhub);
% tem1 = cell(paran,nhub);
% tem2 = cell(paran,nhub);
% tem3 = cell(paran,nhub);
% tem4 = cell(paran,nhub);
% tem5 = cell(paran,nhub);
% tem6 = cell(paran,nhub);
% tem7 = cell(paran,nhub);
% tem8 = cell(paran,nhub);
% tem9 = cell(paran,nhub);
% tem10 = cell(paran,nhub);
tem11 = cell(paran,1);
tem12 = cell(paran,1);
tem13 = cell(paran,1);
tem14 = cell(paran,1);
    
[AD_amy,CN_amy,EMCI_amy,LMCI_amy,SMC_amy]= Proprocess_original_Amyloid_data('Data\signal_Data\Amyloid_SUVR.xlsx');
% [AD_tau,CN_tau,EMCI_tau,LMCI_tau,SMC_tau]= Proprocess_original_Tau_data('Data\signal_Data\Tau_SUVR.xlsx');
% [AD_fdg,CN_fdg,EMCI_fdg,LMCI_fdg,SMC_fdg]= Proprocess_original_FDG_data('Data\signal_Data\FDG_SUVR.xlsx');
AC_num = size(AD_amy.Amyloid,2)+size(CN_amy.Amyloid,2);
LC_num = size(LMCI_amy.Amyloid,2)+size(CN_amy.Amyloid,2);
CE_num = size(CN_amy.Amyloid,2)+size(EMCI_amy.Amyloid,2);
EL_num = size(EMCI_amy.Amyloid,2)+size(LMCI_amy.Amyloid,2);
    
accuracy_ac_fold = cell(paran,4);
accuracy_ce_fold = cell(paran,4);
accuracy_el_fold = cell(paran,4);
accuracy_lc_fold = cell(paran,4);

specificity_ac_fold = cell(paran,4);
specificity_ce_fold = cell(paran,4);
specificity_el_fold = cell(paran,4);
specificity_lc_fold = cell(paran,4);

sensitivity_ac_fold = cell(paran,4);
sensitivity_ce_fold = cell(paran,4);
sensitivity_el_fold = cell(paran,4);
sensitivity_lc_fold = cell(paran,4);

fscore_ac_fold = cell(paran,4);
fscore_ce_fold = cell(paran,4);
fscore_el_fold = cell(paran,4);
fscore_lc_fold = cell(paran,4);

accuracy_ac = zeros(paran,4);
accuracy_ce = zeros(paran,4);
accuracy_el = zeros(paran,4);
accuracy_lc = zeros(paran,4);
accuracy_ac_std = zeros(paran,4);
accuracy_ce_std = zeros(paran,4);
accuracy_el_std = zeros(paran,4);
accuracy_lc_std = zeros(paran,4);

specificity_ac = zeros(paran,4);
specificity_ce = zeros(paran,4);
specificity_el = zeros(paran,4);
specificity_lc = zeros(paran,4);
specificity_ac_std = zeros(paran,4);
specificity_ce_std = zeros(paran,4);
specificity_el_std = zeros(paran,4);
specificity_lc_std = zeros(paran,4);

sensitivity_ac = zeros(paran,4);
sensitivity_ce = zeros(paran,4);
sensitivity_el = zeros(paran,4);
sensitivity_lc = zeros(paran,4);
sensitivity_ac_std = zeros(paran,4);
sensitivity_ce_std = zeros(paran,4);
sensitivity_el_std = zeros(paran,4);
sensitivity_lc_std = zeros(paran,4);

fscore_ac = zeros(paran,4);
fscore_ce = zeros(paran,4);
fscore_el = zeros(paran,4);
fscore_lc = zeros(paran,4);
fscore_ac_std = zeros(paran,4);
fscore_ce_std = zeros(paran,4);
fscore_el_std = zeros(paran,4);
fscore_lc_std = zeros(paran,4);

% multiPhi = cell(paran,1);
% singlePhi = cell(paran,1);

for i = 1:nhub
    sim1 = alpha.* ones(n) - LatentLaplacian; % 
    sim2 = diag(ones(size(u_vector_scale_1{i}))-u_vector_scale_1{i});
    sim3 = diag(ones(size(u_vector{i}))-u_vector{i});
    sim4 = diag(ones(size(u_vector_scale_3{i}))-u_vector_scale_3{i});
    sim5 = Phi_ave * Phi_ave.';
    phi_k{i,4} = [phi_k{i,1}(:,1:p1(i)) phi_k{i,2}(:,1:p2(i)) phi_k{i,3}(:,1:p3(i))];
    inti{i,1} = phi_k{i,4}.'*phi_k{i,4};
    looptime = 1;
    e{looptime,i} = 1000;
%     f{looptime,i} = trace(phi_k{i,4}.'*(sim1-u2{i,j}.*sim5)*phi_k{i,4} - u1{i,j}.*gama1{i}*phi_k{i,4}.'*sim2*phi_k{i,4} - u1{i,j}.*gama2{i}*phi_k{i,4}.'*sim3*phi_k{i,4} - u1{i,j}.*gama3{i}*phi_k{i,4}.'*sim4*phi_k{i,4});
    f{looptime,i} = trace(phi_k{i,4}.'*(sim1-u2.*sim5)*phi_k{i,4} - u1.*gama1{i}*phi_k{i,4}.'*sim2*phi_k{i,4} - u1.*gama2{i}*phi_k{i,4}.'*sim3*phi_k{i,4} - u1.*gama3{i}*phi_k{i,4}.'*sim4*phi_k{i,4});
    item1{looptime,i} = trace(phi_k{i,4}.'*(sim1)*phi_k{i,4});
    item2{looptime,i} = trace(gama1{i}*phi_k{i,4}.'*sim2*phi_k{i,4});
    item3{looptime,i} = trace(gama2{i}*phi_k{i,4}.'*sim3*phi_k{i,4});
    item4{looptime,i} = trace(gama3{i}*phi_k{i,4}.'*sim4*phi_k{i,4});
    item5{looptime,i} = trace(phi_k{i,4}.'*(sim5)*phi_k{i,4});
%     u1{i,j}=0.1*abs(item1{looptime,i})/(abs(item3{looptime,i}));
%     u2{i,j}=100*abs(item1{looptime,i})/(abs(item5{looptime,i}));
    while abs(e{looptime,i})>=1
%         f = trace(phi_k{i,4}.'*(sim1)*phi_k{i,4}+gama1{i}*phi_k{i,4}.'*u1{i,j}.*sim2*phi_k{i,4}+gama2{i}*phi_k{i,4}.'*u1{i,j}.*sim3*phi_k{i,4}+gama3{i}*phi_k{i,4}.'*u1{i,j}.*sim4*phi_k{i,4});
%         M_gra = (sim1 - u2{i,j}.*sim5) * phi_k{i,4} - u1{i,j}.*sim2*phi_k{i,4}*gama1{i}-u1{i,j}.*sim3*phi_k{i,4}*gama2{i}-u1{i,j}.*sim4*phi_k{i,4}*gama3{i};
        M_gra = (sim1 - u2.*sim5) * phi_k{i,4} - u1.*sim2*phi_k{i,4}*gama1{i}-u1.*sim3*phi_k{i,4}*gama2{i}-u1.*sim4*phi_k{i,4}*gama3{i};
        %% full SVD
%         [U,S,V] = svd(M_gra);
%         if q(i,1)>n
%             V(:,149:q(i,1)) = [];
%             phi_k{i,4} = U(:,1:n)*V.';
%         else
%             phi_k{i,4} = U(:,1:q(i,1))*V.';
%         end
        %% compact SVD
        [U,S,V] = svd(M_gra,'econ');
        phi_k{i,4} = U*V.';
        looptime = looptime + 1;
%         f{looptime,i} = trace(phi_k{i,4}.'*(sim1-u2{i,j}.*sim5)*phi_k{i,4} - u1{i,j}.*gama1{i}*phi_k{i,4}.'*sim2*phi_k{i,4} - u1{i,j}.*gama2{i}*phi_k{i,4}.'*sim3*phi_k{i,4} - u1{i,j}.*gama3{i}*phi_k{i,4}.'*sim4*phi_k{i,4});
        f{looptime,i} = trace(phi_k{i,4}.'*(sim1-u2.*sim5)*phi_k{i,4} - u1.*gama1{i}*phi_k{i,4}.'*sim2*phi_k{i,4} - u1.*gama2{i}*phi_k{i,4}.'*sim3*phi_k{i,4} - u1.*gama3{i}*phi_k{i,4}.'*sim4*phi_k{i,4});
        item1{looptime,i} = trace(phi_k{i,4}.'*sim1*phi_k{i,4});
        item2{looptime,i} = trace(gama1{i}*phi_k{i,4}.'*sim2*phi_k{i,4});
        item3{looptime,i} = trace(gama2{i}*phi_k{i,4}.'*sim3*phi_k{i,4});
        item4{looptime,i} = trace(gama3{i}*phi_k{i,4}.'*sim4*phi_k{i,4});
        item5{looptime,i} = trace(phi_k{i,4}.'*sim5*phi_k{i,4});
        e{looptime,i} = f{looptime,i} - f{looptime-1,i};
    end
    iterationnum(i,1) = looptime;
end
    
multiPhi = phi_k{1,4};

for i = 2:nhub
    multiPhi = [multiPhi,phi_k{i,4}];
end


singlePhi = phi_k{1,4}(:,1:6);

for i = 2:nhub
    singlePhi = [singlePhi,phi_k{i,4}(:,1:6)];
end

% right_num_ac_single = zeros(paran,1);
% right_num_ce_single = zeros(paran,1);
% right_num_el_single = zeros(paran,1);
% right_num_lc_single = zeros(paran,1);
% right_num_ac_std_single = zeros(paran,1);
% right_num_ce_std_single = zeros(paran,1);
% right_num_el_std_single = zeros(paran,1);
% right_num_lc_std_single = zeros(paran,1);
% 
% right_num_ac_global = zeros(paran,1);
% right_num_ce_global = zeros(paran,1);
% right_num_el_global = zeros(paran,1);
% right_num_lc_global = zeros(paran,1);
% right_num_ac_std_global = zeros(paran,1);
% right_num_ce_std_global = zeros(paran,1);
% right_num_el_std_global = zeros(paran,1);
% right_num_lc_std_global = zeros(paran,1);
% 
% right_num_ac_original = zeros(paran,1);
% right_num_ce_original = zeros(paran,1);
% right_num_el_original = zeros(paran,1);
% right_num_lc_original = zeros(paran,1);
% right_num_ac_std_original = zeros(paran,1);
% right_num_ce_std_original = zeros(paran,1);
% right_num_el_std_original = zeros(paran,1);
% right_num_lc_std_original = zeros(paran,1);
    
for j = 1:paran
    [accuracy_ac_fold{j,1},sensitivity_ac_fold{j,1},specificity_ac_fold{j,1},fscore_ac_fold{j,1},tem11{j,1}] = SVM_signaldata(multiPhi.'*CN_amy.Amyloid,multiPhi.'*AD_amy.Amyloid,AC_num);
    accuracy_ac(j,1) = mean(accuracy_ac_fold{j,1});
    accuracy_ac_std(j,1) = std(accuracy_ac_fold{j,1});
    sensitivity_ac(j,1) = mean(sensitivity_ac_fold{j,1});
    sensitivity_ac_std(j,1) = std(sensitivity_ac_fold{j,1});
    specificity_ac(j,1) = mean(specificity_ac_fold{j,1});
    specificity_ac_std(j,1) = std(specificity_ac_fold{j,1});
    fscore_ac(j,1) = mean(fscore_ac_fold{j,1});
    fscore_ac_std(j,1) = std(fscore_ac_fold{j,1});
    [accuracy_ce_fold{j,1},sensitivity_ce_fold{j,1},specificity_ce_fold{j,1},fscore_ce_fold{j,1},tem12{j,1}] = SVM_signaldata(multiPhi.'*CN_amy.Amyloid,multiPhi.'*EMCI_amy.Amyloid,CE_num);
    accuracy_ce(j,1) = mean(accuracy_ce_fold{j,1});
    accuracy_ce_std(j,1) = std(accuracy_ce_fold{j,1});
    sensitivity_ce(j,1) = mean(sensitivity_ce_fold{j,1});
    sensitivity_ce_std(j,1) = std(sensitivity_ce_fold{j,1});
    specificity_ce(j,1) = mean(specificity_ce_fold{j,1});
    specificity_ce_std(j,1) = std(specificity_ce_fold{j,1});
    fscore_ce(j,1) = mean(fscore_ce_fold{j,1});
    fscore_ce_std(j,1) = std(fscore_ce_fold{j,1});
    [accuracy_el_fold{j,1},sensitivity_el_fold{j,1},specificity_el_fold{j,1},fscore_el_fold{j,1},tem13{j,1}] = SVM_signaldata(multiPhi.'*EMCI_amy.Amyloid,multiPhi.'*LMCI_amy.Amyloid,EL_num);
    accuracy_el(j,1) = mean(accuracy_el_fold{j,1});
    accuracy_el_std(j,1) = std(accuracy_el_fold{j,1});
    sensitivity_el(j,1) = mean(sensitivity_el_fold{j,1});
    sensitivity_el_std(j,1) = std(sensitivity_el_fold{j,1});
    specificity_el(j,1) = mean(specificity_el_fold{j,1});
    specificity_el_std(j,1) = std(specificity_el_fold{j,1});
    fscore_el(j,1) = mean(fscore_el_fold{j,1});
    fscore_el_std(j,1) = std(fscore_el_fold{j,1});
    [accuracy_lc_fold{j,1},sensitivity_lc_fold{j,1},specificity_lc_fold{j,1},fscore_lc_fold{j,1},tem14{j,1}] = SVM_signaldata(multiPhi.'*CN_amy.Amyloid,multiPhi.'*LMCI_amy.Amyloid,LC_num);
    accuracy_lc(j,1) = mean(accuracy_lc_fold{j,1});
    accuracy_lc_std(j,1) = std(accuracy_lc_fold{j,1});
    sensitivity_lc(j,1) = mean(sensitivity_lc_fold{j,1});
    sensitivity_lc_std(j,1) = std(sensitivity_lc_fold{j,1});
    specificity_lc(j,1) = mean(specificity_lc_fold{j,1});
    specificity_lc_std(j,1) = std(specificity_lc_fold{j,1});
    fscore_lc(j,1) = mean(fscore_lc_fold{j,1});
    fscore_lc_std(j,1) = std(fscore_lc_fold{j,1});
    
    [accuracy_ac_fold{j,2},sensitivity_ac_fold{j,2},specificity_ac_fold{j,2},fscore_ac_fold{j,2},~] = SVM_signaldata(singlePhi.'*CN_amy.Amyloid,singlePhi.'*AD_amy.Amyloid,AC_num);
    accuracy_ac(j,2) = mean(accuracy_ac_fold{j,2});
    accuracy_ac_std(j,2) = std(accuracy_ac_fold{j,2});
    sensitivity_ac(j,2) = mean(sensitivity_ac_fold{j,2});
    sensitivity_ac_std(j,2) = std(sensitivity_ac_fold{j,2});
    specificity_ac(j,2) = mean(specificity_ac_fold{j,2});
    specificity_ac_std(j,2) = std(specificity_ac_fold{j,2});
    fscore_ac(j,2) = mean(fscore_ac_fold{j,2});
    fscore_ac_std(j,2) = std(fscore_ac_fold{j,2});
    [accuracy_ce_fold{j,2},sensitivity_ce_fold{j,2},specificity_ce_fold{j,2},fscore_ce_fold{j,2},~] = SVM_signaldata(singlePhi.'*CN_amy.Amyloid,singlePhi.'*EMCI_amy.Amyloid,CE_num);
    accuracy_ce(j,2) = mean(accuracy_ce_fold{j,2});
    accuracy_ce_std(j,2) = std(accuracy_ce_fold{j,2});
    sensitivity_ce(j,2) = mean(sensitivity_ce_fold{j,2});
    sensitivity_ce_std(j,2) = std(sensitivity_ce_fold{j,2});
    specificity_ce(j,2) = mean(specificity_ce_fold{j,2});
    specificity_ce_std(j,2) = std(specificity_ce_fold{j,2});
    fscore_ce(j,2) = mean(fscore_ce_fold{j,2});
    fscore_ce_std(j,2) = std(fscore_ce_fold{j,2});
    [accuracy_el_fold{j,2},sensitivity_el_fold{j,2},specificity_el_fold{j,2},fscore_el_fold{j,2},~] = SVM_signaldata(singlePhi.'*EMCI_amy.Amyloid,singlePhi.'*LMCI_amy.Amyloid,EL_num);
    accuracy_el(j,2) = mean(accuracy_el_fold{j,2});
    accuracy_el_std(j,2) = std(accuracy_el_fold{j,2});
    sensitivity_el(j,2) = mean(sensitivity_el_fold{j,2});
    sensitivity_el_std(j,2) = std(sensitivity_el_fold{j,2});
    specificity_el(j,2) = mean(specificity_el_fold{j,2});
    specificity_el_std(j,2) = std(specificity_el_fold{j,2});
    fscore_el(j,2) = mean(fscore_el_fold{j,2});
    fscore_el_std(j,2) = std(fscore_el_fold{j,2});
    [accuracy_lc_fold{j,2},sensitivity_lc_fold{j,2},specificity_lc_fold{j,2},fscore_lc_fold{j,2},~] = SVM_signaldata(singlePhi.'*CN_amy.Amyloid,singlePhi.'*LMCI_amy.Amyloid,LC_num);
    accuracy_lc(j,2) = mean(accuracy_lc_fold{j,2});
    accuracy_lc_std(j,2) = std(accuracy_lc_fold{j,2});
    sensitivity_lc(j,2) = mean(sensitivity_lc_fold{j,2});
    sensitivity_lc_std(j,2) = std(sensitivity_lc_fold{j,2});
    specificity_lc(j,2) = mean(specificity_lc_fold{j,2});
    specificity_lc_std(j,2) = std(specificity_lc_fold{j,2});
    fscore_lc(j,2) = mean(fscore_lc_fold{j,2});
    fscore_lc_std(j,2) = std(fscore_lc_fold{j,2});
    
    [accuracy_ac_fold{j,3},sensitivity_ac_fold{j,3},specificity_ac_fold{j,3},fscore_ac_fold{j,3},~] = SVM_signaldata(Phi_ave.'*CN_amy.Amyloid,Phi_ave.'*AD_amy.Amyloid,AC_num);
    accuracy_ac(j,3) = mean(accuracy_ac_fold{j,3});
    accuracy_ac_std(j,3) = std(accuracy_ac_fold{j,3});
    sensitivity_ac(j,3) = mean(sensitivity_ac_fold{j,3});
    sensitivity_ac_std(j,3) = std(sensitivity_ac_fold{j,3});
    specificity_ac(j,3) = mean(specificity_ac_fold{j,3});
    specificity_ac_std(j,3) = std(specificity_ac_fold{j,3});
    fscore_ac(j,3) = mean(fscore_ac_fold{j,3});
    fscore_ac_std(j,3) = std(fscore_ac_fold{j,3});
    [accuracy_ce_fold{j,3},sensitivity_ce_fold{j,3},specificity_ce_fold{j,3},fscore_ce_fold{j,3},~] = SVM_signaldata(Phi_ave.'*CN_amy.Amyloid,Phi_ave.'*EMCI_amy.Amyloid,CE_num);
    accuracy_ce(j,3) = mean(accuracy_ce_fold{j,3});
    accuracy_ce_std(j,3) = std(accuracy_ce_fold{j,3});
    sensitivity_ce(j,3) = mean(sensitivity_ce_fold{j,3});
    sensitivity_ce_std(j,3) = std(sensitivity_ce_fold{j,3});
    specificity_ce(j,3) = mean(specificity_ce_fold{j,3});
    specificity_ce_std(j,3) = std(specificity_ce_fold{j,3});
    fscore_ce(j,3) = mean(fscore_ce_fold{j,3});
    fscore_ce_std(j,3) = std(fscore_ce_fold{j,3});
    [accuracy_el_fold{j,3},sensitivity_el_fold{j,3},specificity_el_fold{j,3},fscore_el_fold{j,3},~] = SVM_signaldata(Phi_ave.'*EMCI_amy.Amyloid,Phi_ave.'*LMCI_amy.Amyloid,EL_num);
    accuracy_el(j,3) = mean(accuracy_el_fold{j,3});
    accuracy_el_std(j,3) = std(accuracy_el_fold{j,3});
    sensitivity_el(j,3) = mean(sensitivity_el_fold{j,3});
    sensitivity_el_std(j,3) = std(sensitivity_el_fold{j,3});
    specificity_el(j,3) = mean(specificity_el_fold{j,3});
    specificity_el_std(j,3) = std(specificity_el_fold{j,3});
    fscore_el(j,3) = mean(fscore_el_fold{j,3});
    fscore_el_std(j,3) = std(fscore_el_fold{j,3});
    [accuracy_lc_fold{j,3},sensitivity_lc_fold{j,3},specificity_lc_fold{j,3},fscore_lc_fold{j,3},~] = SVM_signaldata(Phi_ave.'*CN_amy.Amyloid,Phi_ave.'*LMCI_amy.Amyloid,LC_num);
    accuracy_lc(j,3) = mean(accuracy_lc_fold{j,3});
    accuracy_lc_std(j,3) = std(accuracy_lc_fold{j,3});
    sensitivity_lc(j,3) = mean(sensitivity_lc_fold{j,3});
    sensitivity_lc_std(j,3) = std(sensitivity_lc_fold{j,3});
    specificity_lc(j,3) = mean(specificity_lc_fold{j,3});
    specificity_lc_std(j,3) = std(specificity_lc_fold{j,3});
    fscore_lc(j,3) = mean(fscore_lc_fold{j,3});
    fscore_lc_std(j,3) = std(fscore_lc_fold{j,3});
    
    [accuracy_ac_fold{j,4},sensitivity_ac_fold{j,4},specificity_ac_fold{j,4},fscore_ac_fold{j,4},~] = SVM_signaldata(CN_amy.Amyloid,AD_amy.Amyloid,AC_num);
    accuracy_ac(j,4) = mean(accuracy_ac_fold{j,4});
    accuracy_ac_std(j,4) = std(accuracy_ac_fold{j,4});
    sensitivity_ac(j,4) = mean(sensitivity_ac_fold{j,4});
    sensitivity_ac_std(j,4) = std(sensitivity_ac_fold{j,4});
    specificity_ac(j,4) = mean(specificity_ac_fold{j,4});
    specificity_ac_std(j,4) = std(specificity_ac_fold{j,4});
    fscore_ac(j,4) = mean(fscore_ac_fold{j,4});
    fscore_ac_std(j,4) = std(fscore_ac_fold{j,4});
    [accuracy_ce_fold{j,4},sensitivity_ce_fold{j,4},specificity_ce_fold{j,4},fscore_ce_fold{j,4},~] = SVM_signaldata(CN_amy.Amyloid,EMCI_amy.Amyloid,CE_num);
    accuracy_ce(j,4) = mean(accuracy_ce_fold{j,4});
    accuracy_ce_std(j,4) = std(accuracy_ce_fold{j,4});
    sensitivity_ce(j,4) = mean(sensitivity_ce_fold{j,4});
    sensitivity_ce_std(j,4) = std(sensitivity_ce_fold{j,4});
    specificity_ce(j,4) = mean(specificity_ce_fold{j,4});
    specificity_ce_std(j,4) = std(specificity_ce_fold{j,4});
    fscore_ce(j,4) = mean(fscore_ce_fold{j,4});
    fscore_ce_std(j,4) = std(fscore_ce_fold{j,4});
    [accuracy_el_fold{j,4},sensitivity_el_fold{j,4},specificity_el_fold{j,4},fscore_el_fold{j,4},~] = SVM_signaldata(EMCI_amy.Amyloid,LMCI_amy.Amyloid,EL_num);
    accuracy_el(j,4) = mean(accuracy_el_fold{j,4});
    accuracy_el_std(j,4) = std(accuracy_el_fold{j,4});
    sensitivity_el(j,4) = mean(sensitivity_el_fold{j,4});
    sensitivity_el_std(j,4) = std(sensitivity_el_fold{j,4});
    specificity_el(j,4) = mean(specificity_el_fold{j,4});
    specificity_el_std(j,4) = std(specificity_el_fold{j,4});
    fscore_el(j,4) = mean(fscore_el_fold{j,4});
    fscore_el_std(j,4) = std(fscore_el_fold{j,4});
    [accuracy_lc_fold{j,4},sensitivity_lc_fold{j,4},specificity_lc_fold{j,4},fscore_lc_fold{j,4},~] = SVM_signaldata(CN_amy.Amyloid,LMCI_amy.Amyloid,LC_num);
    accuracy_lc(j,4) = mean(accuracy_lc_fold{j,4});
    accuracy_lc_std(j,4) = std(accuracy_lc_fold{j,4});
    sensitivity_lc(j,4) = mean(sensitivity_lc_fold{j,4});
    sensitivity_lc_std(j,4) = std(sensitivity_lc_fold{j,4});
    specificity_lc(j,4) = mean(specificity_lc_fold{j,4});
    specificity_lc_std(j,4) = std(specificity_lc_fold{j,4});
    fscore_lc(j,4) = mean(fscore_lc_fold{j,4});
    fscore_lc_std(j,4) = std(fscore_lc_fold{j,4});
    
%     [right_num_ac_fold_single,~] = SVM_signaldata(singlePhi{j,1}.'*CN_amy.Amyloid,singlePhi{j,1}.'*AD_amy.Amyloid,AC_num);
%     right_num_ac_single(j,1) = mean(right_num_ac_fold_single)/floor(AC_num/10);
%     right_num_ac_std_single(j,1) = std(right_num_ac_fold_single)/floor(AC_num/10);
%     [right_num_ce_fold_single,~] = SVM_signaldata(singlePhi{j,1}.'*CN_amy.Amyloid,singlePhi{j,1}.'*EMCI_amy.Amyloid,CE_num);
%     right_num_ce_single(j,1) = mean(right_num_ce_fold_single)/floor(CE_num/10);
%     right_num_ce_std_single(j,1) = std(right_num_ce_fold_single)/floor(CE_num/10);
%     [right_num_el_fold_single,~] = SVM_signaldata(singlePhi{j,1}.'*EMCI_amy.Amyloid,singlePhi{j,1}.'*LMCI_amy.Amyloid,EL_num);
%     right_num_el_single(j,1) = mean(right_num_el_fold_single)/floor(EL_num/10);
%     right_num_el_std_single(j,1) = std(right_num_el_fold_single)/floor(EL_num/10);
%     [right_num_lc_fold_single,~] = SVM_signaldata(singlePhi{j,1}.'*CN_amy.Amyloid,singlePhi{j,1}.'*LMCI_amy.Amyloid,LC_num);
%     right_num_lc_single(j,1) = mean(right_num_lc_fold_single)/floor(LC_num/10);
%     right_num_lc_std_single(j,1) = std(right_num_lc_fold_single)/floor(LC_num/10);
%     
%     [right_num_ac_fold_global,~] = SVM_signaldata(Phi_ave.'*CN_amy.Amyloid,Phi_ave.'*AD_amy.Amyloid,AC_num);
%     right_num_ac_global(j,1) = mean(right_num_ac_fold_global)/floor(AC_num/10);
%     right_num_ac_std_global(j,1) = std(right_num_ac_fold_global)/floor(AC_num/10);
%     [right_num_ce_fold_global,~] = SVM_signaldata(Phi_ave.'*CN_amy.Amyloid,Phi_ave.'*EMCI_amy.Amyloid,CE_num);
%     right_num_ce_global(j,1) = mean(right_num_ce_fold_global)/floor(CE_num/10);
%     right_num_ce_std_global(j,1) = std(right_num_ce_fold_global)/floor(CE_num/10);
%     [right_num_el_fold_global,~] = SVM_signaldata(Phi_ave.'*EMCI_amy.Amyloid,Phi_ave.'*LMCI_amy.Amyloid,EL_num);
%     right_num_el_global(j,1) = mean(right_num_el_fold_global)/floor(EL_num/10);
%     right_num_el_std_global(j,1) = std(right_num_el_fold_global)/floor(EL_num/10);
%     [right_num_lc_fold_global,~] = SVM_signaldata(Phi_ave.'*CN_amy.Amyloid,Phi_ave.'*LMCI_amy.Amyloid,LC_num);
%     right_num_lc_global(j,1) = mean(right_num_lc_fold_global)/floor(LC_num/10);
%     right_num_lc_std_global(j,1) = std(right_num_lc_fold_global)/floor(LC_num/10);
% 
%     [right_num_ac_fold_original,~] = SVM_signaldata(CN_amy.Amyloid,AD_amy.Amyloid,AC_num);
%     right_num_ac_original(j,1) = mean(right_num_ac_fold_original)/floor(AC_num/10);
%     right_num_ac_std_original(j,1) = std(right_num_ac_fold_original)/floor(AC_num/10);
%     [right_num_ce_fold_original,~] = SVM_signaldata(CN_amy.Amyloid,EMCI_amy.Amyloid,CE_num);
%     right_num_ce_original(j,1) = mean(right_num_ce_fold_original)/floor(CE_num/10);
%     right_num_ce_std_original(j,1) = std(right_num_ce_fold_original)/floor(CE_num/10);
%     [right_num_el_fold_original,~] = SVM_signaldata(EMCI_amy.Amyloid,LMCI_amy.Amyloid,EL_num);
%     right_num_el_original(j,1) = mean(right_num_el_fold_original)/floor(EL_num/10);
%     right_num_el_std_original(j,1) = std(right_num_el_fold_original)/floor(EL_num/10);
%     [right_num_lc_fold_original,~] = SVM_signaldata(CN_amy.Amyloid,LMCI_amy.Amyloid,LC_num);
%     right_num_lc_original(j,1) = mean(right_num_lc_fold_original)/floor(LC_num/10);
%     right_num_lc_std_original(j,1) = std(right_num_lc_fold_original)/floor(LC_num/10);
end
    
%     figure(1)
% %     iterationnum = 20*ones(nhub,1);
%     for i = 1:nhub
%         x = zeros(iterationnum(i,1)-1,1);
%         y = zeros(iterationnum(i,1)-1,1);
%         for j = 1:iterationnum(i,1)-1
%             x(j,1) = j - 1;
%             y(j,1) = f{j+1,i};
%         end
%         subplot(7,7,i);
%         plot(x,y);
%     end
    
%     r = randi(nhub);
%     
%     tem10{j} = zeros(n,size(AD_amy.Amyloid,2));
%     tem11{j} = zeros(n,size(CN_amy.Amyloid,2));
%     tem12{j} = zeros(n,size(EMCI_amy.Amyloid,2));
%     tem13{j} = zeros(n,size(LMCI_amy.Amyloid,2));
%     tem14{j} = zeros(n,size(SMC_amy.Amyloid,2));

%     for r = 1:nhub
%         ex1 = abs(phi_k{r,4}(:,1:p1(r)).'*phi_k{r,4}(:,1:p1(r))-eye(p1(r)))./ones(p1(r));
%         [tem1{j,r},] = max(ex1(:));
%         ex2 = abs(phi_k{r,4}(:,p1(r)+1:p1(r)+p2(r)).'*phi_k{r,4}(:,p1(r)+1:p1(r)+p2(r))-eye(p2(r)))./ones(p2(r));
%         [tem2{j,r},] = max(ex2(:));
%         ex3 = abs(phi_k{r,4}(:,p1(r)+p2(r)+1:p1(r)+p2(r)+p3(r)).'*phi_k{r,4}(:,p1(r)+p2(r)+1:p1(r)+p2(r)+p3(r))-eye(p3(r)))./ones(p3(r));
%         [tem3{j,r},] = max(ex3(:));
%         ex4 = abs(phi_k{r,4}(:,1:p1(r)).'*phi_k{r,4}(:,p1(r)+1:p1(r)+p2(r)));
%         [tem4{j,r},] = max(ex4(:));
%         fprintf('The hub %d in the %d loop: \n',r,j);
%         fprintf('The error of scale 1 and 2 is %f \n',tem4{j,r});
%         ex5 = abs(phi_k{r,4}(:,1:p1(r)).'*phi_k{r,4}(:,p1(r)+p2(r)+1:p1(r)+p2(r)+p3(r)));
%         [tem5{j,r},] = max(ex5(:));
%         fprintf('The error of scale 1 and 3 is %f \n',tem5{j,r});
%         ex6 = abs(phi_k{r,4}(:,p1(r)+1:p1(r)+p2(r)).'*phi_k{r,4}(:,p1(r)+p2(r)+1:p1(r)+p2(r)+p3(r)));
%         [tem6{j,r},] = max(ex6(:));
%         fprintf('The error of scale 2 and 3 is %f \n',tem6{j,r});
%         ex7 = abs(phi_k{r,4}(:,1:p1(r)).'*Phi_ave);
%         [tem7{j,r},] = max(ex7(:));
%         fprintf('The error of scale 1 and 4 is %f \n',tem7{j,r});
%         ex8 = abs(phi_k{r,4}(:,p1(r)+1:p1(r)+p2(r)).'*Phi_ave);
%         [tem8{j,r},] = max(ex8(:));
%         fprintf('The error of scale 2 and 4 is %f \n',tem8{j,r});
%         ex9 = abs(phi_k{r,4}(:,p1(r)+p2(r)+1:p1(r)+p2(r)+p3(r)).'*Phi_ave);
%         [tem9{j,r},] = max(ex9(:));
%         fprintf('The error of scale 3 and 4 is %f \n',tem9{j,r});
%         fprintf('\n');
%     end

%         for sample = 1:size(AD_amy.Amyloid,2)
%             tem10{j,r}(:,sample) = Phi_ave * Phi_ave.' * (AD_amy.Amyloid(:,sample).*u_vector_scale_1{i}) + phi_k{r,4}(:,1:p1(r))*phi_k{r,4}(:,1:p1(r)).'*AD_amy.Amyloid(:,sample);
% %             tem10{j,r}(:,sample) = Phi_ave * Phi_ave.' * (AD_amy.Amyloid(:,sample).*u_vector{i}) + phi_k{r,4}(:,p1(r)+1:p1(r)+p2(r))*phi_k{r,4}(:,p1(r)+1:p1(r)+p2(r)).'*AD_amy.Amyloid(:,sample);
% %             tem10{j,r}(:,sample) = Phi_ave * Phi_ave.' * (AD_amy.Amyloid(:,sample).*u_vector_scale_3{i}) + phi_k{r,4}(:,p1(r)+p2(r)+1:p1(r)+p2(r)+p3(r))*phi_k{r,4}(:,p1(r)+p2(r)+1:p1(r)+p2(r)+p3(r)).'*AD_amy.Amyloid(:,sample);
%         end
%         for sample = 1:size(CN_amy.Amyloid,2)
%             tem11{j,r}(:,sample) = phi_k{r,4}*phi_k{r,4}.'*CN_amy.Amyloid(:,sample);
%         end
%         for sample = 1:size(EMCI_amy.Amyloid,2)
%             tem12{j,r}(:,sample) = phi_k{r,4}*phi_k{r,4}.'*EMCI_amy.Amyloid(:,sample);
%         end
%         for sample = 1:size(LMCI_amy.Amyloid,2)
%             tem13{j,r}(:,sample) = phi_k{r,4}*phi_k{r,4}.'*LMCI_amy.Amyloid(:,sample);
%         end
%         for sample = 1:size(SMC_amy.Amyloid,2)
%             tem14{j,r}(:,sample) = phi_k{r,4}*phi_k{r,4}.'*SMC_amy.Amyloid(:,sample);
%         end




accuracy_ac = right_num_ac/floor(AC_num/10);
% identify_para = mean(right_num_ac,2)/floor(AC_num/10);
[para_value,para_index] = max(accuracy_ac);

accuracy_ce = right_num_ce/floor(CE_num/10);
% identify_para_lc = mean(right_num_lc,2)/floor(LC_num/10);
% [para_value_lc,para_index_lc] = max(accuracy_lc);

accuracy_el = right_num_el/floor(EL_num/10);
% identify_para_lc = mean(right_num_lc,2)/floor(LC_num/10);
% [para_value_lc,para_index_lc] = max(accuracy_lc);

accuracy_lc = right_num_lc/floor(LC_num/10);
% identify_para_lc = mean(right_num_lc,2)/floor(LC_num/10);
[para_value_lc,para_index_lc] = max(accuracy_lc);

mean(accuracy_ce)
std(accuracy_ce)
mean(accuracy_el)
std(accuracy_el)
mean(accuracy_lc)
std(accuracy_lc)

para_index = 10;
SVM_para = cell(1,nhub);
for i = 1:nhub
    SVM_para{1,i} = tem14{para_index,1}(1+18*(i-1):18*i,1:10);
end

% end
% save('tem10.mat','tem10');
save('results\tem11.mat','tem11');
save('results\tem12.mat','tem12');
save('results\tem13.mat','tem13');
save('results\tem14.mat','tem14');
save('results\SVM_para.mat','SVM_para');
save('results\phi_k.mat','phi_k');
save('results\right_num_ac_fold.mat','right_num_ac_fold');
save('results\right_num_ce_fold.mat','right_num_ce_fold');
save('results\right_num_el_fold.mat','right_num_el_fold');
save('results\right_num_lc_fold.mat','right_num_lc_fold');

writematrix(right_num_ce_fold{para_index,1}.'/floor(LC_num/10),'accuracy_ce.xls','WriteMode','append')
writematrix(right_num_el_fold{para_index,1}.'/floor(LC_num/10),'accuracy_el.xls','WriteMode','append')
writematrix(right_num_lc_fold{para_index,1}.'/floor(LC_num/10),'accuracy_lc.xls','WriteMode','append')

for writexlsx = 1:nhub
    writematrix(SVM_para{1,writexlsx}.','SVM_weight_lc.xls','WriteMode','append')
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
    select_scale{1,i} = mean(abs(SVM_para{1,i}),2);
    select_scale1{1,i} = select_scale{1,i}(1:6,1);
    select_scale2{1,i} = select_scale{1,i}(7:12,1);
    select_scale3{1,i} = select_scale{1,i}(13:18,1);
    [select_m_scale1(1,i),select_wavelets1(1,i)] = max(select_scale1{1,i});
    [select_m_scale2(1,i),select_wavelets2(1,i)] = max(select_scale2{1,i});
    [select_m_scale3(1,i),select_wavelets3(1,i)] = max(select_scale3{1,i});
end
for writexlsx = 1:nhub
    writematrix(select_scale{1,writexlsx}.','select_scale.xls','WriteMode','append')
end
[~,node_index1] = maxk(select_m_scale1,5);
[~,node_index2] = maxk(select_m_scale2,5);
[~,node_index3] = maxk(select_m_scale3,6);
select_m_scale = [select_m_scale1,select_m_scale2,select_m_scale3];
[~,node_index123] = maxk(select_m_scale,15);

CommonNetwork = zeros(n);
for i = 1:m
    CommonNetwork=CommonNetwork+Graph(i).W;
end
CommonNetwork=CommonNetwork/m;
CommonNetwork(CommonNetwork<0.02)=0;
CommonNetwork=(CommonNetwork+CommonNetwork')/2;

d5 = zeros(1,15);
d6 = zeros(1,15);
d7 = zeros(1,15);
d8 = zeros(1,15);
d1 = betweenness_wei(CommonNetwork);
d2 = degrees_und(CommonNetwork);
d3 = efficiency_wei(CommonNetwork,2);
d4 = clustering_coef_wd(CommonNetwork);
for i = 1:15
    d5(1,i) = d1(node_index_index(i),1);
    d6(1,i) = d2(1,node_index_index(i));
    d7(1,i) = d3(node_index_index(i),1);
    d8(1,i) = d4(node_index_index(i),1);
end
writematrix(d5,'betweenness_wei0225.xls','WriteMode','append')
writematrix(d6,'degrees_und0225.xls','WriteMode','append')
writematrix(d7,'efficiency_wei0225.xls','WriteMode','append')
writematrix(d8,'clustering_coef_wd0225.xls','WriteMode','append')

node_index3([4]) = [];
node_index = [node_index1,node_index2,node_index3];
node_index_index = hub_nodes(node_index);
Eglob = zeros(1,15);
Eloc_select = zeros(1,15);
Eloc = efficiency_wei(CommonNetwork,2);
Networkk = zeros(148,148);
for i = 1:15
    for j = 1:size(node_idx_scale_3{node_index(i),1},1)
        Networkk(node_idx_scale_3{node_index(i),1}(j,1),:) = CommonNetwork(node_idx_scale_3{node_index(i),1}(j,1),:);
        Networkk(:,node_idx_scale_3{node_index(i),1}(j,1)) = CommonNetwork(:,node_idx_scale_3{node_index(i),1}(j,1));
    end
    Eglob(1,i) = efficiency_wei(Networkk);
    Eloc_select(1,i) = Eloc(node_index_index(i),1);
end

writematrix(Eglob,'Eglob.xls','WriteMode','append')
writematrix(Eloc_select,'Eloc_select.xls','WriteMode','append')
writematrix(node_index,'node_index0225.xls','WriteMode','append')
writematrix(node_index_index,'node_index_index0225.xls','WriteMode','append')

brain_node = zeros(148,1);
for i = 1:15
    brain_node(node_index_index(i)) = 1;
end
save('results\brain_significant_node.txt','brain_node','-ascii');

CEL_num = size(CN_amy.Amyloid,2)+size(EMCI_amy.Amyloid,2)+size(LMCI_amy.Amyloid,2);
Energy_signdata = zeros(CEL_num,1);
Energy_signdata_square = zeros(CEL_num,1);
Energy_signdata_each = zeros(CEL_num,720);
LC_signdata = [CN_amy.Amyloid,EMCI_amy.Amyloid,LMCI_amy.Amyloid];
for i = 1:CEL_num
    for j = 1:720
        Energy_signdata(i,1) = Energy_signdata(i,1) + dot(LC_signdata(1:148,i),multiPhi{para_index,1}(1:148,j));
        Energy_signdata_square(i,1) = Energy_signdata(i,1)*Energy_signdata(i,1);
        Energy_signdata_each(i,j) = dot(LC_signdata(1:148,i),multiPhi{para_index,1}(1:148,j));
    end
end

colmin = min(Energy_signdata_each);
colmax = max(Energy_signdata_each);
Energy_signdata_norm = rescale(Energy_signdata_each,'InputMin',colmin,'InputMax',colmax);

Energy_signdata_norm_num = sum(Energy_signdata_norm,2);
Energy_signdata_norm_num_square = sum(Energy_signdata_norm.*Energy_signdata_norm,2);

writematrix(Energy_signdata_norm_num_square,'energy.xls','WriteMode','append')

draw_wavelet1 = zeros(148,1);
draw_wavelet2 = zeros(148,1);
draw_wavelet3 = zeros(148,1);
for i = 1:148
    if abs(phi_k{9,4}(i,5)) > 0.1
        draw_wavelet1(i,1) = phi_k{9,4}(i,5);
    end
    if abs(phi_k{39,4}(i,7)) > 0.02
        draw_wavelet2(i,1) = phi_k{39,4}(i,7);
    end
    if abs(phi_k{6,4}(i,13)) > 0.01
        draw_wavelet3(i,1) = phi_k{6,4}(i,13);
    end
end

b720 = mean(abs(tem14{para_index,1}),2);
[b721, b722] = maxk(abs(b720),25);
b723 = mod(b722,18);
b724 = [1,2,3,4,5,8,15,17,25];
b725 = zeros(148,9);
for i = 1:9
    b725(:,i)=multiPhi{para_index,1}(:,b722(b724(i)));
end
for i = 1:148
    for j = 1:9
        if abs(b725(i,j))<0.1
            b725(i,j) = 0;
        end
    end
end

node_order = [5, 1, 2, 3, 4, 6, 9, 10, 7, 8, 11, 15, 13, 12, 14];
node_index_index = node_index_index(node_order);


% for i = 1:paran
%     tem11{i,1}=mean(tem11{i,1},'all');
% end

%% plot wavelets
% figure(2)
% for i = 1:4
%     subplot(2,2,i);
%     imagesc(phi_k{i+20,4}(:,1:p1(i+20)))
% end
% figure(3)
% for i = 1:4
%     subplot(2,2,i);
%     imagesc(phi_k{i+20,4}(:,p1(i+20)+1:p1(i+20)+p2(i+20)))
% end
% figure(4)
% for i = 1:4
%     subplot(2,2,i);
%     imagesc(phi_k{i+20,4}(:,p1(i+20)+p2(i+20)+1:p1(i+20)+p2(i+20)+p3(i+20)))
% end

%% orthogonality
ortho = cell(2,2);
for i = 1:nhub
    ortho{i,1}=phi_k{i,4}(:,1:p1(i)).'*phi_k{i,4}(:,p1(i)+1:p1(i)+p2(i));
    ortho{i,2}=phi_k{i,4}(:,p1(i)+1:p1(i)+p2(i)).'*phi_k{i,4}(:,p1(i)+p2(i)+1:p1(i)+p2(i)+p3(i));
    ortho{i,3}=phi_k{i,4}(:,1:p1(i)).'*phi_k{i,4}(:,p1(i)+p2(i)+1:p1(i)+p2(i)+p3(i));
    ortho{i,4}=phi_k{i,4}(:,1:p1(i)).'*phi_k{i,4}(:,1:p1(i));
    ortho{i,5}=phi_k{i,4}(:,p1(i)+1:p1(i)+p2(i)).'*phi_k{i,4}(:,p1(i)+1:p1(i)+p2(i));
    ortho{i,6}=phi_k{i,4}(:,p1(i)+p2(i)+1:p1(i)+p2(i)+p3(i)).'*phi_k{i,4}(:,p1(i)+p2(i)+1:p1(i)+p2(i)+p3(i));
    ortho{i,7}=phi_k{i,4}(:,1:p1(i)).'*Phi_ave;
    ortho{i,8}=phi_k{i,4}(:,p1(i)+1:p1(i)+p2(i)).'*Phi_ave;
    ortho{i,9}=phi_k{i,4}(:,p1(i)+p2(i)+1:p1(i)+p2(i)+p3(i)).'*Phi_ave;
end

% figure(5)
% for i = 1:9*nhub
%     subplot(nhub,9,i)
%     line = mod(i-1,9)+1;
%     row = (i - line)/9+1;
%     imagesc(ortho{row,line})
% end