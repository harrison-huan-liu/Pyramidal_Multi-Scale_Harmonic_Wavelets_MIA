function [Harmonicwavelets,multiPhi,singlePhi] = Identify_RM_harmonic_wavelets(n,hub_nodes,LatentLaplacian,CommonHarmonics,u_vector)
%% identified paraemter
nhub = size(hub_nodes,1);
u1 = 0.8; %3; %0.2; %0.8; %180;
u2 = 0.7; %6; %0.2; %0.7; %120;
Q = 6;
select_modules = 3;

%% initialize phi_k1, phi_k2, phi_k3
Theta_k = cell(nhub,select_modules);
phi_k = cell(nhub,select_modules);

for i = 1:nhub
    for j = 1:select_modules
        Theta_k{i,j} = LatentLaplacian + u1.*diag(ones(size(u_vector{i,j}))-u_vector{i,j}) + u2.*(CommonHarmonics * CommonHarmonics.');
        [phi_k{i,j},~] = eig(Theta_k{i,j});
    end
end

%% select the number of wavelets
gama = cell(nhub,select_modules);
alpha = zeros(nhub,1);
Harmonicwavelets = cell(nhub,1);

for i = 1:nhub
    for j = 1:select_modules
        gama{i,j} = zeros(3*Q,3*Q);
        for k = 1:Q
            gama{i,j}(k+Q*(j-1),k+Q*(j-1)) = 1;
        end
        Harmonicwavelets{i} = [Harmonicwavelets{i} phi_k{i,j}(:,1:Q)];
    end
end

%% calculate harmonic wavelets
mask_item = cell(nhub,select_modules);
linear_item = cell(nhub,1);
A = cell(nhub,1);
e = cell(1,nhub);
f = zeros(100,nhub);
mask_item_value = cell(nhub,select_modules);
linear_item_value = cell(nhub,1);
iterationnum = zeros(nhub,1);

for i = 1:nhub
    A{i} = LatentLaplacian+u2.*(CommonHarmonics*CommonHarmonics.');
    for j = 1:select_modules
        A{i} = A{i} + u1.*diag(ones(size(u_vector{i,j}))-u_vector{i,j});
        mask_item{i,j} = diag(ones(size(u_vector{i,j}))-u_vector{i,j});
    end
    [~,eigenvalue]=eig(A{i});
    if max(max(eigenvalue))>alpha(i)
        alpha(i) = max(max(eigenvalue));
    end
    linear_item{i,1} = alpha(i).* ones(n) - LatentLaplacian - u2.*CommonHarmonics * CommonHarmonics.';
    
    linear_item_value{i,1} = zeros(100,1);
    looptime = 1;
    e{looptime,i} = 1000;
    linear_item_value{i,1}(looptime,1) = trace(Harmonicwavelets{i}.'*linear_item{i,1}*Harmonicwavelets{i});
    f(looptime,i) = linear_item_value{i,1}(looptime,1);
    for j = 1:select_modules
        mask_item_value{i,j} = zeros(100,1);
        mask_item_value{i,j}(looptime,1) = trace(gama{i,j}*Harmonicwavelets{i}.'*mask_item{i,j}*Harmonicwavelets{i});
        f(looptime,i) = f(looptime,i) - u1.*mask_item_value{i,j}(looptime,1);
    end
    while abs(e{looptime,i})>=0.01&&looptime<100000
        M_gra = linear_item{i,1} * Harmonicwavelets{i};
        for j = 1:select_modules
            M_gra = M_gra - u1.*mask_item{i,j}*Harmonicwavelets{i}*gama{i,j};
        end
        [U,S,V] = svd(M_gra,'econ');
        Harmonicwavelets{i} = U*V.';
        looptime = looptime + 1;
        linear_item_value{i,1}(looptime,1) = trace(Harmonicwavelets{i}.'*linear_item{i,1}*Harmonicwavelets{i});
        f(looptime,i) = linear_item_value{i,1}(looptime,1);
        for j = 1:select_modules
            mask_item_value{i,j}(looptime,1) = trace(gama{i,j}*Harmonicwavelets{i}.'*mask_item{i,j}*Harmonicwavelets{i});
            f(looptime,i) = f(looptime,i) - u1.*mask_item_value{i,j}(looptime,1);
        end
        e{looptime,i} = f(looptime,i) - f(looptime-1,i);
    end
    iterationnum(i,1) = looptime;
    fprintf('Now the code is running in the %d node; \n',i);
end

% save('results\Theta_k.mat','Theta_k')
% save('results\phi_k.mat','phi_k')
% save('results\gama.mat','gama')
% save('results\alpha.mat','alpha')
% save('results\Harmonicwavelets.mat','Harmonicwavelets')
% save('results\mask_item.mat','mask_item')
% save('results\linear_item.mat','linear_item')
% save('results\mask_item_value.mat','mask_item_value')
% save('results\linear_item_value.mat','linear_item_value')

multiPhi = Harmonicwavelets{1};
for i = 2:nhub
    multiPhi = [multiPhi,Harmonicwavelets{i}];
end

singlePhi = Harmonicwavelets{1}(:,1:Q);
for i = 2:nhub
    singlePhi = [singlePhi,Harmonicwavelets{i}(:,1:Q)];
end
% save('results\multiPhi.mat','multiPhi')
% save('results\singlePhi.mat','singlePhi')
end