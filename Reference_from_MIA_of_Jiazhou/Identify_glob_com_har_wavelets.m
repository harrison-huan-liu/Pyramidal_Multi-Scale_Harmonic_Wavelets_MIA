function CommonHarWavelets=Identify_glob_com_har_wavelets(Graph,k,beta_1,beta_2)
%% Identify region-adaptive harmonic wavelets for each sample
% Constructing common network
SubjectNum=size(Graph,2);
NodeNum=size(Graph(1).W,1);
% NodeNum=40;
CommonNetwork=zeros(NodeNum);
for i=1:SubjectNum
    CommonNetwork=CommonNetwork+Graph(i).W;
end
CommonNetwork=CommonNetwork/SubjectNum;
CommonNetwork(CommonNetwork<0.002)=0;
CommonNetwork=(CommonNetwork+CommonNetwork')/2;
% % Initialize parameters 
% beta_1=10;
% beta_2=30;
% % k=10; % k is the number of eignvectors of common harmonic wavelets
gama=0.005;
GlobalComHarmonics=Graph(1).GlobalComHarmonics;
%% Identifying the region-adaptive common harmonic wavelets of each node.
CommonHarWavelets=struct;
% load u_vector
for N_i=1:NodeNum % 40
    % Run Algrithom: Update Common Harmonic wavelets     
    u_vec=zeros(NodeNum,1);
    [~,I_index]=maxk(CommonNetwork(:,N_i),9);
    u_vec(I_index)=1;
    u_vec(N_i)=1;
% N_i = 68;
%     u_vec=u_vector{N_i};
    v=diag(1-u_vec);
    Theta=beta_1*v+beta_2*(GlobalComHarmonics*GlobalComHarmonics');
   
    Phi_k=Graph(1).LocalizedRegion(N_i).Harmonics(:,1:k);
    err=1;       
    iter2=1;
    while err>0.00001&&iter2<100 % 0.0000001
        Phi_increment=zeros(NodeNum,k);
        for S_j=1:SubjectNum
            Phi_increment=Phi_increment+Phi_k*Graph(S_j).LocalizedRegion(N_i).Harmonics(:,1:k)'*Phi_k-Graph(S_j).LocalizedRegion(N_i).Harmonics(:,1:k);
        end
        Phi_increment=Phi_increment+2*(eye(NodeNum)-Phi_k*Phi_k')*Theta*Phi_k;
        Phi_increment=-gama*Phi_increment;
        [Q,R]=qr((eye(NodeNum)-Phi_k*Phi_k')*Phi_increment,0);
        A=Phi_k'*Phi_increment;
        BC=expm([A,-R';R,zeros(k)])*[eye(k);zeros(k)];        
        Phi_k=Phi_k*BC(1:k,:)+Q*BC(k+1:2*k,:);       
        err=norm(Phi_increment,'fro');
        iter2=iter2+1;
    end   
    CommonHarWavelets(N_i).Region_mask=u_vec;
    CommonHarWavelets(N_i).v=v;
    CommonHarWavelets(N_i).Harmonics=Phi_k;
    fprintf('node %d',N_i);
end   
end

function [value,idx]=maxk(x,k)
value = zeros(1,k);
idx   = zeros(1,k);
m = min(x);
for j = 1:k
    [value(j),idx(j)] = max(x);
    x(idx(j))=m;
end
end
