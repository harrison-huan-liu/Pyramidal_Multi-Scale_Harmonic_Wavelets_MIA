function Graph=Construct_individual_harmonic_wavelets(CommonHarmonics,Graph,miu_1,miu_2)
%% Constructing region-adaptive individual harmonic wavelets for each sample
SubjectNum=size(Graph,2);
NodeNum=size(Graph(1).W,1);
% NodeNum=40;
% Initializing parameters
% miu_1=10;
% miu_2=30;
Q=10; % Q is the number of eignvectors of individual harmonic wavelets

%% 20220727
% % load u_vector
% for i = 1:size(Graph(1).W,1)
%     u_vector{i} = CommonHarWavelets(i).Region_mask;
% end

for S_i=1:SubjectNum
    Graph(S_i).GlobalComHarmonics=CommonHarmonics;
    for N_j=1:NodeNum
        u_vec=zeros(NodeNum,1);
        [~,I_index]=maxk(Graph(S_i).W(:,N_j),9);
        u_vec(I_index)=1;
        u_vec(N_j)=1;
%         u_vec=u_vector{N_j};
        v=diag(1-u_vec);
%         Graph(S_i).LocalizedRegion(N_j).Region_mask=u_vec; % subnetwork mask
%         Graph(S_i).LocalizedRegion(N_j).v=v;
        Theta=Graph(S_i).L+miu_1*v+miu_2*(Graph(S_i).GlobalComHarmonics*Graph(S_i).GlobalComHarmonics');      
        [LocalEigenvector,LocalEigenvalue]=eig(Theta);
        Graph(S_i).LocalizedRegion(N_j).Harmonics=LocalEigenvector(:,1:Q);% individual harmonic wavelets
%         Graph(S_i).LocalizedRegion(N_j).Eigenvalue=diag(LocalEigenvalue(1:Q,1:Q)); % corresponding eigenvalue
    end
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