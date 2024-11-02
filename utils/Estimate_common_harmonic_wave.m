function [Graph,LatentLaplacian,CommonHarmonics] = Estimate_common_harmonic_wave(p,m,Graph,MST,CommonNetwork)
temp_D = diag(sum(CommonNetwork,2));
LatentLaplacian = temp_D-CommonNetwork;
[Phi_temp,~] = eig(LatentLaplacian);
Phi_ave = Phi_temp(:,1:p);

for i = 1:m
    Graph(i).W = MST{i};
    Graph(i).D = diag(sum(Graph(i).W,2));
    Graph(i).L = Graph(i).D-Graph(i).W;
    [Phi_temp,value] = eig(Graph(i).L);
    if sum(Phi_temp(:,1))<0
        Phi_temp = -Phi_temp;
    end
    Graph(i).Phi{1} = Phi_temp(:,1:p);
    Graph(i).Eigenvalue = value(1:p,1:p);
end

%% Estimating global common harmonic waves
CommonHarmonics=Phi_ave;
end