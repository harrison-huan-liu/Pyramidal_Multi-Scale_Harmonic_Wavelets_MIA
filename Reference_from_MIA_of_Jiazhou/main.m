close all;clear;clc
%% load brain connectivity network
load('harmonics_all.mat')
Thresholding_Code

for i = 1:size(Graph,2)
    W = Graph(i).W;
    deg = sum(W,1);
    flg = find(deg == 0);
    if ~isempty(flg)
        disp(['Data is error']); 
    end
end
    
nnn = 0;
for k = 1:size(Graph,2)
    graphcol = sum(Graph(k).W,2);
    for i = 1:264
        if graphcol(i) == 0
            nnn = nnn + 1;
        end
    end
end

p=55;% p is the number of eignvectors %55
for i = 1:size(Graph(1).W,1)
    u_vector{i} = CommonHarWavelets(i).Region_mask;
end

Graph=Preprocess_network_data('Data\DataTS.csv','Data\AD-Data\',p);
%% Estimating global common harmonic waves
fprintf('Estimating global common harmonic waves!\n');
CommonHarmonics=Estimate_glob_com_harmonics(Graph,p);
%% Constructing region-adaptive individual harmonic wavelets
fprintf('Constructing region-adaptive individual harmonic wavelets!\n');
Graph=Construct_individual_harmonic_wavelets(CommonHarmonics,Graph);
%% Identifying region-adaptive common harmonic wavelets
fprintf('Identifying region-adaptive common harmonic wavelets!\n');
CommonHarWavelets = cell(15,1);
CommonHarWaveletsLOSS = zeros(264,15);
for k = 1:30
    CommonHarWavelets{k}=Identify_glob_com_har_wavelets(Graph,u_vector,k);
    for i = 1:264
        CommonHarWaveletsLOSS(i,k)=norm(CommonHarWavelets{k}(i).Harmonics*(CommonHarWavelets{k}(i).Harmonics)',1)/sqrt(k);
    end
end
ddd = mean(CommonHarWaveletsLOSS);
for i = 1:19
    ddd(i) = ddd(i)*sqrt(i);
end
fprintf('Finished!\n')

%% save
save('CommonHarmonics.mat','CommonHarmonics');
save('CommonHarWavelets.mat','CommonHarWavelets');
save('Graph.mat','Graph');