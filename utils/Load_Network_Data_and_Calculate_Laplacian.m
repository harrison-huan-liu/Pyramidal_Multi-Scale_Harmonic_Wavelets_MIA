function [Graph,CommonNetwork,LatentLaplacian,Phi_ave,SubjectNum,NodeNum] = Load_Network_Data_and_Calculate_Laplacian(filter_label)
close all;clear;clc

%% load brain connectivity network
% the number of eignvectors % 55
p = 55;
% there are matrices of adjacency, degree, Laplacian and its corresponding eignvectors of brain connectivity network in var Graph
if filter_label
    Graph = Preprocess_data_np_noTransM_01('../Data/DataTS_Label_SMC_to_LMCI.xlsx','../Data/AD-Data/',p);
else
    Graph = Preprocess_network_data('../Data/DataTS.csv','../Data/AD-Data/',p);
end
% the number of subjects
SubjectNum = size(Graph,2);
% the number of node in each subject
NodeNum = size(Graph(1).L,1);

%% calculate European avarage common network
CommonNetwork=zeros(NodeNum);
for i = 1:SubjectNum
    CommonNetwork=CommonNetwork+Graph(i).W;
end
CommonNetwork=CommonNetwork/SubjectNum;
CommonNetwork(CommonNetwork<0.02)=0;
CommonNetwork=(CommonNetwork+CommonNetwork')/2;

%% calculate the Laplacian matrix of the European avarage common network
temp_D=diag(sum(CommonNetwork,2));
LatentLaplacian=temp_D-CommonNetwork;
[Phi_temp,~]=eig(LatentLaplacian);
Phi_ave=Phi_temp(:,1:p);
end