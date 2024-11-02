function Graph = Preprocess_data_np_noTransM_01(path_subject_information,path_subject_matrices,p)
% mkdir results
%RAW=readcell(path_subject_information);
[~,~,RAW]=xlsread(path_subject_information);
Data_profile=RAW;
Graph=struct;
GroupNum=1;
for i=2:size(Data_profile,1)
    str=[path_subject_matrices,Data_profile{i,1},'_fdt_network_matrix'];%/home/ubuntu/tgh/TMI_and_MICCAI/
    temp = load(str);
    temp(temp<0.007)=0;
    v=sum(temp,2);
    if sum(v==0)==0
        D=diag(v); 
        temp=D^-1*temp;
%         temp(temp<0.007)=0;
%         temp(temp>=0.005)=1;
        Graph(GroupNum).W=(temp+temp')/2;
        Graph(GroupNum).D=diag(sum(Graph(GroupNum).W,2));
        Graph(GroupNum).L=Graph(GroupNum).D-Graph(GroupNum).W;
        [Phi_temp,~]=eig(Graph(GroupNum).L);
        if sum(Phi_temp(:,1))<0
            Phi_temp=-Phi_temp;
        end
        Graph(GroupNum).Phi{1}=Phi_temp(:,1:p);
        Graph(GroupNum).DX_bl=Data_profile{i,5};
        Graph(GroupNum).TransformedPhi{1}=Graph(GroupNum).Phi{1};
        GroupNum=GroupNum+1;
    end   
end
end