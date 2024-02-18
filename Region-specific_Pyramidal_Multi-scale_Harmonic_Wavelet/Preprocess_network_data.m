function Graph = Preprocess_network_data(path1,path2,p)
[~,~,RAW]=xlsread(path1);
Data_profile=cell(10,13);

k=1;% The number of baseline brain network 
i=2;
while i<size(RAW,1)
    j=i;
    Status=1;
    while strcmp(RAW{i,2},RAW{j,2})
       if strcmp(RAW{j,4},'bl')
           Data_profile(k,:)=RAW(j,:);
           Status=0;
       end
       j=j+1;
    end
    if Status
        Data_profile(k,:)=RAW(i,:);
    end
    i=j;
    k=k+1;
end
Graph=struct;
GroupNum=1;
for i=2:size(Data_profile,1)
    str=[path2,Data_profile{i,1},'_fdt_network_matrix'];
    temp = load(str);
    if ~isempty(temp)
        temp(temp<0.007)=0;
        v=sum(temp,2);
        if sum(v==0)==0
            D=diag(v); 
            temp=D^-1*temp;
            Graph(GroupNum).W=(temp+temp')/2; %Adjacency matrix
            Graph(GroupNum).D=diag(sum(Graph(GroupNum).W,2)); % Degree matrix
            Graph(GroupNum).L=Graph(GroupNum).D-Graph(GroupNum).W;% Laplacian matrix
            [Phi_temp,value]=eig(Graph(GroupNum).L);
            if sum(Phi_temp(:,1))<0
                Phi_temp=-Phi_temp;
            end
            Graph(GroupNum).Phi{1}=Phi_temp(:,1:p);
            Graph(GroupNum).Eigenvalue=value(1:p,1:p);
            Graph(GroupNum).SubjectID=Data_profile{i,1};
            Graph(GroupNum).PTID=Data_profile{i,2};
            Graph(GroupNum).VISCODE=Data_profile{i,4};
            Graph(GroupNum).DX_bl=Data_profile{i,5};
            Graph(GroupNum).Age=Data_profile{i,6};
            Graph(GroupNum).Gender=Data_profile{i,7};
            GroupNum=GroupNum+1;
        end
    end
end
end