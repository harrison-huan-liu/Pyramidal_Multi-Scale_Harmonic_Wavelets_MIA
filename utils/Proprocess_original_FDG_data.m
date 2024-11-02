function [AD_data,CN_data,EMCI_data,LMCI_data,SMC_data] = Proprocess_original_FDG_data(path)
% mkdir results
FDG_RAW=cell(10,159);
[~,~,FDG_temp]=xlsread(path);
k=1;
i=2;
while i<size(FDG_temp,1)
    j=i;
    while strcmp(FDG_temp{i,1},FDG_temp{j,1})&&j<size(FDG_temp,1)
       j=j+1;
    end
    FDG_RAW(k,:)=FDG_temp(i,:);
    i=j;
    k=k+1;
end

AD_count=1;
CN_count=1;
EMCI_count=1;
LMCI_count=1;
SMC_count=1;
AD_data.FDG=zeros(148,2);
AD_data.SubjectID=cell(2,1);
CN_data.FDG=zeros(148,2);
CN_data.SubjectID=cell(2,1);
EMCI_data.FDG=zeros(148,2);
EMCI_data.SubjectID=cell(2,14);
LMCI_data.FDG=zeros(148,2);
LMCI_data.SubjectID=cell(2,1);
SMC_data.FDG=zeros(148,2);
SMC_data.SubjectID=cell(2,1);

for i=2:size(FDG_RAW,1)
    disease=FDG_RAW{i,4};
    if ~isempty(disease)
        temp_vector=FDG_RAW(i,12:159);
        switch char(disease)
            case 'AD'
                AD_data.FDG(:,AD_count)=cell2mat(temp_vector)';
                AD_data.SubjectID(AD_count,:)=FDG_RAW(i,1);
                AD_count=AD_count+1;
            case 'CN'
                CN_data.FDG(:,CN_count)=cell2mat(temp_vector)';
                CN_data.SubjectID(CN_count,:)=FDG_RAW(i,1);
                CN_count=CN_count+1;            
            case 'EMCI'
                EMCI_data.FDG(:,EMCI_count)=cell2mat(temp_vector);
                EMCI_data.SubjectID(EMCI_count,:)=FDG_RAW(i,1);
                EMCI_count=EMCI_count+1;            
            case 'LMCI'
                LMCI_data.FDG(:,LMCI_count)=cell2mat(temp_vector)';
                LMCI_data.SubjectID(LMCI_count,:)=FDG_RAW(i,1);
                LMCI_count=LMCI_count+1;            
            case 'SMC' 
                SMC_data.FDG(:,SMC_count)=cell2mat(temp_vector)';
                SMC_data.SubjectID(SMC_count,:)=FDG_RAW(i,1);
                SMC_count=SMC_count+1;            
        end
    end
end

end