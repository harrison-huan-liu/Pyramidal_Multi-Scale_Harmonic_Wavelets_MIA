function [AD_data,SMC_data,LMCI_data,EMCI_data,CN_data] = Proprocess_MMSE_FDGData_1(MMSE_path,FDG_path)
% [~,~,FDG_RAW]=xlsread(FDG_path);

FDG_RAW=cell(10,159);
[~,~,FDG_temp]=xlsread(FDG_path);
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

[~,~,MMSE_RAW]=xlsread(MMSE_path);
PTID_MMSE=unique(MMSE_RAW(2:size(MMSE_RAW,1),3));

AD_count=1;
CN_count=1;
EMCI_count=1;
LMCI_count=1;
SMC_count=1;
AD_data.FDG=zeros(148,2);
AD_data.SubjectID=cell(2,18);
AD_data.MMSE=zeros(2,1);
CN_data.FDG=zeros(148,2);
CN_data.SubjectID=cell(2,18);
CN_data.MMSE=zeros(2,1);
EMCI_data.FDG=zeros(148,2);
EMCI_data.SubjectID=cell(2,18);
EMCI_data.MMSE=zeros(2,1);
LMCI_data.FDG=zeros(148,2);
LMCI_data.SubjectID=cell(2,18);
LMCI_data.MMSE=zeros(2,1);
SMC_data.FDG=zeros(148,2);
SMC_data.SubjectID=cell(2,18);
SMC_data.MMSE=zeros(2,1);

for i=1:size(PTID_MMSE)
    temp=MMSE_RAW(strcmp(MMSE_RAW(:,3),PTID_MMSE{i}),:);
    [~, I]=min(str2double(temp(:,4)));
    [~, JJ]=max(str2double(temp(:,4)));
    MMSE_exist=temp(I,15);
    while isnan(MMSE_exist{1})
        I=I+1;
        MMSE_exist=temp(I,15);
    end
    mmse=temp(I,15);
    MMSE_J_exist=temp(JJ,15);
    while isnan(MMSE_J_exist{1}) || (cell2mat(MMSE_J_exist)<15)
        JJ=JJ-1;
        MMSE_J_exist=temp(JJ,15);
    end
    mmse_J=temp(JJ,15);
    temp_vector=FDG_RAW(strcmp(FDG_RAW(:,1),PTID_MMSE{i}),12:159);
    disease=FDG_RAW(strcmp(FDG_RAW(:,1),PTID_MMSE{i}),4);
    disease=cell2mat(disease);
    if ~isempty(temp_vector)
        switch char(disease)
            case 'AD'
                AD_data.MMSE(AD_count,:)=cell2mat(mmse);
                AD_data.FDG(:,AD_count)=cell2mat(temp_vector)';
                AD_data.SubjectID(AD_count,:)=temp(I,:);
                AD_count=AD_count+1;
            case 'CN'
                CN_data.MMSE(CN_count,:)=cell2mat(mmse_J);
                CN_data.FDG(:,CN_count)=cell2mat(temp_vector);
                CN_data.SubjectID(CN_count,:)=temp(I,:);
                CN_count=CN_count+1;            
            case 'EMCI'
                EMCI_data.MMSE(EMCI_count,:)=cell2mat(mmse);
                EMCI_data.FDG(:,EMCI_count)=cell2mat(temp_vector);
                EMCI_data.SubjectID(EMCI_count,:)=temp(I,:);
                EMCI_count=EMCI_count+1;            
            case 'LMCI'
                LMCI_data.MMSE(LMCI_count,:)=cell2mat(mmse_J);
                LMCI_data.FDG(:,LMCI_count)=cell2mat(temp_vector);
                LMCI_data.SubjectID(LMCI_count,:)=temp(I,:);
                LMCI_count=LMCI_count+1;            
            case 'SMC' 
                SMC_data.MMSE(SMC_count,:)=cell2mat(mmse);
                SMC_data.FDG(:,SMC_count)=cell2mat(temp_vector);
                SMC_data.SubjectID(SMC_count,:)=temp(I,:);
                SMC_count=SMC_count+1;            
        end
    end
end
end