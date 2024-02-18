function [AD_data,SMC_data,LMCI_data,EMCI_data,CN_data] = Proprocess_MMSE_TauData_1(MMSE_path,Tau_path)
% [~,~,Tau_RAW]=xlsread(Tau_path);

Tau_RAW=cell(10,159);
[~,~,Tau_temp]=xlsread(Tau_path);
k=1;
i=2;
while i<size(Tau_temp,1)
    j=i;
    while strcmp(Tau_temp{i,1},Tau_temp{j,1})&&j<size(Tau_temp,1)
       j=j+1;
    end
    Tau_RAW(k,:)=Tau_temp(i,:);
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
AD_data.Tau=zeros(148,2);
AD_data.SubjectID=cell(2,18);
AD_data.MMSE=zeros(2,1);
CN_data.Tau=zeros(148,2);
CN_data.SubjectID=cell(2,18);
CN_data.MMSE=zeros(2,1);
EMCI_data.Tau=zeros(148,2);
EMCI_data.SubjectID=cell(2,18);
EMCI_data.MMSE=zeros(2,1);
LMCI_data.Tau=zeros(148,2);
LMCI_data.SubjectID=cell(2,18);
LMCI_data.MMSE=zeros(2,1);
SMC_data.Tau=zeros(148,2);
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
    temp_vector=Tau_RAW(strcmp(Tau_RAW(:,1),PTID_MMSE{i}),12:159);
    disease=Tau_RAW(strcmp(Tau_RAW(:,1),PTID_MMSE{i}),4);
    disease=cell2mat(disease);
    if ~isempty(temp_vector)
        switch char(disease)
            case 'AD'
                AD_data.MMSE(AD_count,:)=cell2mat(mmse);
                AD_data.Tau(:,AD_count)=cell2mat(temp_vector)';
                AD_data.SubjectID(AD_count,:)=temp(I,:);
                AD_count=AD_count+1;
            case 'CN'
                CN_data.MMSE(CN_count,:)=cell2mat(mmse_J);
                CN_data.Tau(:,CN_count)=cell2mat(temp_vector);
                CN_data.SubjectID(CN_count,:)=temp(I,:);
                CN_count=CN_count+1;            
            case 'EMCI'
                EMCI_data.MMSE(EMCI_count,:)=cell2mat(mmse);
                EMCI_data.Tau(:,EMCI_count)=cell2mat(temp_vector);
                EMCI_data.SubjectID(EMCI_count,:)=temp(I,:);
                EMCI_count=EMCI_count+1;            
            case 'LMCI'
                LMCI_data.MMSE(LMCI_count,:)=cell2mat(mmse_J);
                LMCI_data.Tau(:,LMCI_count)=cell2mat(temp_vector);
                LMCI_data.SubjectID(LMCI_count,:)=temp(I,:);
                LMCI_count=LMCI_count+1;            
            case 'SMC' 
                SMC_data.MMSE(SMC_count,:)=cell2mat(mmse);
                SMC_data.Tau(:,SMC_count)=cell2mat(temp_vector);
                SMC_data.SubjectID(SMC_count,:)=temp(I,:);
                SMC_count=SMC_count+1;            
        end
    end
end
end