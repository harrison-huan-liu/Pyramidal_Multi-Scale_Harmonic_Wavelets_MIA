function [AD_data,SMC_data,LMCI_data,EMCI_data,CN_data] = Proprocess_MMSE_AmyloidData(MMSE_path,Thickness_path)
[~,~,Thickness_RAW]=xlsread(Thickness_path);
[~,~,MMSE_RAW]=xlsread(MMSE_path);

PTID_MMSE=unique(MMSE_RAW(2:size(MMSE_RAW,1),3));
AD_count=1;
CN_count=1;
EMCI_count=1;
LMCI_count=1;
SMC_count=1;
AD_data.Thickness=zeros(148,2);
AD_data.SubjectID=cell(2,14);
CN_data.Thickness=zeros(148,2);
CN_data.SubjectID=cell(2,14);
EMCI_data.Thickness=zeros(148,2);
EMCI_data.SubjectID=cell(2,14);
LMCI_data.Thickness=zeros(148,2);
LMCI_data.SubjectID=cell(2,14);
SMC_data.Thickness=zeros(148,2);
SMC_data.SubjectID=cell(2,14);
for i=1:size(PTID_MMSE)
    temp=MMSE_RAW(strcmp(MMSE_RAW(:,3),PTID_MMSE{i}),:);
    [Y, I]=min(str2double(temp(:,5)));
    Subject_Id=temp(I,2);
    mmse_vector=temp(I,15);
    temp_vector=Thickness_RAW(strcmp(Thickness_RAW(:,2),Subject_Id),5);
    disease=temp(I,6);
    if ~isempty(temp_vector)
        switch char(disease)
            case 'AD'
                AD_data.MMSE(AD_count,:)=str2double(mmse_vector);
                AD_data.Thickness(:,AD_count)=str2double(temp_vector);
                AD_data.SubjectID(AD_count,:)=temp(I,:);
                AD_count=AD_count+1;
            case 'CN'
                CN_data.MMSE(CN_count,:)=str2double(mmse_vector);
                CN_data.Thickness(:,CN_count)=str2double(temp_vector);
                CN_data.SubjectID(CN_count,:)=temp(I,:);
                CN_count=CN_count+1;            
            case 'EMCI'
                EMCI_data.MMSE(EMCI_count,:)=str2double(mmse_vector);
                EMCI_data.Thickness(:,EMCI_count)=str2double(temp_vector);
                EMCI_data.SubjectID(EMCI_count,:)=temp(I,:);
                EMCI_count=EMCI_count+1;            
            case 'LMCI'
                LMCI_data.MMSE(LMCI_count,:)=str2double(mmse_vector);
                LMCI_data.Thickness(:,LMCI_count)=str2double(temp_vector);
                LMCI_data.SubjectID(LMCI_count,:)=temp(I,:);
                LMCI_count=LMCI_count+1;            
            case 'SMC' 
                SMC_data.MMSE(SMC_count,:)=str2double(mmse_vector);
                SMC_data.Thickness(:,SMC_count)=str2double(temp_vector);
                SMC_data.SubjectID(SMC_count,:)=temp(I,:);
                SMC_count=SMC_count+1;            
        end
    end
end
end