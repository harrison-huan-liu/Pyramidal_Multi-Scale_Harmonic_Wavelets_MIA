function [AD_data,SMC_data,LMCI_data,EMCI_data,CN_data] = Proprocess_MMSE_ThickData_1(Information_new_path,Thickness_path,Subject_path)
[~,~,Thickness_RAW]=xlsread(Thickness_path);
[~,~,Information_new_RAW]=xlsread(Information_new_path);
[~,~,Subject_RAW]=xlsread(Subject_path);

PTID=unique(Subject_RAW(2:size(Subject_RAW,1),3));
PTID_new=unique(Information_new_RAW(2:size(Information_new_RAW,1),3));
AD_count=1;
CN_count=1;
EMCI_count=1;
LMCI_count=1;
SMC_count=1;
AD_data.Thickness=zeros(148,2);
AD_data.SubjectID=cell(2,14);
AD_data.MMSE=zeros(2,1);
AD_data.MMSE_new=zeros(2,1);
CN_data.Thickness=zeros(148,2);
CN_data.SubjectID=cell(2,14);
CN_data.MMSE=zeros(2,1);
CN_data.MMSE_new=zeros(2,1);
EMCI_data.Thickness=zeros(148,2);
EMCI_data.SubjectID=cell(2,14);
EMCI_data.MMSE=zeros(2,1);
EMCI_data.MMSE_new=zeros(2,1);
LMCI_data.Thickness=zeros(148,2);
LMCI_data.SubjectID=cell(2,14);
LMCI_data.MMSE=zeros(2,1);
LMCI_data.MMSE_new=zeros(2,1);
SMC_data.Thickness=zeros(148,2);
SMC_data.SubjectID=cell(2,14);
SMC_data.MMSE=zeros(2,1);
SMC_data.MMSE_new=zeros(2,1);
for i=1:size(PTID)
    temp=Subject_RAW(strcmp(Subject_RAW(:,3),PTID{i}),:);
    [Y, I]=min(cell2mat(temp(:,5)));
    Subject_Id=temp(I,2);
    Subject_PTID=temp(I,3);
    temp_vector=Thickness_RAW(strcmp(Thickness_RAW(:,2),Subject_Id),5);
%     mmse=temp(I,14);
    mmse_new=Information_new_RAW(strcmp(Information_new_RAW(:,3),Subject_PTID),:);
    [Z, X]=min(cell2mat(mmse_new(:,4)));
    mmse_new_exist=mmse_new(X,15);
    while isempty(mmse_new_exist)
        X=X+1;
    end
    mmse_new_r=mmse_new(X,15);
    disease=temp(I,6);
    if ~isempty(temp_vector)
        switch char(disease)
            case 'AD'
%                 AD_data.MMSE(AD_count,:)=cell2mat(mmse);
                AD_data.MMSE_new(AD_count,:)=cell2mat(mmse_new_r);
                AD_data.Thickness(:,AD_count)=str2double(temp_vector);
                AD_data.SubjectID(AD_count,:)=temp(I,:);
                AD_count=AD_count+1;
            case 'CN'
%                 CN_data.MMSE(CN_count,:)=cell2mat(mmse);
                CN_data.MMSE_new(CN_count,:)=cell2mat(mmse_new_r);
                CN_data.Thickness(:,CN_count)=str2double(temp_vector);
                CN_data.SubjectID(CN_count,:)=temp(I,:);
                CN_count=CN_count+1;            
            case 'EMCI'
%                 EMCI_data.MMSE(EMCI_count,:)=cell2mat(mmse);
                EMCI_data.MMSE_new(EMCI_count,:)=cell2mat(mmse_new_r);
                EMCI_data.Thickness(:,EMCI_count)=str2double(temp_vector);
                EMCI_data.SubjectID(EMCI_count,:)=temp(I,:);
                EMCI_count=EMCI_count+1;            
            case 'LMCI'
%                 LMCI_data.MMSE(LMCI_count,:)=cell2mat(mmse);
                LMCI_data.MMSE_new(LMCI_count,:)=cell2mat(mmse_new_r);
                LMCI_data.Thickness(:,LMCI_count)=str2double(temp_vector);
                LMCI_data.SubjectID(LMCI_count,:)=temp(I,:);
                LMCI_count=LMCI_count+1;            
            case 'SMC' 
%                 SMC_data.MMSE(SMC_count,:)=cell2mat(mmse);
                SMC_data.MMSE_new(SMC_count,:)=cell2mat(mmse_new_r);
                SMC_data.Thickness(:,SMC_count)=str2double(temp_vector);
                SMC_data.SubjectID(SMC_count,:)=temp(I,:);
                SMC_count=SMC_count+1;            
        end
    end
end
end