function [AD_data,CN_data,EMCI_data,LMCI_data,SMC_data] = Proprocess_original_Tau_data(path)
% mkdir results
Tau_RAW=cell(10,159);
[~,~,Tau_temp]=xlsread(path);
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

AD_count=1;
CN_count=1;
EMCI_count=1;
LMCI_count=1;
SMC_count=1;
AD_data.Tau=zeros(148,2);
AD_data.SubjectID=cell(2,1);
CN_data.Tau=zeros(148,2);
CN_data.SubjectID=cell(2,1);
EMCI_data.Tau=zeros(148,2);
EMCI_data.SubjectID=cell(2,14);
LMCI_data.Tau=zeros(148,2);
LMCI_data.SubjectID=cell(2,1);
SMC_data.Tau=zeros(148,2);
SMC_data.SubjectID=cell(2,1);

for i=1:size(Tau_RAW,1)
    disease=Tau_RAW{i,4};
    if ~isempty(disease)
        temp_vector=Tau_RAW(i,12:159);
        switch char(disease)
            case 'AD'
                AD_data.Age(AD_count,:)=Tau_RAW(i,5);
                AD_data.Gender(AD_count,:)=Tau_RAW(i,6);
                AD_data.Tau(:,AD_count)=cell2mat(temp_vector)';
                AD_data.SubjectID(AD_count,:)=Tau_RAW(i,1);
                AD_count=AD_count+1;
            case 'CN'
                CN_data.Age(CN_count,:)=Tau_RAW(i,5);
                CN_data.Gender(CN_count,:)=Tau_RAW(i,6);
                CN_data.Tau(:,CN_count)=cell2mat(temp_vector)';
                CN_data.SubjectID(CN_count,:)=Tau_RAW(i,1);
                CN_count=CN_count+1;            
            case 'EMCI'
                EMCI_data.Age(EMCI_count,:)=Tau_RAW(i,5);
                EMCI_data.Gender(EMCI_count,:)=Tau_RAW(i,6);
                EMCI_data.Tau(:,EMCI_count)=cell2mat(temp_vector);
                EMCI_data.SubjectID(EMCI_count,:)=Tau_RAW(i,1);
                EMCI_count=EMCI_count+1;            
            case 'LMCI'
                LMCI_data.Age(LMCI_count,:)=Tau_RAW(i,5);
                LMCI_data.Gender(LMCI_count,:)=Tau_RAW(i,6);
                LMCI_data.Tau(:,LMCI_count)=cell2mat(temp_vector)';
                LMCI_data.SubjectID(LMCI_count,:)=Tau_RAW(i,1);
                LMCI_count=LMCI_count+1;            
            case 'SMC'
                SMC_data.Age(SMC_count,:)=Tau_RAW(i,5);
                SMC_data.Gender(SMC_count,:)=Tau_RAW(i,6);
                SMC_data.Tau(:,SMC_count)=cell2mat(temp_vector)';
                SMC_data.SubjectID(SMC_count,:)=Tau_RAW(i,1);
                SMC_count=SMC_count+1;            
        end
    end
end

end