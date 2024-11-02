function []=statistic_MRI_data()
male_cn_num=0;
male_smc_num=0;
male_emci_num=0;
male_lmci_num=0;
male_ad_num=0;
female_cn_num=0;
female_smc_num=0;
female_emci_num=0;
female_lmci_num=0;
female_ad_num=0;
male_age=zeros(77,1);
male_num=0;
female_age=zeros(61,1);
female_num=0;
for i=1:138
    if strcmp(Graph(i).Gender,'Male')
        male_num=male_num+1;
        male_age(male_num)=Graph(i).Age;
        if strcmp(Graph(i).DX_bl,'CN')
            male_cn_num=male_cn_num+1;
        elseif strcmp(Graph(i).DX_bl,'SMC')
            male_smc_num=male_smc_num+1;
        elseif strcmp(Graph(i).DX_bl,'EMCI')
            male_emci_num=male_emci_num+1;
        elseif strcmp(Graph(i).DX_bl,'LMCI')
            male_lmci_num=male_lmci_num+1;
        elseif strcmp(Graph(i).DX_bl,'AD')
            male_ad_num=male_ad_num+1;
        end
    else
        female_num=female_num+1;
        female_age(female_num)=Graph(i).Age;
        if strcmp(Graph(i).DX_bl,'CN')
            female_cn_num=female_cn_num+1;
        elseif strcmp(Graph(i).DX_bl,'SMC')
            female_smc_num=female_smc_num+1;
        elseif strcmp(Graph(i).DX_bl,'EMCI')
            female_emci_num=female_emci_num+1;
        elseif strcmp(Graph(i).DX_bl,'LMCI')
            female_lmci_num=female_lmci_num+1;
        elseif strcmp(Graph(i).DX_bl,'AD')
            female_ad_num=female_ad_num+1;
        end
    end
end
end