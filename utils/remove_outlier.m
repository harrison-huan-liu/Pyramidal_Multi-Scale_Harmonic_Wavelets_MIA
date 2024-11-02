function [Energy,outlier_num]=remove_outlier(Energy)
    Q1=prctile(Energy, 25);
    Q3=prctile(Energy, 75);

    low_bound=Q1-1.5*(Q3-Q1);
    up_bound=Q3+1.5*(Q3-Q1);

    len=size(Energy,1);
    outlier_num=zeros(2,1);
    j=1;
    for i=1:len
        if Energy(i)>=up_bound||Energy(i)<=low_bound
            outlier_num(j,1)=i;
            j=j+1;
        end
    end

    len_outlier=size(outlier_num,1);
    for i=1:len_outlier
        Energy(outlier_num(i)-i+1)=[];
    end