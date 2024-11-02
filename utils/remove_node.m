function Energy=remove_node(Energy,outlier_num)
len_outlier=size(outlier_num,1);
for i=1:len_outlier
    Energy(outlier_num(i)-i+1,:)=[];
end