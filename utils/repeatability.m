clear
clc
% SVM_weight_lc=readtable('SVM_weight_lc.xls');
% SVM_weight_lc=load('results/SVM_para.mat');
SVM_weight_lc=load('results/tem14.mat');
select_num=zeros(10000,1);
for i=1:1000
    weight_single=mean(SVM_weight_lc.tem14{i,1},2);
    [~,select_num(1+10*(i-1):10*i,1)]=maxk(weight_single,10);
end
[M,F] = mode(select_num);
