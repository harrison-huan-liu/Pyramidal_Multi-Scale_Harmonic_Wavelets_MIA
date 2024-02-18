function [accuracy_my,sensitivity,specificity,f_score,Beta]=Classification_original(CN,LMCI,sample_num)
Objectnums_CN=size(CN,2);
Objectnums_LMCI=size(LMCI,2);
Objectnums_more=max(Objectnums_CN,Objectnums_LMCI);
Objectnums_less=min(Objectnums_CN,Objectnums_LMCI);
select_CN_LMCI_same=randperm(Objectnums_more,Objectnums_less);
Objectnums_reduce=2*Objectnums_less;

Objectnums=Objectnums_CN+Objectnums_LMCI;
total_data=cell(2,1);

Objectnums_more=max(Objectnums_CN,Objectnums_LMCI);
Objectnums_less=min(Objectnums_CN,Objectnums_LMCI);
% find the big and the small number of Objectnums_CN and Objectnums_LMCI
% select_CN_LMCI_same=randperm(Objectnums_more,Objectnums_less);

if Objectnums_less==Objectnums_CN
    for i=1:Objectnums_less
        total_data{i,1}=CN(1:size(CN,1),i);
        total_data{i+Objectnums_less,1}=LMCI(1:size(CN,1),select_CN_LMCI_same(i));
    end
else
    for i=1:Objectnums_less
        total_data{i,1}=CN(1:size(CN,1),select_CN_LMCI_same(i));
        total_data{i+Objectnums_less,1}=LMCI(1:size(CN,1),i);
    end
end

Objectnums_reduce=2*Objectnums_less;
label=zeros(Objectnums_reduce,1);
% select=randperm(Objectnums_reduce,Objectnums_reduce);

%% SVM Classification
% select_data=select;
% data=cell(2,1);
for i=1:Objectnums_reduce
%     data{i,1}=total_data{select_data(1,i),1};
    if i>Objectnums_less% select_data(1,i)
        label(i,1)=1;
    else
        label(i,1)=-1;
    end
end
data_double=cell2mat(total_data');

% disp(Objectnums_CN);

% % train_num=round(4*Objectnums_reduce/5);
% % test_num=Objectnums_reduce-train_num;
% % 
% % CN_sample=0;
% % LMCI_sample=0;
% % dele=0;
% % mm=1;
% % step_sample=zeros(test_num,1);
% % del_sample=zeros(1,1);
% % while CN_sample+LMCI_sample<test_num
% %     if select_data(1,mm)>Objectnums_less&&LMCI_sample<floor(test_num/2)
% %         LMCI_sample=LMCI_sample+1;
% %         step_sample(mm-dele,1)=mm;
% %     elseif select_data(1,mm)<=Objectnums_less&&CN_sample<test_num-floor(test_num/2)
% %         CN_sample=CN_sample+1;
% %         step_sample(mm-dele,1)=mm;
% %     else
% %         dele=dele+1;
% %         del_sample(dele,1)=mm;
% %     end
% %     mm=mm+1;
% % end
% % test_data=zeros(148,1);
% % test_label=zeros(1,1);
% % train_data=zeros(148,1);
% % train_label=zeros(1,1);
% % for i=1:test_num
% %     test_data(:,i)=data_double(:,step_sample(i,1));
% %     test_label(i,:)=label(step_sample(i,1),:);
% % end
% % if dele~=0
% %     for i=1:size(del_sample,1)% Objectnums_reduce
% %         train_data(:,i)=data_double(:,del_sample(i,1));
% %         train_label(i,:)=label(del_sample(i,1),:);
% %     end
% %     train_data(:,size(del_sample,1)+1:size(del_sample,1)+1+Objectnums_reduce-mm)=data_double(:,mm:Objectnums_reduce);
% %     train_label(size(del_sample,1)+1:size(del_sample,1)+1+Objectnums_reduce-mm,:)=label(mm:Objectnums_reduce,:);
% % else
% %     train_data(:,size(del_sample,1):size(del_sample,1)+Objectnums_reduce-mm)=data_double(:,mm:Objectnums_reduce);
% %     train_label(size(del_sample,1):size(del_sample,1)+Objectnums_reduce-mm,:)=label(mm:Objectnums_reduce,:);
% % end
% train_data=data_double(:,1:train_num);
% test_data=data_double(:,train_num+1:Objectnums_reduce);
% train_label=label(1:train_num,:);
% test_label=label(train_num+1:Objectnums_reduce,:);

% [auc,aupr]=Linear_SVM(train_data',test_data',train_label,test_label);

[accuracy_my,sensitivity,specificity,f_score,Beta]=Linear_SVM(data_double',label);

accuracy_my = accuracy_my(50);
sensitivity = sensitivity(50);
specificity = specificity(50);
f_score = f_score(50);
end
