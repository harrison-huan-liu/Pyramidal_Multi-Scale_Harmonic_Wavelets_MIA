function [accuracy,sensitivity,specificity,f_score,Beta] = SVM_signaldata(CN_signdata,AD_signdata,sample_num)
%% 10-fold
CN_num = size(CN_signdata,2);
AD_num = size(AD_signdata,2);
% right_num = zeros(10,1);
accuracy = zeros(10,1);
sensitivity = zeros(10,1);
specificity = zeros(10,1);
f_score = zeros(10,1);

Beta = zeros(size(AD_signdata,1),10);

index_rand = randperm(sample_num);
total_signdata = [CN_signdata,AD_signdata];
total_signdata_rand = total_signdata(:,index_rand);
CN_signlabel = ones(CN_num,1);
AD_signlabel = zeros(AD_num,1);
total_signlabel = [CN_signlabel;AD_signlabel];
total_signlabel_rand = total_signlabel(index_rand,1);
for loop = 1:10
    start_num = round((loop-1)*sample_num/10)+1;
    end_num = round(loop*sample_num/10);
    total_signdata_test = total_signdata_rand(:,start_num:end_num);
    test_data = total_signdata_test.';
    total_signlabel_test = total_signlabel_rand(start_num:end_num,1);
    test_label = total_signlabel_test;
    
    total_signdata_delete = total_signdata_rand;
    total_signlabel_delete = total_signlabel_rand;
    total_signdata_delete(:,start_num:end_num) = [];
    total_signlabel_delete(start_num:end_num,:) = [];
    total_signdata_train = total_signdata_delete;
    total_signlabel_train = total_signlabel_delete;
    train_data = total_signdata_train.';
    train_label = total_signlabel_train;
    SVM_model = fitcsvm(train_data,train_label,'Standardize',true,'KernelFunction','linear','KernelScale','auto','Solver','L1QP');
    Beta(:,loop) = SVM_model.Beta;
    [label,~] = predict(SVM_model,test_data);
    [cm,~] = confusionmat(test_label,label);
    accuracy(loop) = (cm(1,1)+cm(2,2))/(cm(1,1)+cm(1,2)+cm(2,1)+cm(2,2));
    sensitivity(loop) = cm(1,1)/(cm(1,1)+cm(1,2));
    specificity(loop) = cm(2,2)/(cm(2,2)+cm(2,1));
    precision = cm(1,1)/(cm(1,1)+cm(2,1));
    recall = cm(1,1)/(cm(1,1)+cm(1,2));
    f_score(loop) = 2*precision*recall/(precision+recall);
%     for i = 1:size(test_data,1)
%         if label(i)==test_label(i)
%             right_num(loop) = right_num(loop) + 1;
%         end
%     end
end
end