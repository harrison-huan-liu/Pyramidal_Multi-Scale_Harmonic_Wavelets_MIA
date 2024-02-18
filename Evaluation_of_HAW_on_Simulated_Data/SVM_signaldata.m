function [right_num,Beta] = SVM_signaldata(CN_signdata,AD_signdata,sample_num)
%% 10-fold
right_num = zeros(10,1);
for loop = 1:10
    start_num = (loop-1)*sample_num/10+1;
    end_num = loop*sample_num/10;
    CN_signdata_test = CN_signdata(:,start_num:end_num);
    AD_signdata_test = AD_signdata(:,start_num:end_num);
    test_data = [AD_signdata_test CN_signdata_test];
    test_data = test_data.';
    test_label = zeros(sample_num/5,1);
    for i = 1:sample_num/10
        test_label(i) = 1;
    end
    
    CN_signdata_delete = CN_signdata;
    AD_signdata_delete = AD_signdata;
    CN_signdata_delete(:,start_num:end_num) = [];
    AD_signdata_delete(:,start_num:end_num) = [];
    CN_signdata_train = CN_signdata_delete;
    AD_signdata_train = AD_signdata_delete;
    train_data = [AD_signdata_train CN_signdata_train];
    train_data = train_data.';
    train_label = zeros(sample_num*9/5,1);
    for i = 1:sample_num*9/10
        train_label(i) = 1;
    end
    SVM_model = fitcsvm(train_data,train_label,'Standardize',true,'KernelFunction','linear','KernelScale','auto','Solver','L1QP');
    Beta = SVM_model.Beta;
    [label,~] = predict(SVM_model,test_data);
    for i = 1:sample_num/5
        if label(i)==test_label(i)
            right_num(loop) = right_num(loop) + 1;
        end
    end
end
end