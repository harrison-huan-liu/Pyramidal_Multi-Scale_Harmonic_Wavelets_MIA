function [accuracy,sensitivity,specificity,f_score,Beta] = SVM_signaldata(CN_signdata,AD_signdata,sample_num)
%% 10-fold
CN_num = size(CN_signdata,2);
AD_num = size(AD_signdata,2);
% right_num = zeros(10,1);
Beta = zeros(size(AD_signdata,1),10);

index_rand = randperm(sample_num);
total_signdata = [CN_signdata,AD_signdata];
total_signdata_rand = total_signdata(:,index_rand);
CN_signlabel = ones(CN_num,1);
AD_signlabel = zeros(AD_num,1);
total_signlabel = [CN_signlabel;AD_signlabel];
total_signlabel_rand = total_signlabel(index_rand,1);

train_data = total_signdata_rand.';
train_label = total_signlabel_rand;
CVSVM_model = fitcsvm(train_data,train_label,'KFold',10,'Standardize',true,'KernelFunction','linear');% ,'KernelScale','auto','Solver','L1QP'
ScoreCVSVMModel = fitSVMPosterior(CVSVM_model);
[labelPred,scorePred] = kfoldPredict(ScoreCVSVMModel);
for loop = 1:10
    Beta(:,loop) = ScoreCVSVMModel.Trained{loop,1}.Beta;
end
[cm,~] = confusionmat(train_label,labelPred);
accuracy = (cm(1,1)+cm(2,2))/(cm(1,1)+cm(1,2)+cm(2,1)+cm(2,2));
sensitivity = cm(1,1)/(cm(1,1)+cm(1,2));
specificity = cm(2,2)/(cm(2,2)+cm(2,1));
precision = cm(1,1)/(cm(1,1)+cm(2,1));
recall = cm(1,1)/(cm(1,1)+cm(1,2));
f_score = 2*precision*recall/(precision+recall);
if isnan(f_score)
    f_score = 0;
end
%     for i = 1:size(test_data,1)
%         if label(i)==test_label(i)
%             right_num(loop) = right_num(loop) + 1;
%         end
%     end
end