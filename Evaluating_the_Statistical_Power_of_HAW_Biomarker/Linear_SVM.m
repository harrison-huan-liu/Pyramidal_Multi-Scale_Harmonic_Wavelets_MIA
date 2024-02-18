function [accuracy_my,sensitivity,specificity,f_score,Beta]=Linear_SVM(X_data_b,Y_data_b)
lenn=size(Y_data_b,1);
wid=size(X_data_b,2);
Beta = zeros(size(X_data_b,2),10);
predict_label=zeros(99,lenn);
X_data=zeros(lenn,wid);
Y_data=zeros(lenn,1);
x=1;
y=1;
xx=2;
yy=2;
for i=1:lenn
    if i<=lenn/2
        X_data(x,:)=X_data_b(i,:);
        Y_data(y)=Y_data_b(i);
        x=x+2;
        y=y+2;
    else
        X_data(xx,:)=X_data_b(i,:);
        Y_data(yy)=Y_data_b(i);
        xx=xx+2;
        yy=yy+2;
    end
end
for loop=1:10
    train_num=floor(9*lenn/10);
    test_num=lenn-train_num;
    X_train=zeros(2,wid);
    Y_train=zeros(2,1);
    X_test=zeros(2,wid);
    Y_test=zeros(2,1);
    m=1;
    n=1;
    for i=1:lenn
        if i<=(loop-1)*test_num+test_num&&i>=(loop-1)*test_num+1
            X_test(m,:)=X_data(i,:);
            Y_test(m,:)=Y_data(i,:);
            m=m+1;
        else
            X_train(n,:)=X_data(i,:);
            Y_train(n,:)=Y_data(i,:);
            n=n+1;
        end
    end
    SVMModel = fitcsvm(X_train,Y_train,'KernelFunction','linear','Standardize',true,'Solver','L1QP');
    CompactSVMModel = compact(SVMModel);
    % [aucModel,aucParameters] = fitPosterior(SVMModel,'Holdout',0.10);
    [aucModel,aucParameters] = fitPosterior(CompactSVMModel,X_test,Y_test);
    Beta(:,loop) = aucModel.Beta;
    % [label,postProbs] = resubPredict(aucModel);
    [label,postProbs] = predict(aucModel,X_test);
    nn=size(X_test,1);
    for i=1:99
        for j=1:nn
            if postProbs(j,1)>=i/100
                predict_label(i,j+(loop-1)*test_num)=-1;
            else
                predict_label(i,j+(loop-1)*test_num)=1;
            end
        end
    end
end
TP=zeros(99,1);
TN=zeros(99,1);
FP=zeros(99,1);
FN=zeros(99,1);
for j=1:99
    for i=1:lenn
        if Y_data(i,1)==1&&predict_label(j,i)==1
            TP(j,1)=TP(j,1)+1;
        elseif Y_data(i,1)==-1&&predict_label(j,i)==-1
            TN(j,1)=TN(j,1)+1;
        elseif Y_data(i,1)==-1&&predict_label(j,i)==1
            FP(j,1)=FP(j,1)+1;
        elseif Y_data(i,1)==1&&predict_label(j,i)==-1
            FN(j,1)=FN(j,1)+1;
        end
    end
end
TPR=zeros(99,1);
FPR=zeros(99,1);
P=zeros(99,1);
% R=zeros(99,1);
for i=1:99
    TPR(i,1)=TP(i,1)/(TP(i,1)+FN(i,1));
    FPR(i,1)=FP(i,1)/(FP(i,1)+TN(i,1));
    P(i,1) = TP(i,1)/(TP(i,1)+FP(i,1));
end
X=FPR;
Y=TPR;
R=TPR;
for i=1:99
    if isnan(P(i,1))
        P(i,1)=1;
    end
end
accuracy_my=zeros(99,1);
sensitivity=zeros(99,1);
recall=zeros(99,1);
specificity=zeros(99,1);
f_score=zeros(99,1);
for i = 1:99
    accuracy_my(i)=(TP(i)+TN(i))/lenn;
    sensitivity(i)=TP(i)/(TP(i)+FP(i));
    recall(i)=TP(i)/(TP(i)+FN(i));
    specificity(i)=TN(i)/(TN(i)+FN(i));
    f_score(i)=2*sensitivity(i)*recall(i)/(sensitivity(i)+recall(i));
end
% [X,Y,T,auc] = perfcurve(resp, score(:,2), 'true');
auc=X(1)*Y(1)/2;
for i=2:size(X,1)
    auc=auc+(X(i)-X(i-1))*(Y(i)+Y(i-1))/2;
end
auc=auc+(1-X(size(X,1)))*(1+Y(size(X,1)))/2;
aupr=R(1)*(P(1)+1)/2;
% figure()
% plot(X,Y)
for i=2:size(P,1)
    aupr=aupr+(R(i)-R(i-1))*(P(i)+P(i-1))/2;
end
% figure()
% plot(R,P)
end