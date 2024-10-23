function TransformedPhi=Calculate_IndividualPhi(IndividualPhi,IndividualL,Phi,lambda)
NodeNum=size(IndividualL,1);
p=size(IndividualPhi,2);
Phi_temp=IndividualPhi;
L=IndividualL;
B=lambda*Phi;
E=eig(L);
alpha=E(NodeNum);
L_tilde=alpha*eye(NodeNum,NodeNum)-L;
iter1=1;
Step1Diff=1;
Step1ObjectiveFuncValue=zeros(2,1);
Step1ObjectiveFuncValue(iter1)=trace(Phi_temp'*L*Phi_temp)+lambda*trace((Phi_temp-Phi)'*(Phi_temp-Phi));
while Step1Diff>0.00001&&iter1<100            
    M=L_tilde*Phi_temp+B;
    [U,S,V]=svd(M);
    Phi_temp=U(:,1:p)*V';
    iter1=iter1+1;
    Step1ObjectiveFuncValue(iter1)=trace(Phi_temp'*L*Phi_temp)+lambda*trace((Phi_temp-Phi)'*(Phi_temp-Phi));
    Step1Diff=abs(Step1ObjectiveFuncValue(iter1)-Step1ObjectiveFuncValue(iter1-1));           
end
TransformedPhi=Phi_temp;
if sum(sum(1*(isnan(TransformedPhi))))>0
    display([num2str(sum(sum(1*(isnan(TransformedPhi)))))]);
    fprintf('TransformedPhi Nan!\n')
end
end