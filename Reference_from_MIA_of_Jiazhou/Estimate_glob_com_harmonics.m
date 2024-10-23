function CommonHarmonics=Estimate_glob_com_harmonics(Graph,p)
%% Initialize common Harmonics as the average of all subject's networks
SubjectNum=size(Graph,2);
NodeNum=size(Graph(1).W,1);
CommonNetwork=zeros(NodeNum);
for i=1:SubjectNum
    CommonNetwork=CommonNetwork+Graph(i).W;
end
CommonNetwork=CommonNetwork/SubjectNum;
CommonNetwork(CommonNetwork<0.02)=0;
CommonNetwork=(CommonNetwork+CommonNetwork')/2;
temp_D=diag(sum(CommonNetwork,2));
LatentLaplacian=temp_D-CommonNetwork;
[Phi_temp,~]=eig(LatentLaplacian);
Phi_ave=Phi_temp(:,1:p);

%% Estimating global common harmonic waves
% Initializing parameters
lambda_1=0.005;
gama=0.0004/(lambda_1); %0.005
iter=1;
Phi=cell(5,1);
Diff(1)=100;
ObjectiveFuncValue=zeros(2,1);
Phi{1}=Phi_ave;
% Runing Algrithom   
while Diff(iter)>0.1&&iter<200 % 0.1
    %Step1: Updating individual Phi
    for i=1:SubjectNum
        Graph(i).Phi{iter+1}=Calculate_IndividualPhi(Graph(i).Phi{iter},Graph(i).L,Phi{iter},lambda_1);     
    end

    %Step2: Updating Common Harmonics Phi    
    Phi_k=Graph(1).Phi{iter+1};
    err=1;
    iter2=1;
    Step2ObjectFunction=zeros(2,1);
    for i=1:SubjectNum
        Step2ObjectFunction(iter2)=Step2ObjectFunction(iter2)+lambda_1*(p-trace(Graph(i).Phi{iter+1}'*Phi_k));
    end
    while err>0.0001&&iter2<500 %0.0001
        Phi_increment=zeros(NodeNum,p);
        for i=1:SubjectNum
            Phi_increment=Phi_increment+Phi_k*Graph(i).Phi{iter+1}'*Phi_k-Graph(i).Phi{iter+1};
        end
        Phi_increment=-gama*lambda_1*Phi_increment;
        if sum(sum(1*(isnan(Phi_increment))))>0
            display([num2str(sum(sum(1*(isnan(Phi_increment)))))]);
            fprintf('Phi_increment Nan!\n')
        end
        [Q,R]=qr((eye(NodeNum)-Phi_k*Phi_k')*Phi_increment,0);
        A=Phi_k'*Phi_increment;
        BC=expm([A,-R';R,zeros(p)])*[eye(p);zeros(p)];        
        Phi{iter+1}=Phi_k*BC(1:p,:)+Q*BC(p+1:2*p,:);
        Phi_k=Phi{iter+1};        
        err=norm(Phi_increment,'fro');
        temp=0;
        for i=1:SubjectNum
            temp=temp+lambda_1*(p-trace(Graph(i).Phi{iter+1}'*Phi_k));
        end
        Step2ObjectFunction(iter2+1)=temp;
        iter2=iter2+1;
    end
    % Convergent condition
    temp=0;
    for i=1:SubjectNum
        temp=temp+trace(Graph(i).Phi{iter+1}'*Graph(i).L*Graph(i).Phi{iter+1})+lambda_1*(p-trace(Graph(i).Phi{iter+1}'*Phi{iter+1}));
    end
    ObjectiveFuncValue(iter+1)=temp;
    Diff(iter+1)=abs(ObjectiveFuncValue(iter+1)- ObjectiveFuncValue(iter));   
    fprintf('The iteration No. is %d, the error is %f ....\n', iter, Diff(iter+1));
    iter=iter+1;   
end
CommonHarmonics=Phi{iter};
end