function [fixed_item,limited_item,mask_item,phi_iter,identify_wavelets_orthogonal,loss,function_value,iteration_num,multiPhi,singlePhi]=Multiscale_harmonic_wavelets_calculation(nhub,scale,alpha,LatentLaplacian,NodeNum,Phi_ave,phi_k,u_vector,u1,u2,gama)
    % u1{i,j}=0.1*abs(fixed_item1{i,iter_time})/(abs(mask_item2{i,iter_time}));
    % u2{i,j}=100*abs(fixed_item1{i,iter_time})/(abs(fixed_item2{i,iter_time}));
    fixed_item_grad = alpha.* ones(NodeNum) - LatentLaplacian;
    limited_item_grad = Phi_ave * Phi_ave.';

    fixed_item = cell(nhub,2);
    limited_item = cell(nhub,2);
    mask_item = cell(nhub,2,scale);

    phi_iter = cell(nhub,2);
    identify_wavelets_orthogonal = cell(nhub,2);
    loss = cell(nhub,2);
    function_value = cell(nhub,2);
    iteration_num = zeros(nhub,1);

    multiPhi = [];
    singlePhi = [];

    wb_admm = parwaitbar(nhub,'WaitMessage','Calculating Optimal MultiPhi...','FinalMessage','Done!')
    for i = 1:nhub
        mask_item_grad = cell(scale,1);
        for j = 1:scale
            mask_item_grad{j,1} = diag(ones(size(u_vector{i,j}))-u_vector{i,j});
        end

        phi = phi_k{i,scale+1};
        iter_time = 1;
        phi_iter{i,iter_time} = phi;
        identify_wavelets_orthogonal{i,iter_time} = phi.'*phi;
        loss{i,iter_time} = 1000;

        [fixed_item,limited_item,mask_item,function_value]=Loss_function_value_calculation(fixed_item_grad,limited_item_grad,mask_item_grad,phi,fixed_item,limited_item,mask_item,function_value,i,iter_time,u1,u2)

        while abs(loss{i,iter_time})>=1
            M_gra = (fixed_item_grad - u2.*limited_item_grad) * phi;
            
            for j = 1:scale
                M_gra = M_gra - ...
                u1.*mask_item_grad{j,1}*phi*gama{i,j};
            end

            % % full SVD
            % [U,S,V] = svd(M_gra);
            % if q(i,1)>n
            %     V(:,149:q(i,1)) = [];
            %     phi = U(:,1:n)*V.';
            % else
            %     phi = U(:,1:q(i,1))*V.';
            % end
            % compact SVD
            [U,S,V] = svd(M_gra,'econ');
            phi = U*V.';

            iter_time = iter_time + 1;
            phi_iter{i,iter_time} = phi;
            identify_wavelets_orthogonal{i,iter_time} = phi.'*phi;

            [fixed_item,limited_item,mask_item,function_value]=Loss_function_value_calculation(fixed_item_grad,limited_item_grad,mask_item_grad,phi,fixed_item,limited_item,mask_item,function_value,i,iter_time,u1,u2)

            loss{i,iter_time} = function_value{i,iter_time} - function_value{i,iter_time-1};
        end
        iteration_num(i,1) = iter_time;
        multiPhi = [multiPhi,phi];
        singlePhi = [singlePhi,phi(:,1:6)];
        wb_admm.progress();
    end
end