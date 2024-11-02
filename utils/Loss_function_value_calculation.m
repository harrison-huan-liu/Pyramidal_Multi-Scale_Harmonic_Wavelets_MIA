function [fixed_item,limited_item,mask_item,function_value]=Loss_function_value_calculation(fixed_item_grad,limited_item_grad,mask_item_grad,phi,fixed_item,limited_item,mask_item,function_value,i,iter_time,u1,u2)
    fixed_item{i,iter_time} = trace(phi.'*fixed_item_grad*phi);
    limited_item{i,iter_time} = trace(phi.'*limited_item_grad*phi);

    function_value{i,iter_time} = fixed_item{i,iter_time} - u2.*limited_item{i,iter_time};
    for j = 1:scale
        mask_item{i,iter_time,j} = trace(gama{i,j}*phi.'*mask_item_grad{j,1}*phi);

        function_value{i,iter_time} = function_value{i,iter_time} - ...
        u1.*mask_item{i,iter_time,j};
    end
end