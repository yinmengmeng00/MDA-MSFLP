function result = combine_kernels(weights, kernels)
%核融合（在获取各个核权重之后）
    % length of weights should be equal to length of matrices
    n = length(weights);
    result = zeros(size(kernels(:,:,1)));    
    
    for i=1:n
        result = result + weights(i) * kernels(:,:,i);
    end
end