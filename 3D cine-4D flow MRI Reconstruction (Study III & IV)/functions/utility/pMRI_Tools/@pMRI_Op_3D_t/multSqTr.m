% Squared-Hermitian-Transposed Matrix multiply 
function y = multSqTr(obj,x)

% Use Rank 1 approximation for |A^H|.^2
if obj.uniform_var == 1
    y = ones(obj.N,1)*(obj.Fro2*sum(x,1));

% Use true |A^H|.^2
else
    % Put input on GPU
    if strcmpi(obj.compute,'gpu_off') || strcmpi(obj.compute,'gpu')
        x = gpuArray(x);
        tmp = zeros(obj.frame_size(1),obj.frame_size(2),obj.frame_size(3),size(obj.C,4),obj.Q,'gpuArray');
    else
       tmp = zeros(obj.frame_size(1),obj.frame_size(2),obj.frame_size(3),size(obj.C,4),obj.Q);
    end
 
    % Weigthed Least Squares
    if ~isempty(obj.mask_weights)
        x = (obj.mask_weights.^2).*x;
    end
    
    % Upsample
    tmp(obj.mask_patterns) = x;
    
    % True Square transpose multiply
    tmp = sum(sum(sum(tmp,1),2),3);
    y = bsxfun(@times,obj.CSqTr,tmp);
    y = sum(y,4);
    y = y(:);

    % put array back on cpu
    if strcmpi(obj.compute,'gpu_off')
        y = gather(y);
    end
end

end

