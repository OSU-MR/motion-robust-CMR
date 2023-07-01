% Hermitian 3D Operator
function y = multTr(obj,x)

    % Upsample
    
    
%     % Put input on GPU
%     if strcmpi(obj.compute,'gpu_off')
%         x = gpuArray(x);
%     end
    
    % Use single
    if strcmpi(obj.precision,'single')
        if strcmpi(obj.compute,'gpu')
            y = zeros(obj.frame_size(1),obj.frame_size(2),obj.frame_size(3),size(obj.C,4),obj.Q,'single','gpuArray');
        else
            y = zeros(obj.frame_size(1),obj.frame_size(2),obj.frame_size(3),size(obj.C,4),obj.Q,'single');
        end
    else
        y = zeros(obj.frame_size(1),obj.frame_size(2),obj.frame_size(3),size(obj.C,4),obj.Q);
    end
    
%     % Put on GPU
%     if strcmpi(obj.compute,'gpu_off') || strcmpi(obj.compute,'gpu')
%         y = gpuArray(y);
%     end
    
    % Weigthed Least Squares
    if ~isempty(obj.mask_weights)
        x = obj.mask_weights.*x;
    end
    
    
    
    y(obj.mask_patterns) = x;
    % iFFT
    y = sqrt(size(y,1)*size(y,2)*size(y,3))*ifft(ifft(ifft(y,[],1),[],2),[],3);
    % Coil Multiply
    fun = @(x,y) x.*y;
    
    % Check if multiple sensitivity maps are used for espirit
    arg = bsxfun(@times, y, conj(permute(obj.C, [1,2,3,4,6,5])));
    
%     if size(obj.C,5) > 1
%         arg = cat(5,bsxfun(fun,y(:,:,:,:,1:obj.Q/2),conj(obj.C(:,:,:,:,1))),bsxfun(fun,y(:,:,:,:,(obj.Q/2)+1:obj.Q),conj(obj.C(:,:,:,:,2))));
%     else
%         arg = bsxfun(fun,y,conj(obj.C));
%     end
    
    y = squeeze(sum(arg,4));

    % Make a column vector
    y = y(:);
    
    % put array back on cpu
    if strcmpi(obj.compute,'gpu_off')
        y = gather(y);
    end
end