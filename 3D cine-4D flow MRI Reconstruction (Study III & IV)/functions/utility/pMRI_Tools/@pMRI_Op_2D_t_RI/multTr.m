% Hermitian 3D Operator
function y = multTr(obj,x)
    % Put input on GPU
    if strcmpi(obj.compute,'gpu_off')
        x = gpuArray(x);
    end
    
    % Use single
    if strcmpi(obj.precision,'single')
        y = zeros(obj.frame_size(1),obj.frame_size(2),size(obj.C,3),obj.Q,'single');
    else
        y = zeros(obj.frame_size(1),obj.frame_size(2),size(obj.C,3),obj.Q);
    end
    
    % Put on GPU
    if strcmpi(obj.compute,'gpu_off') || strcmpi(obj.compute,'gpu')
        y = gpuArray(y);
    end
    
    y(obj.mask_patterns) = x;
%     y = bsxfun(@times,reshape(x,[obj.frame_size(1),obj.frame_size(2),size(obj.C,3),obj.Q]),obj.sampPattern);

    % iFFT
    y = 1/sqrt(size(y,1)*size(y,2))*fft(fft(y,[],1),[],2);
    
    % Coil Multiply
    y = squeeze(sum(y.*(obj.C),3));

    % Reshape to a column vector
    y = y(:);
    
    % put array back on cpu
    if strcmpi(obj.compute,'gpu_off')
        y = gather(y);
    end
end