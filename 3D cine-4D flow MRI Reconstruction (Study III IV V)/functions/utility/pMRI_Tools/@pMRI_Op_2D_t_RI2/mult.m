% Forward 3D Operator
function y = mult(obj,x)

    % Put input on GPU
    if strcmpi(obj.compute,'gpu_off')
        x = gpuArray(x);
    end
    
    % Reshape
    x = x(1:length(x)/2) + 1j*x(length(x)/2+1:end);
    x = reshape(x,obj.frame_size(1),obj.frame_size(2),1,obj.Q);

    % Image domain multiply and FFT
    x = 1/sqrt(size(x,1)*size(x,2))*fft2(bsxfun(@times,x,obj.C));
    
    % Downsample
    y = x(obj.mask_patterns);
%     y = bsxfun(@times,x,obj.sampPattern);

    % put array back on cpu
    if strcmpi(obj.compute,'gpu_off')
        y = gather(y);
    end
end