% Forward 3D Operator
function y = mult(obj,x)

    % Put input on GPU
    if strcmpi(obj.compute,'gpu_off')
        x = gpuArray(x);
    end
    
    % Number of sets of sensitivity maps
    nMaps = size(obj.C, 5);
    
    % Reshape
    x = reshape(x,[obj.frame_size(1),obj.frame_size(2),obj.frame_size(3),1,obj.Q,nMaps]);

    % Image domain multiply and FFT
    x = repmat(x,[1,1,1,size(obj.C,4),1,1]);
    
    % If multiple espirit sensitivity maps are used
%     if nMaps > 1
%         arg = bxsfun(@times, 
%         arg = cat(5,bsxfun(@times,x(:,:,:,:,1:obj.Q/2),obj.C(:,:,:,:,1)),bsxfun(@times,x(:,:,:,:,(obj.Q/2)+1:obj.Q),obj.C(:,:,:,:,2)));
%     else
%         arg = bsxfun(@times,x,obj.C);
%     end
%     
%    arg = sum(bsxfun(@times, x, permute(obj.C, [1,2,3,4,6,5])),6);
% % %     
%    x = 1/sqrt(size(x,1)*size(x,2)*size(x,3))*fft(fft2(arg),[],3);
% 
 x = sum(bsxfun(@times, x, permute(obj.C, [1,2,3,4,6,5])),6);
% x=fft(x,[],1);
% x=fft(x,[],2);
% x=fft(x,[],3);
% x=x./sqrt(size(x,1));
 x = 1/sqrt(size(x,1)*size(x,2)*size(x,3))*fft(fft(fft(x,[],1),[],2),[],3); 

    % Downsample
    y = x(obj.mask_patterns);
    
    % Weighted leas squares
    if ~isempty(obj.mask_weights)
        y = obj.mask_weights.*y;
    end
    
    % put array back on cpu
    if strcmpi(obj.compute,'gpu_off')
        y = gather(y);
    end
end