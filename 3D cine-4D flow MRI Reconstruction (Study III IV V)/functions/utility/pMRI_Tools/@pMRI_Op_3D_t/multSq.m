function y = multSq(obj,x)

% Use Rank 1 approximation for |A|.^2
if obj.uniform_var == 1
    y = ones(obj.M,1)*(obj.Fro2*sum(x,1));
    
% Use true |A|.^2
else
    
    % Put input on GPU
    if strcmpi(obj.compute,'gpu_off')
        x = gpuArray(x);
    end
    
    % Perform True Sq multiply
    x = reshape(x,[obj.frame_size(1)*obj.frame_size(2)*obj.frame_size(3),obj.Q]);
    y(1,1,1,:,:) = obj.CSq*x;
    y = repmat(y,[obj.frame_size(1),obj.frame_size(2),obj.frame_size(3),1,1]);

    % Downsample
    y = y(obj.mask_patterns);
 
    % Weighted least squares
    if ~isempty(obj.mask_weights)
        y = (obj.mask_weights.^2).*y;
    end
    
    % put array back on cpu
    if strcmpi(obj.compute,'gpu_off')
        y = gather(y);
    end
end

