function [ResCol,ResLin,ResPar,Rate_s,Rate_w,NImageCols,NImageLins,NImagePars] = constructArray_v2(FlowSG4D_Outputs, opt, name, saveLocation, param, raw_header)


SG_Reveal4D_Data = FlowSG4D_Outputs;
clear FlowSG4D_Outputs


ar_size = SG_Reveal4D_Data.ArraySize;
ar_ind = SG_Reveal4D_Data.Indices;
ar_ind = [ar_ind,(1:size(ar_ind,1)).'];
data = SG_Reveal4D_Data.Data;
rWeights = SG_Reveal4D_Data.rWeights;
rWeights = squeeze(single(rWeights));
    
%     ar_size(end) = 1;

% crop = max([round(ar_size(1) / 4), round((ar_size(1)-96)/2)]); % 42 for patient 3d cine
% crop = max([round(ar_size(1) / 4), round((ar_size(1)-384)/2)]);
crop1 = round(opt.crop(1)*ar_size(1));
crop2 = round(opt.crop(2)*ar_size(1));
asym_size = max(find(data(1:end-1,1,1)==0));
asym_percent = asym_size/size(data,1);
% ar_size(1) = ar_size(1)-crop*2;
ar_size(1) = ar_size(1)-(crop1+crop2);
data = ifftshift(ifft(fftshift(data),[],1));
% data = data(crop+1:end-crop,:,:);
data = data(crop1+1:end-crop2,:,:);
data = ifftshift(fft(fftshift(data),[],1));

    %%
    if ~opt.flow
        ar_size(5) = 1;
    end
    for r_ind = 1:ar_size(end)
        k_space = (zeros(ar_size([1,2,3,4,5,6]),'single'));
        weights = (zeros(1,1,ar_size(3),ar_size(4),ar_size(5),ar_size(6),'single'));
        r1_rows = find(ar_ind(:,5)==r_ind);
%         ar_ind_r1 = ar_ind(r1_rows,:);
        ar_ind_r1 = ar_ind;
        
        replace = zeros(1,length(ar_ind_r1));
        counter = 0;
        counter2 = 0;
        % k_space(:,:,ar_ind(ind,1),ar_ind(ind,2),ar_ind(ind,4)) = data(:,:,ind);
        
        rWeights_tmp = rWeights(:,r_ind);
        for ind = 1:length(ar_ind_r1)
            
            if rWeights(ar_ind_r1(ind,6),r_ind) > weights(:,:,ar_ind_r1(ind,1),ar_ind_r1(ind,2),ar_ind_r1(ind,3),ar_ind_r1(ind,4))
                k_space(:,:,ar_ind_r1(ind,1),ar_ind_r1(ind,2),ar_ind_r1(ind,3),ar_ind_r1(ind,4)) = data(:,:,ar_ind_r1(ind,6));
                weights(:,:,ar_ind_r1(ind,1),ar_ind_r1(ind,2),ar_ind_r1(ind,3),ar_ind_r1(ind,4)) = rWeights_tmp(ar_ind_r1(ind,6));
%                 replace(ar_ind_r1(ind,6),r_ind) = 1;
%                 counter = counter +1;
%                 replaceW(counter,r_ind) = rWeights(ar_ind_r1(ind,6),r_ind);
            else
%                 counter2 = counter2 +1;
%                 skip(counter2,r_ind) = rWeights(ar_ind_r1(ind,6),r_ind);
            end
            
            
        end
    
        k_space = permute(k_space,[1,3,4,2,6,5]);
        kb = k_space(:,:,:,:,:,1);
        if opt.flow
            kx = k_space(:,:,:,:,:,2);
            ky = k_space(:,:,:,:,:,3);
            kz = k_space(:,:,:,:,:,4);
        else
            kx = kb;
            ky = kb;
            kz = kb;
        end
        clear k_space
        
        weights = repmat(weights, ar_size(1), ar_size(2), 1, 1, 1, 1);
        weights = permute(weights, [1,3,4,2,6,5]);
        weightsB = weights(:,:,:,:,:,1);
        if opt.flow
            weightsX = weights(:,:,:,:,:,2);
            weightsY = weights(:,:,:,:,:,3);
            weightsZ = weights(:,:,:,:,:,4);
        else
            weightsX = weightsB;
            weightsY = weightsB;
            weightsZ = weightsB;
        end
        clear weights

        sampB = zeros(size(squeeze(kx(:,:,:,1,:))),'logical');
        sampX = zeros(size(squeeze(kx(:,:,:,1,:))),'logical');
        sampY = zeros(size(squeeze(kx(:,:,:,1,:))),'logical');
        sampZ = zeros(size(squeeze(kx(:,:,:,1,:))),'logical');

        sampB(abs(squeeze(kb(:,:,:,1,:)))>0) =1;
        sampX(abs(squeeze(kx(:,:,:,1,:)))>0) =1;
        sampY(abs(squeeze(ky(:,:,:,1,:)))>0) =1;
        sampZ(abs(squeeze(kz(:,:,:,1,:)))>0) =1;
    
        % Enforce Asymmetric Echo
        sampB(1:round(size(sampB,1)*asym_percent),:,:,:,:) = 0;
        sampX(1:round(size(sampX,1)*asym_percent),:,:,:,:) = 0;
        sampY(1:round(size(sampY,1)*asym_percent),:,:,:,:) = 0;
        sampZ(1:round(size(sampZ,1)*asym_percent),:,:,:,:) = 0;
        
        kb = bsxfun(@times,kb,permute(sampB,[1,2,3,5,4]));
        kx = bsxfun(@times,kx,permute(sampX,[1,2,3,5,4]));
        ky = bsxfun(@times,ky,permute(sampY,[1,2,3,5,4]));
        kz = bsxfun(@times,kz,permute(sampZ,[1,2,3,5,4]));
        
        weightsB = bsxfun(@times,weightsB,permute(sampB,[1,2,3,5,4]));
        weightsX = bsxfun(@times,weightsX,permute(sampX,[1,2,3,5,4]));
        weightsY = bsxfun(@times,weightsY,permute(sampY,[1,2,3,5,4]));
        weightsZ = bsxfun(@times,weightsZ,permute(sampZ,[1,2,3,5,4]));
        
        weightsB = squeeze(weightsB(:,:,:,1,:));
        weightsX = squeeze(weightsX(:,:,:,1,:));
        weightsY = squeeze(weightsY(:,:,:,1,:));
        weightsZ = squeeze(weightsZ(:,:,:,1,:));
        
        % ==========================================================================================
        % MAKE RESOLUTION ISOTROPIC (OR AS CLOSE AS POSSIBLE)
        % (This needs more work)
        % ==========================================================================================
        % Maximum allowable matrix size as contraint
        % ******************************************
        maxCols = 96;
        maxLins = 96;
        maxPars = 72;
        
        NImageCols = maxCols;
        NImageLins = maxLins;
        NImagePars = maxPars;
        
%         % ******************************************
%         % Current resolution
        ResCol = SG_Reveal4D_Data.param.ResCol;
        ResLin = SG_Reveal4D_Data.param.ResLin;
        ResPar = SG_Reveal4D_Data.param.ResPar;
        
        
        
%         target_resolution = 1.25;
%         size_new = zeros(1,5);
%         size_new(1) = round((ResCol/target_resolution) * size(kb, 1));
%         size_new(2) = round((ResLin/target_resolution) * size(kb, 2));
%         size_new(3) = round((ResPar/target_resolution) * size(kb, 3));
%         %         size_new_1 = 77;
% %         size_new_2 = 352;
% %         size_new_3 = 75;
% 
%         pad_post = [mod(size_new(1),2),mod(size_new(2),2),mod(size_new(3),2),0,0];      
%         kb = padarray(kb,pad_post,0,'post');
%         kx = padarray(kx,pad_post,0,'post');
%         ky = padarray(ky,pad_post,0,'post');
%         kz = padarray(kz,pad_post,0,'post');
%         pad_post = squeeze(pad_post(:,:,:,1,:));
%         sampB = padarray(sampB,pad_post,0,'post');
%         sampX = padarray(sampX,pad_post,0,'post');
%         sampY = padarray(sampY,pad_post,0,'post');
%         sampZ = padarray(sampZ,pad_post,0,'post');
%         weightsB = padarray(weightsB,pad_post,0,'post');
%         weightsX = padarray(weightsX,pad_post,0,'post');
%         weightsY = padarray(weightsY,pad_post,0,'post');
%         weightsZ = padarray(weightsZ,pad_post,0,'post');  
%         % pad both ends
%         pad_both = (size_new - size(kb))/2; pad_both(4:end) = 0; crop_both = -min(pad_both, 0);
%         pad_both = max(pad_both, 0);
%         pad_both = max(pad_both, 0);
%         kb = padarray(kb,pad_both,0,'both');
%         kx = padarray(kx,pad_both,0,'both');
%         ky = padarray(ky,pad_both,0,'both');
%         kz = padarray(kz,pad_both,0,'both');
%         pad_both = squeeze(pad_both(:,:,:,1,:));
%         sampB = padarray(sampB,pad_both,0,'both');
%         sampX = padarray(sampX,pad_both,0,'both');
%         sampY = padarray(sampY,pad_both,0,'both');
%         sampZ = padarray(sampZ,pad_both,0,'both');
%         weightsB = padarray(weightsB,pad_both,0,'both');
%         weightsX = padarray(weightsX,pad_both,0,'both');
%         weightsY = padarray(weightsY,pad_both,0,'both');
%         weightsZ = padarray(weightsZ,pad_both,0,'both');
%         % Crop both ends
%         cc = crop_both;
%         kb = kb(cc(1)+1:end-cc(1), cc(2)+1:end-cc(2),cc(3)+1:end-cc(3),cc(4)+1:end-cc(4),cc(5)+1:end-cc(5));
%         kx = kx(cc(1)+1:end-cc(1), cc(2)+1:end-cc(2),cc(3)+1:end-cc(3),cc(4)+1:end-cc(4),cc(5)+1:end-cc(5));
%         ky = ky(cc(1)+1:end-cc(1), cc(2)+1:end-cc(2),cc(3)+1:end-cc(3),cc(4)+1:end-cc(4),cc(5)+1:end-cc(5));
%         kz = kz(cc(1)+1:end-cc(1), cc(2)+1:end-cc(2),cc(3)+1:end-cc(3),cc(4)+1:end-cc(4),cc(5)+1:end-cc(5));
%         sampB = sampB(cc(1)+1:end-cc(1), cc(2)+1:end-cc(2),cc(3)+1:end-cc(3),cc(5)+1:end-cc(5));
%         sampX = sampX(cc(1)+1:end-cc(1), cc(2)+1:end-cc(2),cc(3)+1:end-cc(3),cc(5)+1:end-cc(5));
%         sampY = sampY(cc(1)+1:end-cc(1), cc(2)+1:end-cc(2),cc(3)+1:end-cc(3),cc(5)+1:end-cc(5));
%         sampZ = sampZ(cc(1)+1:end-cc(1), cc(2)+1:end-cc(2),cc(3)+1:end-cc(3),cc(5)+1:end-cc(5));
%         weightsB = weightsB(cc(1)+1:end-cc(1), cc(2)+1:end-cc(2),cc(3)+1:end-cc(3),cc(5)+1:end-cc(5));
%         weightsX = weightsX(cc(1)+1:end-cc(1), cc(2)+1:end-cc(2),cc(3)+1:end-cc(3),cc(5)+1:end-cc(5));
%         weightsY = weightsY(cc(1)+1:end-cc(1), cc(2)+1:end-cc(2),cc(3)+1:end-cc(3),cc(5)+1:end-cc(5));
%         weightsZ = weightsZ(cc(1)+1:end-cc(1), cc(2)+1:end-cc(2),cc(3)+1:end-cc(3),cc(5)+1:end-cc(5));
%         ResCol = target_resolution;
%         ResLin = target_resolution;
%         ResPar = target_resolution;     
        
%         
%         % *** Temporary ****
%         padLength = floor((-size(kb,1) + size_new_1)/2);
%         kb = padarray(kb,[padLength,0,0,0,0],0,'both');
%         kx = padarray(kx,[padLength,0,0,0,0],0,'both');
%         ky = padarray(ky,[padLength,0,0,0,0],0,'both');
%         kz = padarray(kz,[padLength,0,0,0,0],0,'both');
%         sampB = padarray(sampB,[padLength,0,0,0],0,'both');
%         sampX = padarray(sampX,[padLength,0,0,0],0,'both');
%         sampY = padarray(sampY,[padLength,0,0,0],0,'both');
%         sampZ = padarray(sampZ,[padLength,0,0,0],0,'both');
%         weightsB = padarray(weightsB,[padLength,0,0,0],0,'both');
%         weightsX = padarray(weightsX,[padLength,0,0,0],0,'both');
%         weightsY = padarray(weightsY,[padLength,0,0,0],0,'both');
%         weightsZ = padarray(weightsZ,[padLength,0,0,0],0,'both');
%         %
%         padLength = floor((-size(kb,2) + size_new_2)/2);
%         kb = padarray(kb,[0,padLength,0,0,0],0,'both');
%         kx = padarray(kx,[0,padLength,0,0,0],0,'both');
%         ky = padarray(ky,[0,padLength,0,0,0],0,'both');
%         kz = padarray(kz,[0,padLength,0,0,0],0,'both');
%         sampB = padarray(sampB,[0,padLength,0,0],0,'both');
%         sampX = padarray(sampX,[0,padLength,0,0],0,'both');
%         sampY = padarray(sampY,[0,padLength,0,0],0,'both');
%         sampZ = padarray(sampZ,[0,padLength,0,0],0,'both');
%         weightsB = padarray(weightsB,[0,padLength,0,0],0,'both');
%         weightsX = padarray(weightsX,[0,padLength,0,0],0,'both');
%         weightsY = padarray(weightsY,[0,padLength,0,0],0,'both');
%         weightsZ = padarray(weightsZ,[0,padLength,0,0],0,'both');
%         %
%         padLength = floor((-size(kb,3) + size_new_3)/2);
%         kb = padarray(kb,[0,0,padLength,0,0],0,'both');
%         kx = padarray(kx,[0,0,padLength,0,0],0,'both');
%         ky = padarray(ky,[0,0,padLength,0,0],0,'both');
%         kz = padarray(kz,[0,0,padLength,0,0],0,'both');
%         sampB = padarray(sampB,[0,0,padLength,0],0,'both');
%         sampX = padarray(sampX,[0,0,padLength,0],0,'both');
%         sampY = padarray(sampY,[0,0,padLength,0],0,'both');
%         sampZ = padarray(sampZ,[0,0,padLength,0],0,'both');
%         weightsB = padarray(weightsB,[0,0,padLength,0],0,'both');
%         weightsX = padarray(weightsX,[0,0,padLength,0],0,'both');
%         weightsY = padarray(weightsY,[0,0,padLength,0],0,'both');
%         weightsZ = padarray(weightsZ,[0,0,padLength,0],0,'both');
%         

        
        
        
        
        
        
        
        
        
        
        
        
%         
%         % Crop kspace to target resolution
%         tres = 1.5;
%         mx_col = size(kb, 1);
%         mx_lin = size(kb, 2);
%         mx_par = size(kb,3);
%         %
%         mx_col_new = round((0.5*(ResCol*mx_col)/tres))*2;
%         mx_lin_new = round((0.5*(ResLin*mx_lin)/tres))*2;
%         mx_par_new = round((0.5*(ResPar*mx_par)/tres))*2;
%         %
%         ResCol = (mx_col*ResCol)/mx_col_new;
%         ResLin = (mx_lin*ResLin)/mx_lin_new;
%         ResPar = (mx_par*ResPar)/mx_par_new;
%         %
%         crop_col = (mx_col - mx_col_new)/2;
%         crop_lin = (mx_lin - mx_lin_new)/2;
%         crop_par = (mx_par - mx_par_new)/2;
%         %
%         kb = kb(crop_col+1:end-crop_col,crop_lin+1:end-crop_lin,crop_par+1:end-crop_par,:,:);
%         kx = kx(crop_col+1:end-crop_col,crop_lin+1:end-crop_lin,crop_par+1:end-crop_par,:,:);
%         ky = ky(crop_col+1:end-crop_col,crop_lin+1:end-crop_lin,crop_par+1:end-crop_par,:,:);
%         kz = kz(crop_col+1:end-crop_col,crop_lin+1:end-crop_lin,crop_par+1:end-crop_par,:,:);
% 
%         sampB = sampB(crop_col+1:end-crop_col,crop_lin+1:end-crop_lin,crop_par+1:end-crop_par,:);
%         sampX = sampX(crop_col+1:end-crop_col,crop_lin+1:end-crop_lin,crop_par+1:end-crop_par,:);
%         sampY = sampY(crop_col+1:end-crop_col,crop_lin+1:end-crop_lin,crop_par+1:end-crop_par,:);
%         sampZ = sampZ(crop_col+1:end-crop_col,crop_lin+1:end-crop_lin,crop_par+1:end-crop_par,:);
% 
%         weightsB = weightsB(crop_col+1:end-crop_col,crop_lin+1:end-crop_lin,crop_par+1:end-crop_par,:);
%         weightsX = weightsX(crop_col+1:end-crop_col,crop_lin+1:end-crop_lin,crop_par+1:end-crop_par,:);
%         weightsY = weightsY(crop_col+1:end-crop_col,crop_lin+1:end-crop_lin,crop_par+1:end-crop_par,:);
%         weightsZ = weightsZ(crop_col+1:end-crop_col,crop_lin+1:end-crop_lin,crop_par+1:end-crop_par,:);
%         
%         SG_Reveal4D_Data.param.ResCol = ResCol;
%         SG_Reveal4D_Data.param.ResLin = ResLin;
%         SG_Reveal4D_Data.param.ResPar = ResPar;
%         SG_Reveal4D_Data.param.NImageCols = mx_col_new;
%         SG_Reveal4D_Data.param.NImageLins = mx_lin_new;
%         SG_Reveal4D_Data.param.NImagePars = mx_par_new;


        
%         % First make sure Col (FE) doesn't exceed maximum
%         NImageCols = min([SG_Reveal4D_Data.param.NImageCols, maxCols]);
%         padLength = (NImageCols - size(kb,1)) / 2;
%         if padLength >= 0
%             % pad
%             kb = padarray(kb,[padLength,0,0,0,0],0,'both');
%             kx = padarray(kx,[padLength,0,0,0,0],0,'both');
%             ky = padarray(ky,[padLength,0,0,0,0],0,'both');
%             kz = padarray(kz,[padLength,0,0,0,0],0,'both');
%             sampB = padarray(sampB,[padLength,0,0,0],0,'both');
%             sampX = padarray(sampX,[padLength,0,0,0],0,'both');
%             sampY = padarray(sampY,[padLength,0,0,0],0,'both');
%             sampZ = padarray(sampZ,[padLength,0,0,0],0,'both');
%             weightsB = padarray(weightsB,[padLength,0,0,0],0,'both');
%             weightsX = padarray(weightsX,[padLength,0,0,0],0,'both');
%             weightsY = padarray(weightsY,[padLength,0,0,0],0,'both');
%             weightsZ = padarray(weightsZ,[padLength,0,0,0],0,'both');
%         else
%             % crop
%             padLength = -padLength;
%             kb = kb(padLength+1:end-padLength,:,:,:,:);
%             kx = kx(padLength+1:end-padLength,:,:,:,:);
%             ky = ky(padLength+1:end-padLength,:,:,:,:);
%             kz = kz(padLength+1:end-padLength,:,:,:,:);
%             sampB = sampB(padLength+1:end-padLength,:,:,:);
%             sampX = sampX(padLength+1:end-padLength,:,:,:);
%             sampY = sampY(padLength+1:end-padLength,:,:,:);
%             sampZ = sampZ(padLength+1:end-padLength,:,:,:);
%             weightsB = weightsB(padLength+1:end-padLength,:,:,:);
%             weightsX = weightsX(padLength+1:end-padLength,:,:,:);
%             weightsY = weightsY(padLength+1:end-padLength,:,:,:);
%             weightsZ = weightsZ(padLength+1:end-padLength,:,:,:);
%         end
%         ResCol = ResCol/((size(kb,1)+2*padLength)/size(kb,1));
%         
%         
%         % Next, make in-plane resolution isotropic, Lin (PE) if possible
%         NImageLins = min([NImageCols, maxLins]);
%         padLength = (NImageLins - size(kb,2)) / 2;
%         if padLength >= 0
%             % pad
%             kb = padarray(kb,[0,padLength,0,0,0],0,'both');
%             kx = padarray(kx,[0,padLength,0,0,0],0,'both');
%             ky = padarray(ky,[0,padLength,0,0,0],0,'both');
%             kz = padarray(kz,[0,padLength,0,0,0],0,'both');
%             sampB = padarray(sampB,[0,padLength,0,0],0,'both');
%             sampX = padarray(sampX,[0,padLength,0,0],0,'both');
%             sampY = padarray(sampY,[0,padLength,0,0],0,'both');
%             sampZ = padarray(sampZ,[0,padLength,0,0],0,'both');
%             weightsB = padarray(weightsB,[0,padLength,0,0],0,'both');
%             weightsX = padarray(weightsX,[0,padLength,0,0],0,'both');
%             weightsY = padarray(weightsY,[0,padLength,0,0],0,'both');
%             weightsZ = padarray(weightsZ,[0,padLength,0,0],0,'both');
%         else
%             % crop
%             padLength = - padLength;
%             kb = kb(:,padLength+1:end-padLength,:,:,:);
%             kx = kx(:,padLength+1:end-padLength,:,:,:);
%             ky = ky(:,padLength+1:end-padLength,:,:,:);
%             kz = kz(:,padLength+1:end-padLength,:,:,:);
%             sampB = sampB(:,padLength+1:end-padLength,:,:);
%             sampX = sampX(:,padLength+1:end-padLength,:,:);
%             sampY = sampY(:,padLength+1:end-padLength,:,:);
%             sampZ = sampZ(:,padLength+1:end-padLength,:,:);
%             weightsB = weightsB(:,padLength+1:end-padLength,:,:);
%             weightsX = weightsX(:,padLength+1:end-padLength,:,:);
%             weightsY = weightsY(:,padLength+1:end-padLength,:,:);
%             weightsZ = weightsZ(:,padLength+1:end-padLength,:,:);
%         end
%         ResLin = ResLin/((size(kb,2)+2*padLength)/size(kb,2));     
%         
% %         % Finally, make Par (partitions/slices) isotropic if possible
% %         ResPar_tmp = min([ResCol, ResLin]);
% %         if (size(kb,3)*ResPar)/ResPar_tmp <= maxPars
% %             NImagePars = (size(kb,3)*ResPar)/ResPar_tmp;
% %             padLength = floor((NImagePars - size(kb,3)) / 2);
% %             if padLength >= 0
% %                 % pad
% %                 kb = padarray(kb,[0,0,padLength,0,0],0,'both');
% %                 kx = padarray(kx,[0,0,padLength,0,0],0,'both');
% %                 ky = padarray(ky,[0,0,padLength,0,0],0,'both');
% %                 kz = padarray(kz,[0,0,padLength,0,0],0,'both');
% %                 sampB = padarray(sampB,[0,0,padLength,0],0,'both');
% %                 sampX = padarray(sampX,[0,0,padLength,0],0,'both');
% %                 sampY = padarray(sampY,[0,0,padLength,0],0,'both');
% %                 sampZ = padarray(sampZ,[0,0,padLength,0],0,'both');
% %                 weightsB = padarray(weightsB,[0,0,padLength,0],0,'both');
% %                 weightsX = padarray(weightsX,[0,0,padLength,0],0,'both');
% %                 weightsY = padarray(weightsY,[0,0,padLength,0],0,'both');
% %                 weightsZ = padarray(weightsZ,[0,0,padLength,0],0,'both');
% %             else
% %                 % crop
% %                 kb = kb(:,:,padLength+1:end-padLength,:,:);
% %                 kx = kx(:,:,padLength+1:end-padLength,:,:);
% %                 ky = ky(:,:,padLength+1:end-padLength,:,:);
% %                 kz = kz(:,:,padLength+1:end-padLength,:,:);
% %                 sampB = sampB(:,:,padLength+1:end-padLength,:);
% %                 sampX = sampX(:,:,padLength+1:end-padLength,:);
% %                 sampY = sampY(:,:,padLength+1:end-padLength,:);
% %                 sampZ = sampZ(:,:,padLength+1:end-padLength,:);
% %                 weightsB = weightsB(:,:,padLength+1:end-padLength,:);
% %                 weightsX = weightsX(:,:,padLength+1:end-padLength,:);
% %                 weightsY = weightsY(:,:,padLength+1:end-padLength,:);
% %                 weightsZ = weightsZ(:,:,padLength+1:end-padLength,:);
% %             end
% %             ResPar = ResPar/((size(kb,3)+2*padLength)/size(kb,3));
% %             
% %         else
% %             ResPar_tmp = max([ResCol, ResLin]);
% %             if (size(kb,3)*ResPar)/ResPar_tmp <= maxPars
% %                 NImagePars = (size(kb,3)*ResPar)/ResPar_tmp;
% %                 padLength = floor((NImagePars - size(kb,3)) / 2);
% %                 if padLength >= 0
% %                     % pad
% %                     kb = padarray(kb,[0,0,padLength,0,0],0,'both');
% %                     kx = padarray(kx,[0,0,padLength,0,0],0,'both');
% %                     ky = padarray(ky,[0,0,padLength,0,0],0,'both');
% %                     kz = padarray(kz,[0,0,padLength,0,0],0,'both');
% %                     sampB = padarray(sampB,[0,0,padLength,0],0,'both');
% %                     sampX = padarray(sampX,[0,0,padLength,0],0,'both');
% %                     sampY = padarray(sampY,[0,0,padLength,0],0,'both');
% %                     sampZ = padarray(sampZ,[0,0,padLength,0],0,'both');
% %                     weightsB = padarray(weightsB,[0,0,padLength,0],0,'both');
% %                     weightsX = padarray(weightsX,[0,0,padLength,0],0,'both');
% %                     weightsY = padarray(weightsY,[0,0,padLength,0],0,'both');
% %                     weightsZ = padarray(weightsZ,[0,0,padLength,0],0,'both');
% %                 else
% %                     % crop
% %                     kb = kb(:,:,padLength+1:end-padLength,:,:);
% %                     kx = kx(:,:,padLength+1:end-padLength,:,:);
% %                     ky = ky(:,:,padLength+1:end-padLength,:,:);
% %                     kz = kz(:,:,padLength+1:end-padLength,:,:);
% %                     sampB = sampB(:,:,padLength+1:end-padLength,:);
% %                     sampX = sampX(:,:,padLength+1:end-padLength,:);
% %                     sampY = sampY(:,:,padLength+1:end-padLength,:);
% %                     sampZ = sampZ(:,:,padLength+1:end-padLength,:);
% %                     weightsB = weightsB(:,:,padLength+1:end-padLength,:);
% %                     weightsX = weightsX(:,:,padLength+1:end-padLength,:);
% %                     weightsY = weightsY(:,:,padLength+1:end-padLength,:);
% %                     weightsZ = weightsZ(:,:,padLength+1:end-padLength,:);
% %                 end
% %                 ResPar = ResPar/((size(kb,3)+2*padLength)/size(kb,3));
% %                 
% %                 
% %             else
%                 NImagePars = maxPars;
%                 padLength = floor((NImagePars - size(kb,3)) / 2);
%                 if padLength >= 0
%                     % pad
%                     kb = padarray(kb,[0,0,padLength,0,0],0,'both');
%                     kx = padarray(kx,[0,0,padLength,0,0],0,'both');
%                     ky = padarray(ky,[0,0,padLength,0,0],0,'both');
%                     kz = padarray(kz,[0,0,padLength,0,0],0,'both');
%                     sampB = padarray(sampB,[0,0,padLength,0],0,'both');
%                     sampX = padarray(sampX,[0,0,padLength,0],0,'both');
%                     sampY = padarray(sampY,[0,0,padLength,0],0,'both');
%                     sampZ = padarray(sampZ,[0,0,padLength,0],0,'both');
%                     weightsB = padarray(weightsB,[0,0,padLength,0],0,'both');
%                     weightsX = padarray(weightsX,[0,0,padLength,0],0,'both');
%                     weightsY = padarray(weightsY,[0,0,padLength,0],0,'both');
%                     weightsZ = padarray(weightsZ,[0,0,padLength,0],0,'both');
%                 else
%                     % crop
%                     kb = kb(:,:,padLength+1:end-padLength,:,:);
%                     kx = kx(:,:,padLength+1:end-padLength,:,:);
%                     ky = ky(:,:,padLength+1:end-padLength,:,:);
%                     kz = kz(:,:,padLength+1:end-padLength,:,:);
%                     sampB = sampB(:,:,padLength+1:end-padLength,:);
%                     sampX = sampX(:,:,padLength+1:end-padLength,:);
%                     sampY = sampY(:,:,padLength+1:end-padLength,:);
%                     sampZ = sampZ(:,:,padLength+1:end-padLength,:);
%                     weightsB = weightsB(:,:,padLength+1:end-padLength,:);
%                     weightsX = weightsX(:,:,padLength+1:end-padLength,:);
%                     weightsY = weightsY(:,:,padLength+1:end-padLength,:);
%                     weightsZ = weightsZ(:,:,padLength+1:end-padLength,:);
%                 end
%                 ResPar = ResPar/((size(kb,3)+2*padLength)/size(kb,3));
%                 
% %             end  
% %         end
%         
%         
%         
% %         % First make in-plane resolution isotropic:
% %         ResCol = SG_Reveal4D_Data.param.ResCol;
% %         ResLin = SG_Reveal4D_Data.param.ResLin;
% %         ResPar = SG_Reveal4D_Data.param.ResPar;
% %         
% %         NImageLins = min(SG_Reveal4D_Data.param.NImageLins, maxLin);
% %         padLength = (NImageLins - size(kb,2)) / 2;
% %         
% %         kb = padarray(kb,[0,padLength,0,0,0],0,'both');
% %         kx = padarray(kx,[0,padLength,0,0,0],0,'both');
% %         ky = padarray(ky,[0,padLength,0,0,0],0,'both');
% %         kz = padarray(kz,[0,padLength,0,0,0],0,'both');
% %         
% %         sampB = padarray(sampB,[0,padLength,0,0],0,'both');
% %         sampX = padarray(sampX,[0,padLength,0,0],0,'both');
% %         sampY = padarray(sampY,[0,padLength,0,0],0,'both');
% %         sampZ = padarray(sampZ,[0,padLength,0,0],0,'both');
% %         
% %         weightsB = padarray(weightsB,[0,padLength,0,0],0,'both');
% %         weightsX = padarray(weightsX,[0,padLength,0,0],0,'both');
% %         weightsY = padarray(weightsY,[0,padLength,0,0],0,'both');
% %         weightsZ = padarray(weightsZ,[0,padLength,0,0],0,'both');
% %         
% %         ResLin = ResCol;
% %         
% %         % Second, check if even padding can be acheived for partitions
% %         
% %         padLength1 = 0;
% %         padLength2 = 0;
% %         padLength3 = 0;
% %         
% %         if ResCol < ResPar
% %             tmp = round((ResPar/ResCol)*size(kb,3)) - size(kb,3);
% %             if mod(tmp,2) == 0
% %                 padLength3 = tmp/2;
% %             else
% %                 padLength3 = (tmp-1)/2;
% %             end
% %             % Update resolution
% %             ResPar = ResPar/((size(kb,3)+2*padLength3)/size(kb,3));
% %             % Pad arrays
% %             kb = padarray(kb,[0,0,padLength3,0,0],0,'both');
% %             kx = padarray(kx,[0,0,padLength3,0,0],0,'both');
% %             ky = padarray(ky,[0,0,padLength3,0,0],0,'both');
% %             kz = padarray(kz,[0,0,padLength3,0,0],0,'both');
% %             sampB = padarray(sampB,[0,0,padLength3,0],0,'both');
% %             sampX = padarray(sampX,[0,0,padLength3,0],0,'both');
% %             sampY = padarray(sampY,[0,0,padLength3,0],0,'both');
% %             sampZ = padarray(sampZ,[0,0,padLength3,0],0,'both');
% %             weightsB = padarray(weightsB,[0,0,padLength3,0],0,'both');
% %             weightsX = padarray(weightsX,[0,0,padLength3,0],0,'both');
% %             weightsY = padarray(weightsY,[0,0,padLength3,0],0,'both');
% %             weightsZ = padarray(weightsZ,[0,0,padLength3,0],0,'both');
% %             
% %         elseif ResCol > ResPar
% %             tmp1 = round((ResCol/ResPar)*size(kb,1)) - size(kb,1);
% %             tmp2 = round((ResLin/ResPar)*size(kb,2)) - size(kb,2);
% %             if mod(tmp1,2) == 0
% %                 padLength1 = tmp1/2;
% %                 if mod(tmp2,2) == 0
% %                     padLength2 = tmp2/2;
% %                 else
% %                     padLength2 = (tmp2-1)/2;
% %                 end
% %             else
% %                 padLength1 = (tmp1-1)/2;
% %                 if mod(tmp2,2) == 0
% %                     padLength2 = tmp2/2;
% %                 else
% %                     padLength2 = (tmp2-1)/2;
% %                 end
% %             end
% %             % Update resolution
% %             ResCol = ResCol/((size(kb,1)+2*padLength1)/size(kb,1));
% %             ResLin = ResLin/((size(kb,2)+2*padLength2)/size(kb,2));
% %             % Pad arrays
% %             kb = padarray(kb,[padLength1,padLength2,0,0,0],0,'both');
% %             kx = padarray(kx,[padLength1,padLength2,0,0,0],0,'both');
% %             ky = padarray(ky,[padLength1,padLength2,0,0,0],0,'both');
% %             kz = padarray(kz,[padLength1,padLength2,0,0,0],0,'both');
% %             sampB = padarray(sampB,[padLength1,padLength2,0,0],0,'both');
% %             sampX = padarray(sampX,[padLength1,padLength2,0,0],0,'both');
% %             sampY = padarray(sampY,[padLength1,padLength2,0,0],0,'both');
% %             sampZ = padarray(sampZ,[padLength1,padLength2,0,0],0,'both');
% %             weightsB = padarray(weightsB,[padLength1,padLength2,0,0],0,'both');
% %             weightsX = padarray(weightsX,[padLength1,padLength2,0,0],0,'both');
% %             weightsY = padarray(weightsY,[padLength1,padLength2,0,0],0,'both');
% %             weightsZ = padarray(weightsZ,[padLength1,padLength2,0,0],0,'both');
% %             
% %         end

        scanParam = SG_Reveal4D_Data.param;
        scanParam.totalbeats=opt.numofbeats;
        scanParam.totaltime= opt.totalscantime;
        scanParam.flow= opt.flow;
        scanParam.NImageCols = size(kb,1);
        scanParam.NImageLins = size(kb,2);
        scanParam.NImagePars = size(kb,3);
        scanParam.ResCol = ResCol;
        scanParam.ResLin = ResLin;
        scanParam.ResPar = ResPar;
        
        
        dirName = [saveLocation,'\',name];
        if ~exist(dirName,'dir')
            mkdir(dirName);
        end
        
                
        % estimate acceleration rate
        % do not include asymmetric echo in rate calculation
        Rate_s = (size(sampB,2)*size(sampB,3)*size(sampB,4)) / sum(sum(sum(sampB(round(size(sampB,1)/2),:,:,:),2),3),4);
        Rate_w = (size(sampB,2)*size(sampB,3)*size(sampB,4)) / sum(sum(sum(weightsB(round(size(sampB,1)/2),:,:,:).^2,2),3),4);
        
        scanParam.AccelerationRate_S = Rate_s;
        scanParam.AccelerationRate_W = Rate_w;

        convert_binning_output(dirName, name, kb, kx, ky, kz, sampB, sampX, sampY, sampZ, weightsB, weightsX, weightsY, weightsZ, scanParam);

%         Rate_s = prod(size(sampB)) / sum(sampB(:));
%         Rate_w = prod(size(sampB)) / sum(weightsB(:).^2);
        disp(['Acceleration rate...']);
        disp(['Samples: ',num2str(Rate_s)]);
        disp(['Weighted: ',num2str(Rate_w)]);        
        
%         % estimate acceleration rate
%         % do not include asymmetric echo in rate calculation
%         Rate_s = (size(sampB,2)*size(sampB,3)*size(sampB,4)) / sum(sum(sum(sampB(round(size(sampB,1)/2),:,:,:),2),3),4);
%         Rate_w = (size(sampB,2)*size(sampB,3)*size(sampB,4)) / sum(sum(sum(weightsB(round(size(sampB,1)/2),:,:,:).^2,2),3),4);
%         
%         scanParam.AccelerationRate_S = Rate_s;
%         scanParam.AccelerationRate_W = Rate_w;
% 
% %         Rate_s = prod(size(sampB)) / sum(sampB(:));
% %         Rate_w = prod(size(sampB)) / sum(weightsB(:).^2);
%         disp(['Acceleration rate...']);
%         disp(['Samples: ',num2str(Rate_s)]);
%         disp(['Weighted: ',num2str(Rate_w)]);
%         
%         save([dirName,'\',name,'_resp',num2str(r_ind),'_dataB'],'kb','-v7.3');
%         save([dirName,'\',name,'_resp',num2str(r_ind),'_dataX'],'kx','-v7.3');
%         save([dirName,'\',name,'_resp',num2str(r_ind),'_dataY'],'ky','-v7.3');
%         save([dirName,'\',name,'_resp',num2str(r_ind),'_dataZ'],'kz','-v7.3');
% 
%         save([dirName,'\',name,'_resp',num2str(r_ind),'_sampB'],'sampB','-v7.3');
%         save([dirName,'\',name,'_resp',num2str(r_ind),'_sampX'],'sampX','-v7.3');
%         save([dirName,'\',name,'_resp',num2str(r_ind),'_sampY'],'sampY','-v7.3');
%         save([dirName,'\',name,'_resp',num2str(r_ind),'_sampZ'],'sampZ','-v7.3');
%         
%         save([dirName,'\',name,'_resp',num2str(r_ind),'_weightsB'],'weightsB','-v7.3');
%         save([dirName,'\',name,'_resp',num2str(r_ind),'_weightsX'],'weightsX','-v7.3');
%         save([dirName,'\',name,'_resp',num2str(r_ind),'_weightsY'],'weightsY','-v7.3');
%         save([dirName,'\',name,'_resp',num2str(r_ind),'_weightsZ'],'weightsZ','-v7.3');
%         
%         save([dirName,'\',name,'_resp',num2str(r_ind),'_scanParam'],'scanParam','-v7.3');
%         save([dirName,'\',name,'_resp',num2str(r_ind),'_raw_header'],'raw_header','-v7.3');

    end
    
    clear kx ky kz
    
%     k_space = padarray(k_space,[padLength1,padLength+padLength2,padLength3,0,0,0],0,'both');
    k_space_n = kb;
    k_space_n(abs(k_space_n)==0) = NaN;
    k_space_n = mean(k_space_n,5,'omitnan');
%     k_space_n = mean(k_space_n,6,'omitnan');
    k_space_n(isnan(k_space_n))=0;
    k_space_n = squeeze(k_space_n);
    im = ifft3_shift(k_space_n);
    im_c = im;
    im = sos_combine(im);

    %%
    figure;
    mn = min(abs(im(:)));
    mx = max(abs(im(:)))/5;
    imagesc(abs(im(:,:,round(size(im,3)/2))),[mn,mx]); axis('off','image'); colormap('gray')
    figure;
    imagesc(squeeze(abs(im(round(size(im,1)/2),:,:))),[mn,mx]); axis('off','image'); colormap('gray');
    figure;
    imagesc(squeeze(abs(im(:,round(size(im,2)/2),:))),[mn,mx]); axis('off','image'); colormap('gray');
    
    figure;
    nC = size(im_c);
    mn = min(abs(im_c(:)))/2; mx = max(abs(im_c(:)))/3;
    dim1 = 6; dim2 = 6;
    for i = 1 : size(im_c, 4)
        subplot(dim1, dim2, i)
        imagesc(abs(im_c(:,:,round(size(im,3)/2),i)),[mn, mx]); axis('off','image'); colormap('gray');
    end
        
        
        
        
        
        
    
end



