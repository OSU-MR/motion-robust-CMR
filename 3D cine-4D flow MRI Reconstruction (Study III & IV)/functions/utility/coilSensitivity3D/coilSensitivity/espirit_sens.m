function smaps = espirit_sens(kdata,idx,kN,thresh)


%Input:
%   kdata = k-space data [Nx Ny Ncoils Nframes]
%   idx   = sampling pattern [Nx Ny Nframes]
% Output
%   sos combined image
% Notes
%   This is just meant to be an example recon. It is based on:
%      An Eigen-Vector Approach to AutoCalibrating Parallel MRI, Where
%      SENSE Meets GRAPPA. Lustig et al. ISMRM 11. pg 479
 




[Nx, Ny, Ncoils, Nframes] = size(kdata);

%%%Get Coordinates of Kernel
disp('   ESPIRIT::Collect Kernel');
rx =  [-kN:kN];
ry =  [-kN:kN];

f = ones(2*kN+1,2*kN+1);
f = f/numel(f);
mask = convn(idx,f,'same')> 0.9999;


% Image
A = zeros(sum(mask(:)),(2*kN+1)^2*Ncoils);
count=0;
for frame=1:Nframes
    for posx = 1:size(mask,1)
        for posy = 1:size(mask,2)
            if mask(posx,posy,frame)
                source_vals = kdata(rx+posx,ry+posy,:,frame);
                count = count+1;
                A(count,:) = conj(source_vals(:)');
            end
        end
    end
end

%% SVD
disp('   ESPIRIT::SVD');
[u,s,v] = svd(A,0);
s= diag(s);
sN =s/max(abs(s));
idxS = sN > thresh;
vpp = v(:,idxS);
nV = sum(idxS);

%% Put Kernel in Image Domain
disp('   ESPIRIT::FFT Kernels');
G = zeros([numel(rx) numel(ry) Ncoils nV]);
M = zeros([numel(rx) numel(ry) Ncoils]);
for v=1:nV
    M(:)=vpp(:,v);
    G(:,:,:,v)=M;
end

Cx = floor(Nx/2)+1;
Cy = floor(Ny/2)+1;
for c1 =1:Ncoils
    for c2 =1:nV
        kTemp= zeros(Nx,Ny);
        kTemp(Cx+rx,Cy+ry)=G(:,:,c1,c2);
        Gimage(:,:,c1,c2)=chop_fft2( kTemp );
    end
end

%% How Do Eigen Value Decom
disp('   ESPIRIT::Eigen Decomp (slow)')
S = zeros([Nx Ny]);
U = zeros([Nx Ny Ncoils Ncoils]);
for x=1:Nx
    for y=1:Ny
        A = squeeze(Gimage(x,y,:,:));
        [u,s,v] = svd(A,0);
        V(x,y,:,:)=s(1,1);
        U(x,y,:,:)=conj(u);
    end
end

%%% Smaps  TV(phase) + (Mag-I)
smaps=  ( squeeze(U(:,:,:,1)));
smaps = smaps./ repmat(  sqrt(sum(abs(smaps).^2,3)),[1 1 Ncoils]);
smaps = flipdim(smaps,1);
smaps = flipdim(smaps,2);


%% Kdata to Image
% function im = forward( kdata , smaps)
% 
% [Nx, Ny, Ncoils, Nframes] = size(kdata);
% im=single(zeros(Nx,Ny,Nframes));
% for coil=1:Ncoils
%     for frame = 1:Nframes
%         im(:,:,frame) = im(:,:,frame) + conj(smaps(:,:,coil)).*chop_fft2( kdata(:,:,coil,frame) );
%     end
% end

function im = chop_fft2( im)
[Nx, Ny] = size(im);
% persistent chop;
% if isempty(chop)
    [x,y] = meshgrid(1:Ny,1:Nx);
    chop = (-1).^(x+y);
% end
im = chop.*ifft2(chop.*im) * sqrt(Nx*Ny);



%% Image to Kdata
% function kdata = backward( im, smaps)
% 
% [Nx Ny Nframes] = size(im);
% Ncoils = size(smaps,3);
% 
% kdata =single(zeros(Nx,Ny,Ncoils,Nframes));
% for frame = 1:Nframes
%     for coil=1:Ncoils
%         kdata(:,:,coil,frame) =  chop_ifft2( smaps(:,:,coil).*im(:,:,frame) );
%     end
% end


% function im = chop_ifft2( im)
% [Nx, Ny] = size(im);
% persistent chop;
% if isempty(chop)
%     [x,y] = meshgrid(1:Ny,1:Nx);
%     chop = (-1).^(x+y);
% end
% im = chop.*fft2(chop.*im) * 1/sqrt(Nx*Ny);















% %Input:
% %   kdata = k-space data [Nx Ny Ncoils Nframes]
% %   idx   = sampling pattern [Nx Ny Nframes]
% % Output
% %   sos combined image
% % Notes
% %   This is just meant to be an example recon. It is based on:
% %      An Eigen-Vector Approach to AutoCalibrating Parallel MRI, Where
% %      SENSE Meets GRAPPA. Lustig et al. ISMRM 11. pg 479
%  
% 
% 
% [Nx Ny Ncoils Nframes]=size(kdata);
% 
% 
% 
% %Get sensitivity map
% disp('Get coil sensitivities');
% smaps = espirit(kdata,idx,3,0.01);
% 
% 
% 
% %% Espirit
% function smaps = espirit(kdata,idx,kN,thresh)
% 
% [Nx Ny Ncoils Nframes] = size(kdata);
% 
% %%%Get Coordinates of Kernel
% disp('   ESPIRIT::Collect Kernel');
% rx =  [-kN:kN];
% ry =  [-kN:kN];
% 
% f = ones(2*kN+1,2*kN+1);
% f = f/numel(f);
% mask = convn(idx,f,'same')> 0.9999;
% 
% 
% % Image
% A = zeros(sum(mask(:)),(2*kN+1)^2*Ncoils);
% count=0;
% for frame=1:Nframes
%     for posx = 1:size(mask,1)
%         for posy = 1:size(mask,2)
%             if mask(posx,posy,frame)
%                 source_vals = kdata(rx+posx,ry+posy,:,frame);
%                 count = count+1;
%                 A(count,:) = conj(source_vals(:)');
%             end
%         end
%     end
% end
% 
% %% SVD
% disp('   ESPIRIT::SVD');
% [u,s,v] = svd(A,0);
% s= diag(s);
% sN =s/max(abs(s));
% idxS = sN > thresh;
% vpp = v(:,idxS);
% nV = sum(idxS);
% 
% %% Put Kernel in Image Domain
% disp('   ESPIRIT::FFT Kernels');
% G = zeros([numel(rx) numel(ry) Ncoils nV]);
% M = zeros([numel(rx) numel(ry) Ncoils]);
% for v=1:nV
%     M(:)=vpp(:,v);
%     G(:,:,:,v)=M;
% end
% 
% Cx = floor(Nx/2)+1;
% Cy = floor(Ny/2)+1;
% for c1 =1:Ncoils
%     for c2 =1:nV
%         kTemp= zeros(Nx,Ny);
%         kTemp(Cx+rx,Cy+ry)=G(:,:,c1,c2);
%         Gimage(:,:,c1,c2)=chop_fft2( kTemp );
%     end
% end
% 
% %% How Do Eigen Value Decom
% disp('   ESPIRIT::Eigen Decomp (slow)')
% S = zeros([Nx Ny]);
% U = zeros([Nx Ny Ncoils Ncoils]);
% for x=1:Nx
%     for y=1:Ny
%         A = squeeze(Gimage(x,y,:,:));
%         [u,s,v] = svd(A,0);
%         V(x,y,:,:)=s(1,1);
%         U(x,y,:,:)=conj(u);
%     end
% end
% 
% %%% Smaps  TV(phase) + (Mag-I)
% smaps=  ( squeeze(U(:,:,:,1)));
% smaps = smaps./ repmat(  sqrt(sum(abs(smaps).^2,3)),[1 1 Ncoils]);
% smaps = flipdim(smaps,1);
% smaps = flipdim(smaps,2);
% 
% 
% % %%Coil Compression
% % function kdata = coil_compression(kdata,channels)
% % [Nx Ny Ncoils Nframes] = size(kdata);
% % 
% % A = permute(kdata,[1 2 4 3]);
% % A = reshape(A,[Nx*Ny*Nframes Ncoils]);
% % [u,s,v] = svd(A,0);
% % vpp = v(:,1:channels);
% % A = A*vpp;
% % kdata = reshape(A,[Nx Ny Nframes channels]);
% % kdata = permute(kdata,[1 2 4 3]);
% % 
% % 
% % %% Kdata to Image
% % function im = forward( kdata , smaps)
% % 
% % [Nx Ny Ncoils Nframes] = size(kdata);
% % im=single(zeros(Nx,Ny,Nframes));
% % for coil=1:Ncoils
% %     for frame = 1:Nframes
% %         im(:,:,frame) = im(:,:,frame) + conj(smaps(:,:,coil)).*chop_fft2( kdata(:,:,coil,frame) );
% %     end
% % end
% % 
% function im = chop_fft2( im)
% [Nx Ny] = size(im);
% persistent chop;
% if isempty(chop)
%     [x,y] = meshgrid(1:Ny,1:Nx);
%     chop = (-1).^(x+y);
% end
% im = chop.*fft2(chop.*im);
% 
% % 
% % 
% % %% Image to Kdata
% % function kdata = backward( im, smaps)
% % 
% % [Nx Ny Nframes] = size(im);
% % Ncoils = size(smaps,3);
% % 
% % kdata =single(zeros(Nx,Ny,Ncoils,Nframes));
% % for frame = 1:Nframes
% %     for coil=1:Ncoils
% %         kdata(:,:,coil,frame) =  chop_ifft2( smaps(:,:,coil).*im(:,:,frame) );
% %     end
% % end
% % 
% % 
% function im = chop_ifft2( im)
% [Nx Ny] = size(im);
% persistent chop;
% if isempty(chop)
%     [x,y] = meshgrid(1:Ny,1:Nx);
%     chop = (-1).^(x+y);
% end
% im = chop.*ifft2(chop.*im);
% 
% 
% 




