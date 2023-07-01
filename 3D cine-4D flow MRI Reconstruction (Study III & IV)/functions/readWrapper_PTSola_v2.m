function [data, samp, param, cMaps, rawData, pilot_tone_raw] = readWrapper_PTSola_v2(p)
% ===============================
% Rizwan Ahmad (ahmad.46@osu.edu)
% A wrapper for reading Siemens data
% Last Modified: 12/12/2014
% ===============================



%% Read user defined parameters
if isfield(p,'Nc'),      Nc       = p.Nc;         else, Nc       = [100, 0.5];  end % Number of compressed coils 
if isfield(p,'fd'),      fd       = p.fd;         else, fd       = [0.0,0.0];   end % fraction of FE discarded 
if isfield(p,'nAmp'),    nAmp     = p.nAmp;       else, nAmp     = 1;           end % noise amplification factor
if isfield(p,'nStd0'),   nStd     = p.nStd0;      else, nStd     = 1e-5;        end % target noise std
if isfield(p,'nNrm'),    nNrm     = p.nNrm;       else, nNrm     = 0;           end % use the same noise std for each coil? 1:yes, 0:no
if isfield(p,'dis'),     dis      = p.dis;        else, dis      = 1;           end % display time-averaged images?
if isfield(p,'samp'),    samp     = p.samp;       else, samp     = [];          end % sampling pattern
if isfield(p,'sDim'),    sDim     = p.sDim;       else, sDim     = [];          end % dimensions along with subsampling takes place
if isfield(p,'sgAvg'),   sgAvg    = p.sgAvg;      else, sgAvg    = 1;           end % to average the "seg" dimension or not
if isfield(p,'mxFlg'),   mxFlg    = p.mxFlg;      else, mxFlg    = 0;           end % to find Maxwell correctio maps or not
if isfield(p,'kNoise'),  kNoise   = p.kNoise;     else, kNoise   = 0;           end % ignore prescan and use k-space for noise power estimation
if isfield(p,'tap'),     tap      = p.tap;        else, tap      = 0;           end % taper the boundaries along read out   
if isfield(p,'peSft'),   peSft    = p.peSft;      else, peSft    = 0;           end % circular shifting in pe direction   
if isfield(p,'sgI'),     sgI      = p.sgI;        else, sgI      = 0;           end % Self-gating interval
if isfield(p,'flow'),    flow     = p.flow;       else, flow     = 1;           end % flow dataset or not


if isfield(p,'fName'),   fName    = p.fName;      else, error('file name missing'); end  % file name

%% VB or VD data file from Siemens Scanner
rawData = mapVBVDVE(fName,'ignoreSeg'); % Ignore segment here to save memory
disp('---------------- Data dimensions ---------------- ');
disp(rawData{end}.image.sqzDims);

data = 0;
for i = 1:1
    if i > 1 && i*rawData{end}.hdr.Config.NSeg > size(data,3)
        break;
    elseif ~isfield(rawData{end}.hdr.Config, 'NSeg')
%         data = squeeze(rawData{2}.image(''));
        data =  squeeze(rawData{end}.image(:,:,:,1,1,1,:,1,1,1,1,1,1,1,1,1)); % To read only one slice
        break;
    end
    if sgI
        data = rawData{end}.image.unsorted;
    else
        data = squeeze(rawData{end}.image(''));
%     tmp = squeeze(rawData{2}.image(:,:,:,1,1,1,:,1,1,:,i,1,1,1,1,1)); % Read partial matrix to save memory
%     data = data+tmp; % ???????????? only for segmented acquisiton
    end
end

% Read noise if present
if isfield(rawData{2}, 'noise')
    nScan = rawData{2}.noise();
    % Remove piot tone from last line of noise scan
    nScan(end,:,:) = 0;
    % Also remove first 3 points as there is some extraneous signal present
    nScan(1:3,:,:) = 0;
else
    nScan =[];
end

rawData{2}.image

% Acquisition parameters
param.dim       = rawData{2}.image.sqzDims;                                     % Matrix dimension
param.slPos     = unique((rawData{2}.image.slicePos)', 'rows','stable')';       % Slice posiition
param.slThck    = rawData{2}.hdr.MeasYaps.sSliceArray.asSlice{1}.dThickness;    % Slice thickness

if isempty(rawData{2}.hdr.Dicom.flPhaseOS)
    param.FOV       = [(rawData{2}.hdr.Dicom.flReadoutOSFactor*rawData{2}.hdr. ...
        Config.RoFOV(1))*(1-sum(p.fd)), rawData{2}.hdr.Config.PeFOV(1)];    % FOV in mm ??? is oversampling always 2
else
    param.FOV       = [(rawData{2}.hdr.Dicom.flReadoutOSFactor*rawData{2}.hdr. ...
        Config.RoFOV)*(1-sum(p.fd)), rawData{2}.hdr.Config.PeFOV*(1+rawData{2}.hdr.Dicom.flPhaseOS)];    % FOV in mm ??? is oversampling always 2
end
param.matSz     = [2*(rawData{2}.image.centerCol(1)-1), 2*(rawData{2}.image.centerLin(1)-1), rawData{2}.image.NPar]; %[rawData{2}.image.NCol, rawData{2}.image.NLin, rawData{2}.image.NPar]; % image size
if isfield(rawData{2}.hdr.Meas, 'aflMaxwellCoefficients')
    param.mxCoef    = rawData{2}.hdr.Meas.aflMaxwellCoefficients;                   % Maxwell coefficients
else
    param.mxCoef = [];
end
param.patPos    = rawData{2}.hdr.Config.PatientPosition;                        % Patient position
param.TE        = rawData{2}.hdr.MeasYaps.alTE;
param.TRes      = rawData{2}.hdr.MeasYaps.alTR{1};
param.TR        = param.TRes/rawData{2}.hdr.Config.NSeg/rawData{2}.hdr.Config.NSetMeas; 
% param.TR        = rawData{2}.hdr.Meas.lEchoSpacing;
param.FA        = rawData{2}.hdr.Dicom.adFlipAngleDegree;
param.cCol      = mean(rawData{2}.image.centerCol);                             % central column; used to fix the center along FE
param.cRow      = mean(rawData{2}.image.centerLin);                             % central column; used to fix the center along FE

if isfield(rawData{2}.hdr.Config, 'DwellTime')
    param.BW = 1/(rawData{2}.hdr.Config.DwellTime * 1e-9); % Note, the param.BW = BW_scanner * Npixels; where Npixels is the number of pixels (with oversampling) in FE
elseif isfield(rawData{2}.hdr.Meas, 'alDwellTime')
    param.BW = 1/(rawData{2}.hdr.Meas.alDwellTime(1) * 1e-9);
elseif isfield(rawData{2}.hdr.MeasYaps.sRXSPEC, 'alDwellTime')
    param.BW = 1/(rawData{2}.hdr.MeasYaps.sRXSPEC.alDwellTime{1} * 1e-9);
elseif isfield(rawData{2}.hdr.Phoenix.sRXSPEC, 'alDwellTime')
    param.BW = 1/(rawData{2}.hdr.Phoenix.sRXSPEC.alDwellTime{1} * 1e-9);
else
    warning('BW information cannot be located; using a generic value of 400 kHz');
    param.BW = 400000;
end


% header = hdr();

%% Flow/phase-contrast
% ??????????????????????? uncomment this flow
if flow
tmp = rawData{2}.hdr.MeasYaps.sAngio.sFlowArray.lSize;
for i=1:tmp
    rawData{2}.hdr.MeasYaps.sAngio.sFlowArray.asElm{i};
    param.venc(i) = rawData{2}.hdr.MeasYaps.sAngio.sFlowArray.asElm{i}.nVelocity; % Venc values
    param.nDir(i) = rawData{2}.hdr.MeasYaps.sAngio.sFlowArray.asElm{i}.nDir;      % Direction of encoding
end
end
% ???????????????????????????

%% Error checks
if ~isempty(samp) && numel(sDim) ~= ndims(samp)
    error('Dimension mismatch between samp and sDim');
elseif numel(setdiff(sDim, param.dim))
    error('one of the elemens of sDim is not presnet in data');
end

%% Extract pilot tone from last kspace sample in for each readout and channel
% Extract pilot raw pilot tone
pilot_tone_raw = squeeze(data(end,:,:));
% Replace last sample with zero
data(end,:,:) = 0;
% Also remove first 3 points as there is some extraneous signal present
data(1:3,:,:) = 0;
%% Do you want to average out the 'seg' dimension. Using ignoreseg in mapVBVD makes this redundant
% if sgAvg == 1 && ~isempty(find(strcmp(param.dim, 'Seg'),1))
%     sgDim = find(strcmp(param.dim, 'Seg'));
%     data = squeeze(sum(data,sgDim));
%     param.dim(sgDim) = [];
%     disp('---------- Dimension "seg" is averaged ----------');
% end


%% Adjust zero-pading in FE direction so that peak occurs in the center
disp('Padding asymmetric echo...');
z1 = 0; z2 = 0;
MaxInd = double(param.cCol); % the location of peak in k-space
if MaxInd<(floor(size(data,1)/2)+1)
    z1 = 2*ceil((size(data,1)+2-2*MaxInd)/2);
    data = padarray(data,z1,0,'pre');
elseif MaxInd>(floor(size(data,1)/2)+1)
    z2 = 2*ceil((2*MaxInd-2-size(data,1))/2);
    data = padarray(data,z2,0,'post');
end

% ********
% data(1:52,:,:) = 0;
% ********

%% Adjust the PE dimension
% lnDim = find(strcmp(param.dim, 'Lin'));
% z3 = 2*rawData{2}.image.centerLin(1)-1-2*(ceil(rawData{2}.image.centerLin(1)/2)-rawData{2}.image.centerLin(1)/2) - size(data,lnDim);
% if z3 > 0
%     tmp = zeros(1,ndims(data));
%     tmp(lnDim) = 1; 
%     data = padarray(data, tmp*double(z3), 0, 'post');    % missing "unacquired" PE
% end

%% Find Maxwell Correction Maps
cMaps = 0;
N = size(data);
if mxFlg == 1
    disp('Finding Maxwell maps...');
    cMaps = mxCorr(data,param);
    cMaps(1:round(N(1)*fd(1)),:,:,:,:,:,:,:,:,:)=[]; % I am assuming, ndims(data) <= 10 ??
    cMaps(end-round(N(1)*fd(2))+1:end,:,:,:,:,:,:,:,:,:)=[]; %??
end

%% Cropping/tappering FE directions (optional)
if fd(1) ~= 0 || fd(2) ~= 0
disp(['Cropping FE dimension by ', num2str(100*(fd(1)+fd(2))), '%']);
data = single(data);
I = fftshift(ifft(ifftshift(data,1),[],1),1)* sqrt(N(1));
I(1:round(N(1)*fd(1)),:,:,:,:,:,:,:,:,:)=[]; % I am assuming, ndims(data) <= 10 ??
I(end-round(N(1)*fd(2))+1:end,:,:,:,:,:,:,:,:,:)=[]; % ??
% tmp = (cos(linspace(-pi, pi, 2*tap)')+1)/2; % These three lines for tapering the FE
% tmp = [tmp(1:end/2); ones(size(I,1)-numel(tmp),1); tmp(end/2+1:end)];
% I = bsxfun(@times, I, tmp);
data = fftshift(fft(ifftshift(I,1),[],1),1)/sqrt(size(I,1));% / sqrt(N(1)/(N(1)-z1-z2));
disp(['FE dimension cropped from ', num2str(N(1)), ' to ', num2str(size(data,1))]);
N=size(data);
end

%% Find noise std
disp('Finding noise standard deviation...');
chDim = find(strcmp(param.dim, 'Cha'));
if kNoise %|| isempty(nScan) % if separate noise scan does not exist
    warning('Noise is estimated from the k-space data');
    tmpStd = data([1:min(16,round(size(data,1)/8)), end-min(16,round(size(data,1)/8))+1:end],:,:,:,:,:,:,:,:,:);
    tmpStd = tmpStd(:,:,[1:round(size(data,3)/3), end-round(size(data,3)/3)+1:end],:,:,:,:,:,:,:);
    mStd = zeros(size(data,2),1);
    for i = 1:size(data,2)
        tmp = tmpStd(:,i,:,:,:);
        tmp = std(tmp(tmp~=0)); 
        mStd(i) = tmp; % Noise in the ith coil
    end
    tmp = prctile(sort(mStd), 25)/1.2; % pick a coil close to the bottom in terms of signal power
    data = data*nStd/tmp;
%     mStd = median(abs(mStd))/sqrt(0.675) / 2; % This factor of "1.5 to 3" is included to compensate for the signal that is present in outer parts of the kspace
%     mStd = repmat(mStd,[size(data,chDim),1]);
else % if separate noise scan does exist
    if nNrm == 0
        nScan = squeeze(nScan).* sqrt(param.BW/(size(nScan,1)/7.68e-3)); % Note, for prescan, readout is 7.68 ms long;
        NS=size(nScan);
        I = fftshift(ifft(ifftshift(nScan,1),[],1),1)*sqrt(NS(1));
        I(1:round(NS(1)*fd(1)),:,:,:,:,:,:,:,:,:)=[]; % I am assuming, ndims(data) <= 10 ??
        I(end-round(NS(1)*fd(2))+1:end,:,:,:,:,:,:,:,:,:)=[]; % ??
        nScan = fftshift(fft(ifftshift(I,1),[],1),1)*1/sqrt(size(I,1));% / sqrt(N(1)/(N(1)-z1-z2));
        
        nScan = permute(nScan, [2,1,3,4,5,6,7,8]);
        nScan = reshape(nScan, [size(nScan,1), numel(nScan)/size(nScan,1)]);
%         mStd = std(nScan,[],2);
        tmp = 1/(size(nScan,2)-1) * (nScan*nScan'); % covriance matrix
        tmp = (tmp + tmp')/2; % To make the matrix complex hermitian
%         tmp = diag(tmp); tmp = diag(tmp,0); % activate this line to ignore cross-correlation among coil elements
        L = chol(tmp, 'lower');
%         figure; imagesc(abs(L)); colobar;
        L_inv = inv(L);
        mStd = diag(sqrt(abs(tmp)));
        data =  permute(data, [2,1,3:numel(N)]);
        data = reshape(data, [size(data,1), numel(data)/size(data,1)]);
        data = L_inv * data;
        data = reshape(data, N([2,1,3:numel(N)]));
        data = permute(data, [2,1,3:numel(N)]);
        data = data * nStd; %bsxfun(@times, data, nStd./permute(mStd,circshift((1:ndims(data))', chDim-1))); 

    elseif nNrm == 1
        mStd = nScan.* sqrt(param.BW/(size(nScan,1)/7.68e-3)); % Note, for prescan, readout is 7.68 ms long;
        mStd = std(mStd(:));
%         mStd = mStd .* sqrt(param.BW/(size(data,1)/7.68e-3));
        mStd = repmat(mStd,[size(data,find(strcmp(param.dim, 'Cha'))),1]);
        data = bsxfun(@times, data, nStd./permute(mStd,circshift((1:ndims(data))', chDim-1))); 
    else
        error('wrong value assigned to p.nNrm');
    end
end
clear NS I;
disp(['Noise Std: ' sprintf('%0.2e  ',mStd')]);
%     [mean(mStd); median(mStd)]


%% Normalizing noise in each channel
% data = bsxfun(@times, data, nStd./permute(mStd,circshift((1:ndims(data))', chDim-1))); 
if nAmp > 1
	nzInd = data ~=0; % index with nonzero entries
	data(nzInd) = data(nzInd) + sqrt(nAmp^2-1)* nStd/sqrt(2)*single(   randn(sum(nzInd(:)),1)); % broken down into two steps to save memory
	data(nzInd) = data(nzInd) + sqrt(nAmp^2-1)* nStd/sqrt(2)*single(1j*randn(sum(nzInd(:)),1));
	data = data/nAmp;
elseif nAmp < 1
	error('value assigned to nAmp should be greater or equal to one');
end


%% Zero padding
% Circularly shift the PE direction to reverse wraparound
% This is creating small, non-zero entries at unsampled locations
% ???????????????????????????????
% peDim = find(strcmp(param.dim, 'Lin'));
% data = fftshift(ifft(ifftshift(data,peDim),[],peDim),peDim);
% tmp = zeros(1,8); tmp(peDim) = peSft;
% data = circshift(data, tmp);
% data = fftshift(fft(ifftshift(data,peDim),[],peDim),peDim);

% z1 = 0; z2 = 0;
% MaxInd = round(double(param.cCol)*(1-sum(fd))); % the location of peak in k-space
% if MaxInd<(floor(size(data,1)/2)+1)
%     z1 = 2*ceil((size(data,1)+2-2*MaxInd)/2);
%     data = padarray(data,z1,0,'pre');
% elseif MaxInd>(floor(size(data,1)/2)+1)
%     z2 = 2*ceil((2*MaxInd-2-size(data,1))/2);
%     data = padarray(data,z2,0,'post');
% end


%% Match the sampling matrix dimensions with that of data, and apply downsampling
disp('Matching sampling matrix to data...')
sDimInd = zeros(numel(sDim),1);
for i=1:numel(sDim), sDimInd(i) = find(strcmp(param.dim, sDim(i))); end

if ~isempty(samp)
    pDim = zeros(ndims(data),1); % Permutation used to match 'samp' with 'data'
    for i=1:numel(sDim), pDim(sDimInd(i))=i; end
    for i=find(pDim==0)', pDim(i) = max(pDim(:))+1; end
    samp = permute(samp,pDim);
    samp = bsxfun(@times, true(size(data)), samp);
    samp = logical(sum(samp,chDim));
    data = bsxfun(@times, data, samp);
else
    dataTmp = sum(abs(data),chDim);
    samp = dataTmp>0;
end
clear dataTmp;

%% Channel compression ====================================================
% Performing SVD for Channel Compression (optional) 
disp('Compressing channels...');
N=size(data);
% N(chDim),
Nc = ceil(max(Nc(1), Nc(2)*N(chDim)));
dataC = dimData(data,[],chDim, [1:min(Nc,N(chDim))],'r'); % Initialization
slDim = find(strcmp(param.dim, 'Sli')); % Slice dimension
if isempty(slDim), slDim = ndims(data)+1;end
nPer = 1:max(ndims(data),slDim);nPer(chDim)=nPer(end); nPer(end)=chDim;
if Nc < N(chDim) % only compress if virtual channels are lesser than actual channels
    for i=1:size(data,slDim)
        tmp = dimData(data,[],slDim,i,'r'); % Read data along 'slDim'
        tmp = permute(tmp,nPer);
        NTmp=size(tmp);
        tmp = reshape(tmp,[prod(NTmp(1:end-1)),NTmp(end)]);
        [W,S,~] = svd(tmp,0);
        ST = S*eye(size(data,chDim),Nc);            % Truncated Sigma

        % First "Nc" orthogonal channels
        tmp = W*ST; % Make coil orthonormal by W*(ST>0)
        tmp = reshape(tmp,[NTmp(1:end-1),Nc]);
        tmp = permute(tmp, nPer);
        dataC = dimData(dataC,tmp, slDim,i,'w'); % Write data along 'slDim'
    end
    data = dataC;
end
% size(data)
% adfa
clear dataC;
% adfa

%% Re-enforcing the partial Fourier which may have gotten lost due to lack of numerical precision
disp('Reinforcing partial Fourier...');
if z1>0
%     data = padarray(data,2*round(z1*(1-fd(1)-fd(2))/2),0,'pre');
%     samp = padarray(samp,2*round(z1*(1-fd(1)-fd(2))/2),0,'pre');
    data(1:ceil(z1),:,:,:,:,:,:,:,:,:) = 0;
    samp(1:ceil(z1),:,:,:,:,:,:,:,:,:) = 0;
elseif z2>0
%     data = padarray(data,2*round(z2*(1-fd(1)-fd(2))/2),0,'post');
%     samp = padarray(samp,2*round(z2*(1-fd(1)-fd(2))),0,'post');
    data(end:-1:end-ceil(z2)+1,:,:,:,:,:,:,:,:,:) = 0;
    samp(end:-1:end-ceil(z2)+1,:,:,:,:,:,:,:,:,:) = 0;
end
data = bsxfun(@times, data, samp); % apply sampling again to kill small non-zero 
                                   % values created due to limited precision of SVD

rawData = rawData{end};


disp('Done!');


%% Plotting ===============================================================
% % This part is hard-coded, may not work if the structure is different
% if dis==1
%     for i=1:size(data,slDim)
%         dataTmp = dimData(data,[],slDim,i,'r');
%         k_ind = max(abs(dataTmp),[],1)>eps; 
%         k_ind = repmat(k_ind,[size(dataTmp,1),1,1,1]);
%         k_ind = sum(sum(sum(sum(sum(sum(sum(k_ind,4),5),6),7),8),9),10);
%         kAvg =  sum(sum(sum(sum(sum(sum(sum(dataTmp,4),5),6),7),8),9),10)./(k_ind+eps);
%         permSq = [find(strcmp(param.dim, 'Col')), find(strcmp(param.dim, 'Lin')), find(strcmp(param.dim, 'Cha'))];
%         kAvg = permute(kAvg,permSq);
%         slicer(sqrt(abs(fftshift(ifft2(ifftshift(kAvg))))),1,1.0, ['time-averaged coil, slice ' num2str(i)]);
%         figure; imagesc(squeeze(sqrt(sum(abs(fftshift(ifft2(ifftshift(kAvg)))).^2,3)))); 
%         axis('image'); title(['time-averaged SoS, slice ' num2str(i)]); %colormap(gray);
%     end
% end


% %% Error checks
% if ~isempty(samp) && numel(sDim) ~= ndims(samp)
%     error('Dimension mismatch between samp and sDim');
% elseif numel(setdiff(sDim, param.dim))
%     error('one of the elemens of sDim is not presnet in data');
% end
% 
% %% Do you want to average out the 'seg' dimension. Using ignoreseg in mapVBVD makes this redundant
% if sgAvg == 1 && ~isempty(find(strcmp(param.dim, 'Seg'),1))
%     sgDim = find(strcmp(param.dim, 'Seg'));
%     data = squeeze(sum(data,sgDim));
%     param.dim(sgDim) = [];
%     disp('---------- Dimension "seg" is averaged ----------');
% end
% 
% 
% %% Adjust zero-pading in FE direction so that peak occurs in the center
% z1 = 0; z2 = 0;
% MaxInd = double(param.cCol); % the location of peak in k-space
% if MaxInd<(floor(size(data,1)/2)+1)
%     z1 = 1*ceil((size(data,1)+2-2*MaxInd)/1);
%     data = padarray(data,z1,0,'pre');
% elseif MaxInd>(floor(size(data,1)/2)+1)
%     z2 = 1*ceil((2*MaxInd-2-size(data,1))/1);
%     data = padarray(data,z2,0,'post');
% end
% 
% MaxInd = double(param.cRow); % the location of peak in k-space
% if MaxInd<(floor(size(data,3)/2)+1)
%     z3 = 1*ceil((size(data,3)+2-2*MaxInd)/1);
%     data = padarray(data,[0,0,z3,0,0],0,'pre');
% elseif MaxInd>(floor(size(data,3)/2)+1)
%     z4 = 1*ceil((2*MaxInd-2-size(data,3))/1);
%     data = padarray(data,[0,0,z4,0,0],0,'post');
% end
% 
% % break;
% 
% %% Find Maxwell Correction Maps
% cMaps = 0;
% N = size(data);
% if mxFlg == 1 && ~isempty(param.mxCoef)
%     cMaps = mxCorr(data,param);
%     cMaps(1:2*round(N(1)*fd(1)/2),:,:,:,:,:,:,:,:,:)=[]; % I am assuming, ndims(data) <= 10
%     cMaps(end-2*round(N(1)*fd(2)/2)+1:end,:,:,:,:,:,:,:,:,:)=[];
% end
% 
% %% Find noise std
% chDim = find(strcmp(param.dim, 'Cha'));
% if isempty(nScan) || kNoise% if separate noise scan does not exist
%     warning('Separate noise scan does not exist');
%     mStd = data([1:min(10,round(size(data,1)/10)), end-min(10,round(size(data,1)/10))+1:end],:,:,:,:,:,:,:,:,:);
%     mStd = mStd(:,:,[1:max(8,round(size(data,3)/3)), end-max(8,round(size(data,3)/3))+1:end],:,:,:,:,:,:,:);
%     mStd = mStd(mStd~=0);
%     mStd = median(abs(mStd))/sqrt(0.675) / 2; % This factor of "1.5 to 3" is included to compensate for the signal that is present in outer parts of the kspace
%     mStd = repmat(mStd,[size(data,chDim),1]);
% else % if separate noise scan does exist
%     if nNrm == 0
%         nScan = squeeze(nScan);
%         nScan = permute(nScan, [2,1,3,4,5,6,7,8]);
%         nScan = reshape(nScan, [size(nScan,1), numel(nScan)/size(nScan,1)]);
%         mStd = std(nScan,[],2);
%         nInd = find(mStd>3*median(mStd));
%         if numel(nInd)>0
%             warning([num2str(numel(nInd)), ' noisy coils dicarded']);
%             data(:, nInd, :, :, :, :, :, :) = [];
%             mStd(nInd) = [];
%         end
%     elseif nNrm == 1
%         mStd = nScan;
%         mStd = std(mStd(:));
%         mStd = repmat(mStd,[size(data,find(strcmp(param.dim, 'Cha'))),1]);
%     else
%         error('wrong value assigned to p.nNrm');
%     end
%     disp('Noise std in channels:');
%     mStd',
% %     [mean(mStd); median(mStd)]
% end
% 
% 
% %% Normalizing noise in each channel
% data = bsxfun(@times, data, nStd./permute(mStd,circshift((1:ndims(data))', chDim-1))); 
% if nAmp > 1
% 	nzInd = data ~=0; % index with nonzero entries
% 	data(nzInd) = data(nzInd) + sqrt(nAmp^2-1)* nStd/sqrt(2)*single(   randn(sum(nzInd(:)),1)); % broken down into two steps to save memory
% 	data(nzInd) = data(nzInd) + sqrt(nAmp^2-1)* nStd/sqrt(2)*single(1j*randn(sum(nzInd(:)),1));
% 	data = data/nAmp;
% elseif nAmp < 1
% 	error('value assigned to nAmp should be greater or equal to one');
% end
% 
% %% Cropping FE directions (optional)
% data = single(data);
% N=size(data);
% I = fftshift(ifft2(ifftshift(data)));
% I(1:2*round(N(1)*fd(1)/2),:,:,:,:,:,:,:,:,:)=[]; % I am assuming, ndims(data) <= 10
% I(end-2*round(N(1)*fd(2)/2)+1:end,:,:,:,:,:,:,:,:,:)=[];
% data = fftshift(fft2(ifftshift(I)));
% % N=size(data),
% 
% 
% % size(data),
% % sdfasdf
% %% Match the sampling matrix dimensions with that of data, and apply downsampling
% sDimInd = zeros(numel(sDim),1);
% for i=1:numel(sDim), sDimInd(i) = find(strcmp(param.dim, sDim(i))); end
% 
% if ~isempty(samp)
%     pDim = zeros(ndims(data),1); % Permutation used to match 'samp' with 'data'
%     for i=1:numel(sDim), pDim(sDimInd(i))=i; end
%     for i=find(pDim==0)', pDim(i) = max(pDim(:))+1; end
%     samp = permute(samp,pDim);
%     samp = bsxfun(@times, true(size(data)), samp);
%     samp = logical(sum(samp,chDim));
%     data = bsxfun(@times, data, samp);
% else
%     dataTmp = sum(abs(data),chDim);
%     samp = dataTmp>0;
% end
% clear dataTmp;
% 
% %% Channel compression ====================================================
% % Performing SVD for Channel Compression (optional) 
% N=size(data);
% % N(chDim),
% dataC = dimData(data,[],chDim, [1:min(Nc,N(chDim))],'r'); % Initialization
% slDim = find(strcmp(param.dim, 'Sli')); % Slice dimension
% if isempty(slDim), slDim = ndims(data)+1;end
% nPer = 1:max(ndims(data),slDim);nPer(chDim)=nPer(end); nPer(end)=chDim;
% if Nc < N(chDim) % only compress if virtual channels are lesser than actual channels
%     for i=1:size(data,slDim)
%         tmp = dimData(data,[],slDim,i,'r'); % Read data along 'slDim'
%         tmp = permute(tmp,nPer);
%         NTmp=size(tmp);
%         tmp = reshape(tmp,[prod(NTmp(1:end-1)),NTmp(end)]);
%         [W,S,~] = svd(tmp,0);
%         ST = S*eye(size(data,chDim),Nc);            % Truncated Sigma
% 
%         % First "Nc" orthogonal channels
%         tmp = W*ST; % Make coil orthonormal by W*(ST>0)
%         tmp = reshape(tmp,[NTmp(1:end-1),Nc]);
%         tmp = permute(tmp, nPer);
%         dataC = dimData(dataC,tmp, slDim,i,'w'); % Write data along 'slDim'
%     end
%     data = dataC;
% end
% % size(data)
% % adfa
% clear dataC;
% % adfa
% 
% %% Re-enforcing the partial Fourier which may have gotten lost due to lack of numerical precision
% if z1>0
% %     data = padarray(data,2*round(z1*(1-fd(1)-fd(2))/2),0,'pre');
% %     samp = padarray(samp,2*round(z1*(1-fd(1)-fd(2))/2),0,'pre');
%     data(1:ceil(z1*(1-fd(1)-fd(2))),:,:,:,:,:,:,:,:,:) = 0;
%     samp(1:ceil(z1*(1-fd(1)-fd(2))),:,:,:,:,:,:,:,:,:) = 0;
% elseif z2>0
% %     data = padarray(data,2*round(z2*(1-fd(1)-fd(2))/2),0,'post');
% %     samp = padarray(samp,2*round(z2*(1-fd(1)-fd(2))),0,'post');
%     data(end:-1:end-ceil(z2*(1-fd(1)-fd(2)))+1,:,:,:,:,:,:,:,:,:) = 0;
%     samp(end:-1:end-ceil(z2*(1-fd(1)-fd(2)))+1,:,:,:,:,:,:,:,:,:) = 0;
% end
% data = bsxfun(@times, data, samp); % apply sampling again to kill small non-zero 
%                                    % values created due to limited precision of SVD
% 
% 
% 
% %% Plotting ===============================================================
% % This part is hard-coded, may not work if the structure is different
% if dis==1
%     for i=1:size(data,slDim)
%         dataTmp = dimData(data,[],slDim,i,'r');
%         k_ind = max(abs(dataTmp),[],1)>eps; 
%         k_ind = repmat(k_ind,[size(dataTmp,1),1,1,1]);
%         k_ind = sum(sum(sum(sum(sum(sum(sum(k_ind,4),5),6),7),8),9),10);
%         kAvg =  sum(sum(sum(sum(sum(sum(sum(dataTmp,4),5),6),7),8),9),10)./(k_ind+eps);
%         permSq = [find(strcmp(param.dim, 'Col')), find(strcmp(param.dim, 'Lin')), find(strcmp(param.dim, 'Cha'))];
%         kAvg = permute(kAvg,permSq);
%         slicer(sqrt(abs(fftshift(ifft2(ifftshift(kAvg))))),1,1.0, ['time-averaged coil, slice ' num2str(i)]);
%         figure; imagesc(squeeze(sqrt(sum(abs(fftshift(ifft2(ifftshift(kAvg)))).^2,3)))); 
%         axis('image'); title(['time-averaged SoS, slice ' num2str(i)]); %colormap(gray);
%     end
% end

