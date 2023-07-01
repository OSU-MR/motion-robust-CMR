function [Outputs] = binSelfGating(Dat, rawData, param, opt)

Lin(1,1,:) = rawData.image.Lin; % Phase encoding indices
Par(1,1,:) = rawData.image.Par; % Partition encoding indices
Set(1,1,:) = rawData.image.Set; % Velocity encoding indices

% ***************
% Only for 3D Cine
if ~opt.flow
Dat = Dat(:,:,7:end);
Lin = Lin(1,1,7:end);
Par = Par(1,1,7:end);
Set = Set(1,1,7:end);
end
% *******************
opt.timestamp = rawData.image.timestamp;
% *****************************************
% debugSamplingPattern(Lin, Par, Set, opt, rawData);
% *****************************************

ArraySize(1) = size(Dat,1);
ArraySize(2) = rawData.image.sqzSize(2);
ArraySize(3) = rawData.image.sqzSize(3);
ArraySize(4) = rawData.image.sqzSize(4);
ArraySize(5) = max(Set(:));
ArraySize(6) = opt.nCPhases + opt.viewSharing*opt.nCPhases;
ArraySize(7) = opt.nRPhases;

opt.ArraySize = ArraySize;

%% EXTRACT RESPIRATORY AND CARDIAC MOTION SURROGATES
% Extract self-gating signal
 if strcmp(opt.seq, 'cart')
     SG  = Dat(:,:,opt.sgI:opt.sgI:end);
     opt.TR = param.TR*opt.sgI/1e6;
    % Self-gating sampling frequency in Hz
    opt.fs = 1/opt.TR;   
    opt.timeFull = [0:1:size(Dat,3)-1]*(opt.TR/opt.sgI);
    opt.timeFull = opt.timeFull(opt.sgI:opt.sgI:end);
    opt.timestamp = opt.timestamp(opt.sgI:opt.sgI:end);
 else
     SG  = Dat(:,:,1:opt.sgI:end);
     % Self-gating temporal resolution in milliseconds
    opt.TR = param.TR*opt.sgI/1e6;
    % Self-gating sampling frequency in Hz
    opt.fs = 1/opt.TR;   
    opt.timeFull = [0:1:size(Dat,3)-1]*(opt.TR/opt.sgI);
    opt.timeFull = opt.timeFull(1:opt.sgI:end);
    opt.timestamp = opt.timestamp(1:opt.sgI:end);
 end


opt.cardiac_full = zeros(1, length(SG));
opt.resp_full = zeros(1, length(SG));
opt.signal_index = zeros(1,length(SG));

% Trim trailing edge of data, Lin, Par, and Set in case a full cycle wasn't completed
Dat = Dat(:,:,1:size(SG,3)*opt.sgI);
Lin = Lin(:,:,1:size(SG,3)*opt.sgI);
Par = Par(:,:,1:size(SG,3)*opt.sgI);
Set = Set(:,:,1:size(SG,3)*opt.sgI);

% Trim samples during steady-state ramp-up (maybe make this based on time?)
SG  = SG(:,:,opt.rampSS+1:end);
Dat = Dat(:,:,1+opt.sgI*opt.rampSS:end);
Lin = Lin(:,:,1+opt.sgI*opt.rampSS:end);
Par = Par(:,:,1+opt.sgI*opt.rampSS:end);
Set = Set(:,:,1+opt.sgI*opt.rampSS:end);

% Extract cardiac and respiratory motion
% disp('Extracting motion surrogates...');
% [R, C, respRange] = extractMotionSignalsWrapper(SG,opt);

% ****************************
% Temporary for 3D Cine (Sola)
% Solve for quadratic drift in self-gating respiratory signal
% t = 0 : length(R)-1;
% A(:,1) = ones(length(R), 1);
% A(:,2) = t;
% A(:,3) = t.^2;
% A(:,4) = t.^3;
% x = (A'*A)\(A'*R);
% R = R - A*x;
% %*****************************

% Only reconstruct 'recTime' worth of data (e.g. 5, 10, 20 minutes)
if ~isinf(opt.recTime)
    if ~isnumeric(opt.recStart)
        switch opt.recStart
            case 'start'
                firstSample = 1;
                lastSample = min(ceil(opt.recTime / opt.TR), size(SG,3));
            case 'center'
                firstSample = max(round(size(SG,3)/2) - ceil(0.5 * opt.recTime / opt.TR), 1);
                lastSample = min(round(size(SG,3)/2) + ceil(0.5 * opt.recTime / opt.TR), size(SG,3));
            case 'end'
                firstSample = max(size(SG,3) - ceil(opt.recTime / opt.TR), 1);
                lastSample = size(SG,3);
        end
    else
        firstSample = max(ceil(opt.recStart / opt.TR), 1);
        lastSample = min(firstSample + ceil(opt.recTime / opt.TR), size(SG,3));
    end
    
    opt.firstSample = firstSample;
    opt.lastSample = lastSample;
    FlowSG4D_Outputs.firstSample = (firstSample*opt.sgI)-(opt.sgI-1);
    FlowSG4D_Outputs.lastSample = lastSample*opt.sgI;
    
    SG = SG(:,:,firstSample:lastSample);
%     R = R(firstSample:lastSample,:);
%     C = C(firstSample:lastSample,:);
    Dat = Dat(:,:,(firstSample*opt.sgI)-(opt.sgI-1):lastSample*opt.sgI);
    Lin = Lin(:,:,(firstSample*opt.sgI)-(opt.sgI-1):lastSample*opt.sgI);
    Par = Par(:,:,(firstSample*opt.sgI)-(opt.sgI-1):lastSample*opt.sgI);
    Set = Set(:,:,(firstSample*opt.sgI)-(opt.sgI-1):lastSample*opt.sgI);
    
    opt.signal_index([opt.rampSS + firstSample, opt.rampSS + lastSample]) = 1;
    
end     

% Extract cardiac and respiratory motion
disp('Extracting motion surrogates...');
[R, C, respRange] = extractMotionSignalsWrapper(SG,opt);

% Ensure self-gating and data are lined up
if Dat(:,:,opt.sgI:opt.sgI:end) ~= SG
    warning('Data not sychronized with self-gating');
end

%% K-SPACE BINNING
% Total number of lines collected per coil
lengthDat = size(Dat, 3);          

% % Respiratory phase indices
% disp(['Binning into ',num2str(opt.nRPhases), ' respiratory phases...']);
% % [rBins(1,1,:) ] = binReadoutLinesByRespiratoryPhase(lengthDat,R,opt);
% % [rBins(1,1,:), rWeights(1,1,:)] = binRespiratoryPhases(lengthDat,R,opt);
% [rBins(1,1,:), rWeights(1,1,:)] = binRespWrapper(lengthDat, R, opt);


% =========================================================================
% Either here or in the above function, want to implement a bulk motion
% rejection criteria. First need to validate bulk motion occurs. 
% =========================================================================

% Cardiac phase indices
disp(['Binning k-space into ',num2str(ArraySize(6)), ' cardiac phases...']);
% [cBins, meanRR, stdRR] = binReadoutLinesByCardiacPhase(lengthDat,C,opt);
[cBins(1,1,:), param.meanRR, param.stdRR, param.triggerTime, triggers] = binCardiacWrapper(lengthDat, C, opt);
opt.NumberOfCardiacPhasesToRemove = length(find(isnan(cBins)));
shift = opt.TR*(opt.rampSS);
% triggers = triggers + shift;
C_full = zeros(1,length(opt.timeFull)); 
C_full(end-length(C)+1:end) = C;

R_full = zeros(1,length(opt.timeFull)); 
R_full(end-length(R)+1:end) = R;

start = min(find(opt.signal_index)); stop = max(find(opt.signal_index));
C_tmp = orientCardiac_2(C, opt);
opt.cardiac_full(start:stop) = C_tmp;
opt.resp_full(start:stop) = R;
opt.triggers = triggers + opt.timeFull(start);
opt.trigger_pks = interp1(opt.timeFull, opt.cardiac_full, opt.triggers);

figure;
subplot(211); plot(opt.timeFull, opt.resp_full);
subplot(212); plot(opt.timeFull, opt.cardiac_full); hold on; plot(opt.triggers, opt.trigger_pks, '.r');

% *******************************************************
% Estimate number of repeat samples
Lintmp = Lin; Lintmp(:,:,opt.sgI:opt.sgI:end) = [];
Partmp = Par; Partmp(:,:,opt.sgI:opt.sgI:end) = [];
Settmp = Set; Settmp(:,:,opt.sgI:opt.sgI:end) = [];
cBinstmp = cBins; cBinstmp(:,:,opt.sgI:opt.sgI:end) = [];
TrimIndxC = ~isnan(cBinstmp);
Lintmp   = Lintmp(:,:,TrimIndxC);
Partmp   = Partmp(:,:,TrimIndxC);
Settmp   = Settmp(:,:,TrimIndxC);
cBinstmp = cBinstmp(:,:,TrimIndxC);

TMP(:,1) = Lintmp;
TMP(:,2) = Partmp;
TMP(:,3) = Settmp;
TMP(:,4) = cBinstmp;

[U, unique_ind] = unique(TMP, 'rows');
opt.repeatFrac = size(TMP,1)/size(U,1);
repeat_percent = 100*(size(TMP,1)-size(U,1))/size(TMP,1);
% *******************************************************

% param.respRange = (param.FOV(1) / rawData.hdr.Config.NColMeas)*respRange;
param.respRange = opt.resp_range;
disp(['Respiratory range is ',num2str(param.respRange), ' mm']);



% Respiratory phase indices
disp(['Binning into ',num2str(opt.nRPhases), ' respiratory phases...']);
% [rBins(1,1,:) ] = binReadoutLinesByRespiratoryPhase(lengthDat,R,opt);
% [rBins(1,1,:), rWeights(1,1,:)] = binRespiratoryPhases(lengthDat,R,opt);
[rBins(1,1,:), rWeights(1,1,:,:),respEfficiency] = binRespWrapper(lengthDat, R, opt, param);

% Before moving on, ensure sizes are all equal
if min([size(Dat,3),size(Lin,3),size(Par,3),size(Set,3),size(rBins,3),size(cBins,3)]) ~= ...
   max([size(Dat,3),size(Lin,3),size(Par,3),size(Set,3),size(rBins,3),size(cBins,3)])
    warning('Data dimensions are inconsistent');
end

% **************************************************************
% Estimate number of repeat samples (weights)
Weighttmp = rWeights; Weighttmp(:,:, opt.sgI:opt.sgI:end) = [];
Weighttmp = Weighttmp(:,:,TrimIndxC);
Unique_Weights = Weighttmp(:,:,unique_ind);
repeat_percent_weights = 100*(sum(Weighttmp(:).^2) - sum(Unique_Weights(:).^2))/sum(Weighttmp(:));
% **************************************************************



%% OUTPUT
% Trim out self-gating signal
Dat(:,:,opt.sgI:opt.sgI:end)    = [];
Lin(:,:,opt.sgI:opt.sgI:end)    = [];
Par(:,:,opt.sgI:opt.sgI:end)    = [];
Set(:,:,opt.sgI:opt.sgI:end)    = [];
rBins(:,:,opt.sgI:opt.sgI:end)  = [];
for i = 1 : opt.nRPhases
    tmp = rWeights(:,:,:,i);
    tmp(:,:,opt.sgI:opt.sgI:end) = [];
    rWeights0(:,:,:,i) = tmp;
end
rWeights = rWeights0;
% rWeights(:,:,opt.sgI:opt.sgI:end) = [];
cBins(:,:,opt.sgI:opt.sgI:end)  = [];

% Trim out unused points from cardiac binning and respiratory binning
TrimIndxR = ~isnan(rBins);
TrimIndxC = ~isnan(cBins);
TrimIndx = logical(TrimIndxR .* TrimIndxC);
Dat   = Dat(:,:,TrimIndx);
Lin   = Lin(:,:,TrimIndx);
Par   = Par(:,:,TrimIndx);
Set   = Set(:,:,TrimIndx);
rBins = rBins(:,:,TrimIndx);
rWeights = rWeights(:,:,TrimIndx,:);
cBins = cBins(:,:,TrimIndx);



% % Remove cardiac bins 12-20 (test for GPU recon.)
% cc = opt.c(c,:);
% TrimIndx = find(ismember(cBins,cc(1):cc(2)));
% Dat   = Dat(:,:,TrimIndx);
% Lin   = Lin(:,:,TrimIndx);
% Par   = Par(:,:,TrimIndx);
% Set   = Set(:,:,TrimIndx);
% rBins = rBins(:,:,TrimIndx);
% cBins = cBins(:,:,TrimIndx);
% opt.nCPhases = cc(2)-cc(1)+1;

% 
Dat = Dat;
Ind(:,1) = Lin(:);      % Phase encoding
Ind(:,2) = Par(:);      % Partition encoding
Ind(:,3) = Set(:);      % Velocity encoding
Ind(:,4) = cBins(:);    % Cardiac phase
Ind(:,5) = rBins(:);    % Respiratory phase

ArraySize(1) = size(Dat,1);
ArraySize(2) = rawData.image.sqzSize(2);
ArraySize(3) = rawData.image.sqzSize(3);
ArraySize(4) = rawData.image.sqzSize(4);
ArraySize(5) = max(Set(:));
ArraySize(6) = opt.nCPhases + opt.viewSharing*opt.nCPhases;
ArraySize(7) = opt.nRPhases;

ArrayNames{1} = 'Readout';
ArrayNames{2} = 'Channels';
ArrayNames{3} = 'Phase encode';
ArrayNames{4} = 'Partitions';
ArrayNames{5} = 'Velocity encodings';
ArrayNames{6} = 'Cardiac phases';
ArrayNames{7} = 'Respiratory phases';

% Make robust for both VD and VE
% *************************************************************************
FlowSG4D_Outputs.Cardiac = C;
FlowSG4D_Outputs.Respiratory = R;
FlowSG4D_Outputs.Time = [-length(C)/2 : 1 : length(C)/2-1] * opt.TR;
FlowSG4D_Outputs.Time_Full = opt.timeFull;
FlowSG4D_Outputs.C_Full = C_full;
FlowSG4D_Outputs.R_Full = R_full;
FlowSG4D_Outputs.triggers = triggers;

FlowSG4D_Outputs.Signals.Time = opt.timeFull;
FlowSG4D_Outputs.Signals.Stamp = opt.timestamp;
FlowSG4D_Outputs.Signals.C = opt.cardiac_full;
FlowSG4D_Outputs.Signals.R = opt.cardiac_full;
FlowSG4D_Outputs.Signals.TriggerTime = opt.triggers;
FlowSG4D_Outputs.Signals.TriggerPeaks = opt.trigger_pks;
FlowSG4D_Outputs.Signals.NumberofPeaks = numel(opt.trigger_pks);
FlowSG4D_Outputs.Signals.RR = param.meanRR;
FlowSG4D_Outputs.Signals.RR_std = param.stdRR;
FlowSG4D_Outputs.Signals.TRes = param.meanRR/opt.nCPhases;
FlowSG4D_Outputs.Signals.REff = opt.respEff;

FlowSG4D_Outputs.ArraySize = ArraySize;
FlowSG4D_Outputs.ArrayNames = ArrayNames;
FlowSG4D_Outputs.Data = Dat;
FlowSG4D_Outputs.Indices = Ind;
FlowSG4D_Outputs.rWeights = rWeights;
FlowSG4D_Outputs.param = param;
if strcmp(opt.version,'VE') || ~opt.flow
    FlowSG4D_Outputs.param.NColMeas = rawData.hdr.Config.NColMeas;
    FlowSG4D_Outputs.param.NImageCols = rawData.hdr.Config.NImageCols;
    FlowSG4D_Outputs.param.NParMeas = rawData.hdr.Config.NParMeas;
    FlowSG4D_Outputs.param.NImagePars = rawData.hdr.Config.NImagePar;
    FlowSG4D_Outputs.param.NLinMeas = rawData.hdr.Config.NLinMeas;
elseif strcmp(opt.version,'VA')
    FlowSG4D_Outputs.param.NColMeas = rawData.hdr.Config.NRoFTLenOriginal;
    FlowSG4D_Outputs.param.NImageCols = rawData.hdr.Config.NImageCols;
    FlowSG4D_Outputs.param.NParMeas = rawData.hdr.Config.NParMeas;
    FlowSG4D_Outputs.param.NImagePars = rawData.hdr.Config.NImagePars;
    FlowSG4D_Outputs.param.NLinMeas = rawData.hdr.Config.NLinMeas;
    FlowSG4D_Outputs.param.NImageLins = rawData.hdr.Config.NImageLins;
elseif strcmp(opt.version,'VD')
    FlowSG4D_Outputs.param.NColMeas = rawData.hdr.Meas.NRoFTLenOriginal;
    FlowSG4D_Outputs.param.NImageCols = rawData.hdr.Meas.NImageCols;
    FlowSG4D_Outputs.param.NParMeas = rawData.hdr.Meas.NParMeas;
    FlowSG4D_Outputs.param.NImagePars = rawData.hdr.Meas.NImagePars;
    FlowSG4D_Outputs.param.NLinMeas = rawData.hdr.Meas.NLinMeas;
    FlowSG4D_Outputs.param.NImageLins = rawData.hdr.Meas.NImageLins;
end
  
FlowSG4D_Outputs.param.NImageLins = rawData.hdr.Config.NImageLins;
FlowSG4D_Outputs.param.ResCol = param.FOV(1) / FlowSG4D_Outputs.param.NColMeas;
FlowSG4D_Outputs.param.ResLin = param.FOV(2) / FlowSG4D_Outputs.param.NLinMeas;
FlowSG4D_Outputs.param.ResPar = param.slThck / FlowSG4D_Outputs.param.NImagePars;
FlowSG4D_Outputs.param.RR.mean = param.meanRR;
FlowSG4D_Outputs.param.RR.std = param.stdRR;
FlowSG4D_Outputs.param.RR.Tres = param.meanRR / opt.nCPhases;
FlowSG4D_Outputs.param.triggerTime = param.triggerTime;
FlowSG4D_Outputs.param.respRange = respRange * FlowSG4D_Outputs.param.ResCol;
FlowSG4D_Outputs.param.respEfficiency = respEfficiency;
% *************************************************************************

Outputs = FlowSG4D_Outputs;

return
% ****************************
% IndiceNames = {'FE','PE1','PE2','Coils','Time','Vel'};
% lengthIndice = numel(Dat);
% FERange = [1:ArraySize(1)]';
% Indices.FE = repmat(FERange,lengthIndice/length(FERange),1);
% Indices.PE1 = 

counter = 0;
for i = 1 : size(Dat,3)
    tmpLin = Lin(i);
    tmpPar = Par(i);
    tmpcBins = cBins(i);
    tmpSet = Set(i);
    tmpWeights = rWeights(i);
    for j = 1 : size(Dat,2)
        for k = 1 : size(Dat,1)
            counter = counter + 1;
            FE(counter) = k;
            PE1(counter) = tmpLin;
            PE2(counter) = tmpPar;
            Coils(counter) = j;
            Cardiac(counter) = tmpcBins;
            Vel(counter) = tmpSet;
            Weights(counter) = tmpWeights;
            
        end
    end
end
Indices(:,1) = FE;
Indices(:,2) = PE1;
Indices(:,3) = PE2;
Indices(:,4) = Coils;
Indices(:,5) = Cardiac;
Indices(:,6) = Vel;

% *************************************************************************
% Want to look at repeat samples
% I.e. build matrix/array so ReVEAL4D can handle repeat samples
% kSpace: [FE x PE1 x PE2 x Coil x Time x VEL]
% samps: [1 x 1 x PE1 x PE2 x Time x Repeat] (logical)

samps = zeros(1, 1, ArraySize(3), ArraySize(4), ArraySize(5), ArraySize(6));
D = zeros(1, 1, ArraySize(3), ArraySize(4), ArraySize(5), ArraySize(6),'logical');

tmp = zeros(1, 1, ArraySize(3), ArraySize(4), ArraySize(5), ArraySize(6));
for i = 1 : size(Dat, 3)
    tmp(1, 1, Lin(1,1,i), Par(1,1,i), Set(1,1,i), cBins(1,1,i)) = tmp(1, 1, Lin(1,1,i), Par(1,1,i), Set(1,1,i), cBins(1,1,i)) + 1;
end
maxRep = max(max(max(max(max(max(tmp))))));


% Background encoding
% Put D in the form: [FE x PE1 x PE2 x Coils x Time x Rep]
DB =  zeros(ArraySize(1), ArraySize(3), ArraySize(4), ArraySize(2), ArraySize(6), maxRep,'uint32');
Test = zeros(ArraySize(1), ArraySize(3), ArraySize(4), ArraySize(2), ArraySize(6), maxRep,'logical');
indsB = find(Indices(:,6) == 1);
IndicesB = Indices(indsB,:);
DatB = Dat(:);
DatB = DatB(indsB);
WeightsB = Weights(indsB);

for i = 1 : size(IndicesB,1)
    
    for rep = 1 : maxRep
        if DB(IndicesB(i,1),IndicesB(i,2),IndicesB(i,3),IndicesB(i,4),IndicesB(i,5),rep) == 0
           DB(IndicesB(i,1),IndicesB(i,2),IndicesB(i,3),IndicesB(i,4),IndicesB(i,5),rep) = i;
           Test(IndicesB(i,1),IndicesB(i,2),IndicesB(i,3),IndicesB(i,4),IndicesB(i,5),rep) = 1;
           break;
        else
            continue
        end
    end
end

for i = 1 : maxRep
    
    blah(i) = length(find(Test(:,:,:,:,:,i)));
    
end

percentRepeat = 100* (sum(blah(2:end)) / sum(blah(1:end)));

% *******************
indsMask = find(DB);

mInds = DB(indsMask);
kSpace = DatB(mInds);
W = Weights(mInds);

product = ArraySize(1) * ArraySize(3) * ArraySize(4) * ArraySize(2) * ArraySize(6);
indsMask = mod(indsMask, tmp);
indsMask(indsMask==0) = product;


% % D = zeros(1, 1, ArraySize(3), ArraySize(4), ArraySize(5), ArraySize(6), maxRep,'single');    
% for i = 1 : size(Dat, 3)
%     
%     for rep = 1 : maxRep
%         if D(1, 1, Lin(1,1,i), Par(1,1,i), Set(1,1,i), cBins(1,1,i), rep) == 0
%            D(1, 1, Lin(1,1,i), Par(1,1,i), Set(1,1,i), cBins(1,1,i), rep) = i;
%            break;
%         else
%             continue
%         end
%     end
% end
% 
% indsMask = find(D);

% mInds = D(indsMask);
% kSpace = Dat(:,:,mInds);
% weights = rWeights(:,:,mInds);
% 
% 
% tmp = ArraySize(3) * ArraySize(4) * ArraySize(5) * ArraySize(6);
% 
% indsMask = mod(indsMask, tmp);
% indsMask(inds==0) = tmp;


PlaceHolder = [];



