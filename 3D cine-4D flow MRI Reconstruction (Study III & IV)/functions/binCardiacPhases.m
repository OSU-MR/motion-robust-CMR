function [cBins, meanRR, stdRR, triggerTime, triggers] = binCardiacPhases(lengthDat,C,opt)

% % Estimate cardiac period from frequency spectrum
% C_freq = fftshift(fft(C));
% df = [-length(C)/2:1:length(C)/2-1]*(opt.fs / length(C));
% C_freq_mag = abs(C_freq);
% cardiac_freq = max(df(find(C_freq_mag == max(C_freq_mag))));
% cardiac_per = 1 / cardiac_freq;
% 
% % Find typical bin width
% binWidth = cardiac_per / opt.nCPhases;

% Typical length of systole ~ 1/3 * HR 
% (Find some way to make this HR dependant) ???????????????????????????????
nPhases_systole = ceil(opt.nCPhases/3);
% nPhases_diastole = opt.nCPhases - nPhases_systole;

% Ensure cardiac signal is oriented correctly
% (Output estimate of dias/sys duration) ??????????????????????????????????
C = orientCardiac_2( C, opt);



% C = resample(C,2,1);
% interp_factor = opt.sgI;
interp_factor = 2;
if interp_factor ~= 1
    xq = (1/interp_factor):(1/interp_factor):length(C);
    C = interp1(C, xq, 'spline');
end
% % Flipped first derivative of cardiac signal
% dC = -gradient( C );
% % Threshold out values less than 0
% dC(dC < 0) = 0;
% % Nonlinear scaling to eliminate smaller peaks
% dC = dC.^2;
% dC(dC<0.15*max(dC)) = 0;
% dC = sqrt(dC);

% *************************************************
% C = resample(C,2,1);
% Flipped first derivative of cardiac signal
dC = -gradient(C);
% dC = gradient(C); % !!!!! Change back !!!!!!!
% C = -C;
% Threshold values less than 0
dC(dC < 0) = 0;

% % find envelopes
% wl = round(1*interp_factor/opt.TR);
% env_1 = envelope(dC, wl, 'peak');
% cutOff_1 = env_1 - 0.25*prctile(dC,95);
% 
% [pks, loc] = findpeaks(dC);
% indx = find(pks > cutOff_1(loc));
% loc(indx) = []; pks(indx) = [];
% 
% env_2 = interp1(loc, pks, 1:1:length(dC), 'spline');
% 
% cutOff = mean(cat(1,env_1,env_2), 1);
% plot(dC); hold on; plot(cutOff);

% wl = round(20*interp_factor/opt.TR);
% kernal = ones(1,wl) * (1/wl);
% dC_sq = dC.^2;
% dC_rms = sqrt(conv(dC_sq, kernal, 'same'));
% cutOff = dC_rms;
% figure;
% plot(dC); hold on; plot(cutOff);


% Find peaks initially
[pks, loc] = findpeaks(dC);
% cutOff = median(pks);
cutOff = prctile(dC, 95)*0.5;
figure;
plot(dC); hold on;
plot([1,length(dC)],[cutOff,cutOff],'--k');

dC(dC < cutOff) = 0;



% indx = kmeans(pks, 2);
% tmp = [mean(pks(indx==1)),mean(pks(indx==2))];
% maxGroup = find(tmp == max(tmp));
% minGroup = find(tmp == min(tmp));

% cutOff = mean([max(pks(indx==minGroup)), min(pks(indx==maxGroup))]);

% %****************
% plot(dC); hold on;
% plot([1,length(dC)],[max(pks(indx==minGroup)),max(pks(indx==minGroup))],'--r');
% plot([1,length(dC)],[min(pks(indx==maxGroup)),min(pks(indx==maxGroup))],'--r');
% plot([1,length(dC)],[cutOff,cutOff],'--k');
% 
% % ***************
% dC(dC < cutOff) = 0;


% Find peaks of gradient - corresponds to peak systole/ejection
[pks, loc] = findpeaks(dC);
nPeaks = length(pks);

[pks_2, loc_2] = findpeaks(C);

dC_tmp = -gradient(C);
% Trim first peak for robustness
loc = loc(2:end);
loc_temp = zeros(length(loc),1);
loc_true = zeros(length(loc),1);
for i = 1 : length(loc)
    counter = 0;
    stop = 0;
    while ~stop
        counter = counter + 1;
        if (loc(i) - counter) < 1
            stop = 1;
        else
            peak(1) = dC_tmp(loc(i) - counter + 1);
            peak(2) = dC_tmp(loc(i) - counter);
        
%             if peak(2) <= 0
%                 % Check sample is closest to zero (peak)
%                 if abs(dC_tmp(loc(i) - counter)) <= abs(dC_tmp(loc(i) - counter + 1))
%                     loc_temp(i) = loc(i) - counter;
%                     m = dC_tmp(loc(i)-counter+1)-dC_tmp(loc(i)-counter);
%                     b = dC_tmp(loc(i)-counter) - m*(loc(i)-counter);
%                     loc_true(i) = -b/m;
%                     
%                 else
%                     loc_temp(i) = loc(i) - counter + 1;
%                     m = dC_tmp(loc(i)-counter+1)-dC_tmp(loc(i)-counter);
%                     b = dC_tmp(loc(i)-counter) - m*(loc(i)-counter);
%                     loc_true(i) = -b/m;
%                 end
%                 clear peak
%                 stop = 1;
%             end
            if counter >= 3
                if counter == 3
                    previous_location = loc(i)-counter;
                    previous_dif = inf;
                else
                    previous_location = current_location;
                    previous_dif = current_dif;
                end
                
                % Linear regression
                m = dC_tmp(loc(i)-counter+1)-dC_tmp(loc(i)-counter);
                b = dC_tmp(loc(i)-counter) - m*(loc(i)-counter);
                loc_temp(i) = loc(i) - counter + 1;
                current_location = -b/m;
                current_dif = current_location - previous_location;
                
                if abs(current_dif) >= abs(previous_dif)
                    loc_true(i) = previous_location;
                    clear peak;
                    stop = 1;
                end
            end
        end
    end
end

loc_temp(loc_temp == 0) = [];

% for i = 1 : length(loc)
%     if loc(i) > loc_2(1)
%         loc(1:i-1) = [];
%         break;
%     end
% end
% 
% for i = 1 : length(loc_2)
%     if loc_2(i) > loc(1)
%         index = i;
%         break;
%     end
% end
% 
% if index > 2
%     loc_2(1:index - 2) = [];
% end
% 
% gradLoc = gradient(loc);
% trimIndx = find(abs((gradLoc - mean(gradLoc))/std(gradLoc)) > 3);
% loc(trimIndx) = [];
% 
% loc_temp = zeros(length(loc), 1);
% % ind = zeros(length(loc));
% for i = 1 : length(loc)
%     if i == 1
%         dif = loc(1) - loc_2(1);
%         loc_temp(1) = loc_2(1);
%         loc = loc - dif;
% %         ind(1) = 1;
%     else
%         Dif = abs(loc(i)-loc_2);
%         indx = find(Dif == min(Dif));
%         if indx > 1
%             loc_2(i:indx-1) = [];
%         end
%         dif = loc(i) - loc_2(i);
%         loc_temp(i) = loc_2(i);
%         loc = loc - dif;
% %         ind(i) = 1;
%     end
% end

loc = loc_temp;
nPeaks = length(loc);
nPeaks = length(loc_true);


% Time
timeSG = [1:1:length(C)]*opt.TR*(1/interp_factor);
timePeaks = timeSG(loc);
timePeaks_true = interp1([1:length(timeSG)], timeSG, loc_true);
timePeaks = timePeaks_true;
timePeaks(isnan(timePeaks)) = [];
nPeaks = length(timePeaks);

tmp = [0:1:length(C)-1]*opt.TR*(1/interp_factor);
tmp_true = interp1([1:length(tmp)], tmp, loc_true);
triggers = tmp_true;
% % TEST to check systolic/diastolic intervals wrt RR interval
% % **********************************************************
% C_tmp = -C; C_tmp(C_tmp<0) = 0;
% [pks_tmp, loc_tmp] = findpeaks(C_tmp);
% nPeaks_tmp = length(loc_tmp);
% timePeaks_tmp = timeSG(loc_tmp);
% timePeaks_tmp(1) = []; loc_tmp(1) = [];
% % timePeaks_tmp = timePeaks_tmp(1:end-1);
% 
% % Find variability statistics for arrhythmia rejection
% nPeaks = length(timePeaks);
% RR = zeros(1,nPeaks-1);
% for i = 1 : nPeaks-1
%     RR(i) = timePeaks(i+1) - timePeaks(i); 
% end
% RR(125) = [];
% timePeaks_tmp(126) = [];
% timePeaks(end) = [];
% sys = timePeaks_tmp - timePeaks; %sys = sys(1:end-1);
% dia = RR - sys; 
% % RR(124) = [];
% % sys(124:125) = [];
% % dia(124:125) = [];
% % 
% % RR(end) = [];
% % sys(1) = [];
% % dia(1) = [];
% 
% z = abs(RR-mean(RR))/std(RR);
% inds_z = find(z > 1.5);
% RR(inds_z) = [];
% sys(inds_z) = [];
% dia(inds_z) = [];
% % **********************************************************

% Find variability statistics for arrhythmia rejection
RR = zeros(1,nPeaks-1);
for i = 1 : nPeaks-1
    RR(i) = timePeaks(i+1) - timePeaks(i); 
end


% Remove outliers
tmp_RR = RR;
tmp_RR(abs((tmp_RR-mean(tmp_RR))/std(tmp_RR)) > opt.arrhyRej) = [];


medRR = median(tmp_RR); meanRR = mean(tmp_RR); stdRR = std(tmp_RR);
percentile_95 = prctile(tmp_RR, 95);
percentile_25 = prctile(tmp_RR, 25);

binWidth = meanRR/opt.nCPhases;
% binWidth = percentile_25/opt.nCPhases;
nPhases_systole = opt.nCPhases;
% nPhases_diastole = opt.nCPhases - nPhases_systole;


% Define which are systolic bins:
% sBefore = ceil(nPhases_systole/2);
% sAfter = nPhases_systole - sBefore;
binsS = zeros(1,opt.nCPhases);
% binsS(1:sBefore) = 1; binsS(end-sAfter+1:end) = 1; 
binsS(1:nPhases_systole) = 1;

% *********************************************************************
% % Test of Coherent Point Drift algorithm for cardiac cycle registration
% timeSG = timeSG';
% clear cycle
% for i = 1 : length(loc)-1
%     clear tmp
%     tmp(:,1) = timeSG(loc(i):loc(i+1)-1) - timeSG(loc(i)); 
%     tmp(:,2) = C(loc(i):loc(i+1)-1); tmp(:,2) = (tmp(:,2)-min(tmp(:,2)))./(max(tmp(:,2)-min(tmp(:,2))));
%     tmp(:,3) = zeros(size(tmp,1),1);
%     cycle{i} = tmp;
% end
% 
% for i = 2 : length(cycle)
%     clear tmp tmp2 xq vq1
%     tmp = cycle{i};
%     xq = linspace(0,tmp(end,1),100);
%     vq1 = interp1(tmp(:,1),tmp(:,2),xq);
%     tmp2(:,1) = xq;
%     tmp2(:,2) = vq1;
%     tmp2(:,3) = zeros(100,1);
%     cycle2{i} = tmp2;
% end
%     
% meanT = zeros(100,1);
% meanC = zeros(100,1);
% for i = 2 : length(cycle2)
%     meanT = (meanT + cycle2{i}(:,1));
%     meanC = (meanC + cycle2{i}(:,2));
% end
% meanT = meanT/308;
% meanC = meanC/308;
% 
% 
% fixedCloud = pointCloud([meanT, meanC, zeros(100,1)]);
% for i = 2 : length(cycle2)
%     movingCloud = pointCloud(cycle2{i});
%     tform = pcregistercpd(movingCloud,fixedCloud);
%     movingReg = pctransform(movingCloud,tform);    
%     cycle3{i} = movingReg.Location;   
% end
% 
% figure;
% for i = 2 : length(cycle2)
% %     plot(cycle3{i}(:,1), cycle3{i}(:,2)); hold on;
% %     plot(cycle2{i}(:,2)); hold on;
%     plot(cycle3{i}(:,1)); hold on;
% end
% 
% for i = 2 : length(cycle3)
%     for j = 1 : 99
%         FD_b(j,i) = cycle2{i}(j+1,1) - cycle2{i}(j,1);
%         FD_a(j,i) = cycle3{i}(j+1,1) - cycle3{i}(j,1);
%     end
%     dif(:,i) = cycle3{i}(:,1) - cycle2{i}(:,1);
% end
% ratio = FD_a./FD_b;
% 
% ratio = ratio(:,2:end);
% RR_new = RR(2:end);
% 
% [RR_sort, sortInds] = sort(RR_new, 'ascend');
% for i = 1 : size(ratio,1)
%     ratioTmp(i,:) = ratio(i,sortInds);
% end
% ratioNew = ratioTmp;
% p1 = ratioNew(20,:);
% p2 = ratioNew(40,:);
% p3 = ratioNew(60,:);
% p4 = ratioNew(80,:);
% 
%     
% 
% 
% figure
% pcshowpair(movingCloud,fixedCloud,'MarkerSize',50)
% xlabel('X')
% ylabel('Y')
% zlabel('Z')
% title('Point clouds before registration')
% legend({'Moving point cloud','Fixed point cloud'},'TextColor','w')
% legend('Location','southoutside')
% 
% tform = pcregistercpd(movingCloud,fixedCloud);
% movingReg = pctransform(movingCloud,tform);
% 
% figure
% pcshowpair(movingReg,fixedCloud,'MarkerSize',50)
% xlabel('X')
% ylabel('Y')
% zlabel('Z')
% title('Point clouds after registration')
% legend({'Moving point cloud','Fixed point cloud'},'TextColor','w')
% legend('Location','southoutside')
% *********************************************************************

% Create bins.
counter = 0;
for i = 1 : nPeaks-1
    dT = timePeaks(i+1) - timePeaks(i);
%     dScaleWidth = (dT - (nPhases_systole*binWidth)) / nPhases_diastole;
%     diastoleWidth(i) = dScaleWidth;
    for n = 1 : opt.nCPhases
        counter = counter + 1;
        % Arrhythmia rejection (simple Z-score threshold)
        if abs((dT-meanRR) / stdRR) > opt.arrhyRej
            bin(counter, 1) = NaN;
            bin(counter, 2) = NaN;
        else
            if n == 1
                binStart = timePeaks(i) - (binWidth/2) - binWidth;
                bin(counter, 1) = binStart;
                bin(counter, 2) = binStart + binWidth;
            elseif n > 1
                if binsS(n)
                    bin(counter, 1) = bin(counter-1, 2);
                    bin(counter, 2) = min(bin(counter-1, 2) + binWidth, timePeaks(i+1) - (binWidth/2) - binWidth);
                else
                    bin(counter, 1) = bin(counter-1, 2);
                    bin(counter, 2) = bin(counter-1, 2) + dScaleWidth;
                end
            end
        end
    end
end

if opt.viewSharing
    % Create bins offset bins
    counter = 0;
    for i = 1 : nPeaks-1
        dT = timePeaks(i+1) - timePeaks(i);
        %     dScaleWidth = (dT - (nPhases_systole*binWidth)) / nPhases_diastole;
        %     diastoleWidth(i) = dScaleWidth;
        for n = 1 : opt.nCPhases
            counter = counter + 1;
            % Arrhythmia rejection (simple Z-score threshold)
            if abs((dT-meanRR) / stdRR) > opt.arrhyRej
                binV(counter, 1) = NaN;
                binV(counter, 2) = NaN;
            else
                if n == 1
                    binStart = timePeaks(i) - (binWidth/2) - binWidth + (binWidth/2);
                    binV(counter, 1) = binStart;
                    binV(counter, 2) = binStart + binWidth;
                elseif n > 1
                    if binsS(n)
                        binV(counter, 1) = binV(counter-1, 2);
                        binV(counter, 2) = min(binV(counter-1, 2) + binWidth, timePeaks(i+1) - (binWidth/2) - binWidth);
                    else
                        binV(counter, 1) = binV(counter-1, 2);
                        binV(counter, 2) = binV(counter-1, 2) + dScaleWidth;
                    end
                end
            end
        end
    end
    
    tmp(:,:,1) = bin;
    tmp(:,:,2) = binV;
    
    bin = tmp;
    
%     bin_tmp = zeros(2*size(bin,1),2);
%     num = 1:size(bin,1);
%     bin_tmp(2*num-1,:) = bin(num,:);
%     bin_tmp(2*num,:) = binV(num,:);
%     
%     bin = bin_tmp;
%     
    opt.nCPhases = 2*opt.nCPhases;
    
end

% =========================================================================
% plot bins
C_triggers = interp1(timeSG, C, timePeaks);
figure;
for j = 1 : length(bin)
    plot([bin(j,1,1),bin(j,1,1)],[min(C),max(C)],'--k'); hold on; 
end
plot(timeSG,C,'b'); hold on; plot(timePeaks, C_triggers,'.r','MarkerSize',12); hold off;
% =========================================================================


% Now we need to divide the k-space into the bins
% Timestamp for each line
if ~strcmp(opt.sigExt, 'PTSola')
    timeDat = [1:1:lengthDat]*(opt.TR/opt.sgI);
else
    timeDat = [1:1:lengthDat]*(opt.TR);
end

% Trim out SG signal
% timeDat(:,opt.sgI:opt.sgI:end) = [];
% Divide k-space lines into the cardiac bins
DatBin = zeros(1,length(timeDat));

for i = 1 : length(timeDat)
    for j = 1 : length(bin)
        if j == 1
            binStart = timePeaks(1) - (binWidth/2);
%             if timeDat(i) < bin(j, 2) && timeDat(i) >= binStart
            if timeDat(i) < bin(j, 2) && timeDat(i) >= bin(j, 1)
                DatBin(i) = j;
            end
        elseif j > 1
%             if timeDat(i) < bin(j) && timeDat(i) >= bin(j-1)
            if timeDat(i) < bin(j, 2) && timeDat(i) >= bin(j, 1)
                DatBin(i) = j;
            end
        end
    end
end

% Regions to trim later
DatBin(DatBin ==0) = NaN;
% offset 1st bin into late diastole
DatBin = DatBin + opt.binOffset;

% modulo operator by number of cardic phases
DatBin_Phases = mod(DatBin,opt.nCPhases);
DatBin_Phases(DatBin_Phases == 0) = opt.nCPhases;

% Output
cBins = DatBin_Phases;

triggerTime = [];
% triggerTime = zeros(1,opt.nCPhases);
% triggerTime(1:nPhases_systole) = [0:1:nPhases_systole-1]*binWidth;
% triggerTime(nPhases_systole+1:end) = [1:nPhases_diastole]*mean(diastoleWidth) + triggerTime(nPhases_systole);


end

