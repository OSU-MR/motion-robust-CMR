function [ C ] = orientCardiac_2( C0, opt )

disp('Determining cardiac orientation...');

% =========================================================================
% Test 1 - Peak gradients
% Assumption is that rate of change of ventricular contraction is greater
% than relaxation
flipC.test1 = 0;

% Find first derivative of C
dC0 = gradient(C0);

% Threshold derivative so that smaller peaks are eliminated
% range_dC0 = prctile(dC0,95) - prctile(dC0,5);
% dC0(abs(dC0) < 0.3*range_dC0) = 0; % Pick better threshold later???????????

% Find peaks
% Negative
dC0_N = -dC0;
dC0_N(dC0_N<0) = 0;
% dC0_N = dC0_N .^2;
% dC0_N(dC0_N<0.15*prctile(dC0_N,99)) = 0;
% dC0_N = sqrt(dC0_N);

% cutoff = prctile(dC0_N, 95)*0.4;

wl = round(10*(1/opt.TR));
kernal = ones(1,wl) * (1/wl);
dC0_N_sq = dC0_N.^2;
dC0_N_rms = sqrt(conv(dC0_N_sq, kernal, 'same'));
cutoff = dC0_N_rms;

dC0_N(dC0_N<cutoff) = 0;
[pk_dCn, loc_dCn] = findpeaks(dC0_N);

% Positive
dC0_P = dC0;
dC0_P(dC0_P<0) = 0;
% dC0_P = dC0_P .^2;
% dC0_P(dC0_P<0.15*prctile(dC0_P,99)) = 0;
% dC0_P = sqrt(dC0_P);

% cutoff = prctile(dC0_P, 95)*0.4;

wl = round(10*(1/opt.TR));
kernal = ones(1,wl) * (1/wl);
dC0_P_sq = dC0_P.^2;
dC0_P_rms = sqrt(conv(dC0_P_sq, kernal, 'same'));
cutoff = dC0_P_rms;

dC0_P(dC0_P<cutoff) = 0;
[pk_dCp, loc_dCp] = findpeaks(dC0_P);

% Results from test 1
if mean(abs(pk_dCn)) >= mean(abs(pk_dCp))
    flipC.test1 = 0;
elseif mean(abs(pk_dCn)) < mean(abs(pk_dCp))
    flipC.test1 = 1;
end

% =========================================================================
% Test 2 - Compare relative distances b/w peaks of derivative
% Assumption that the duration of diastole is greater than the duration of
% systole

% Already have the positive and negative peaks from test 1...
% Distance b/w postive peak and next negative peak
loc_dCn_shift = loc_dCn;
counter = 0;
while loc_dCp(1) >= loc_dCn_shift(1)
    loc_dCn_shift = circshift(loc_dCn_shift,-1);
    counter = counter + 1;
end
tmp = min(length(loc_dCp),length(loc_dCn_shift));
a = abs(mean(loc_dCp(1:tmp-counter) - loc_dCn_shift(1:tmp-counter)));

% Distance b/w negative peak and next positive peak
loc_dCp_shift = loc_dCp;
counter = 0;
while loc_dCn(1) >= loc_dCp_shift(1)
    loc_dCp_shift = circshift(loc_dCp_shift,-1);
    counter = counter + 1;
end
tmp = min(length(loc_dCn),length(loc_dCp_shift));
b = abs(mean(loc_dCn(1:tmp-counter) - loc_dCp_shift(1:tmp-counter)));

% Results from test 2
if b < a
    flipC.test2 = 0;
elseif a <= b
    flipC.test2 = 1;
end

% =========================================================================
% Finally, compare the two tests
if ~flipC.test1 && ~flipC.test2
    % Keep same orientation
    C = C0;
    disp('Correct orientation detected.');
elseif flipC.test1 && flipC.test2
    % Flip orientation
    C = -C0;
    disp('Incorrect orientation detected...flipping signal');
else
    % Default to result from test 1
    if ~flipC.test1
        disp('Conflicting tests...keeping same orientation');
        C = C0;
    else
        disp('Conflicting tests...flipping signal');
        C = -C0;
    end
end


