% helper function to find peaks of autocorrelation, while taking into account
% the harmonics inherent in autocorrelation

% Author: Max Fisher

function comb = autocorrelation_comb(length, lag_samples)
% make impulse matrix to estimate 

%prior tempo weighting function (rayleigh distribution)
% tempo lag in samples
beta = 60/120; % prior has maximum at 120BPM
tempo_weighting = @(tempo_lag) tempo_lag/beta^2*exp(-tempo_lag^2/(2*beta^2));
comb = zeros(length, 1);
NUM_COMBS = 4;
COMB_BASE_WIDTH = 2;
for p = 1:NUM_COMBS
    pulse_min_index = round((lag_samples*p - COMB_BASE_WIDTH*p + 1));
    pulse_max_index = round((lag_samples*p + COMB_BASE_WIDTH*p - 1));
    comb(pulse_min_index:pulse_max_index) = 1/(NUM_COMBS*(1+pulse_max_index - pulse_min_index));
    % make sure it's the right length;
    comb = comb(1:length);
end;

% put a peak of width 2*p-1, height 1/(2*p-1), centred 
% at the (tempo_lag*p)th sample, for p=1:4;
% tempo lag varies from 1:FEATURE_WIN_LENGTH/4 samples

%stem(comb);
end
