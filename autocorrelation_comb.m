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
NUM_COMBS = 3;

% tempo lag varies from 1:FEATURE_WIN_LENGTH/4 samples
for p = 1:NUM_COMBS
	comb(1+lag_samples*p) = 2;
	comb(1+lag_samples*p - 1) = 1;
	comb(1+lag_samples*p + 1) = 1;
end

% normalise to have sum 1
comb = comb/sum(comb);


%stem(comb);
end
