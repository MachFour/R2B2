% helper function to find peaks of autocorrelation, while taking into account
% the harmonics inherent in autocorrelation

% Author: Max Fisher

function comb = autocorrelation_comb(len, lag_samples)
% make impulse matrix to estimate 

%prior tempo weighting function (rayleigh distribution)
% tempo lag in samples
beta = 60/120; % prior has maximum at 120BPM
tempo_weighting = @(tempo_lag) tempo_lag/beta^2*exp(-tempo_lag^2/(2*beta^2));
comb = zeros(len, 1);
num_combs = floor((len-3)/lag_samples);

if num_combs < 2
	warning('less than 2 combs can fit');
end

for p = 1:num_combs
	comb(1+lag_samples*p -2) = 1;
	comb(1+lag_samples*p -1) = 4;
	comb(1+lag_samples*p +0) = 9;
	comb(1+lag_samples*p +1) = 4;
	comb(1+lag_samples*p +2) = 1;
end

% ensure comb has correct length (should not be needed)

if length(comb) > len
	warning('correcting comb length from %d to %d', ...
		length(comb), len);
	comb = comb(1:len);
end

% normalise to have sum 1
comb = comb/sum(comb);


%stem(comb);
end
