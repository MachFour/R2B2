% Creates a filter bank (in the frequency domain) of triangular filters which are
% spaced linearly in the Mel frequency scale, which is designed to mirror the
% way the ear responds to frequency in a logarithmic way
% The matrix returned by this function can be multiplied by the power spectrum of a signal
% (with given Nfft) in order to get the energy output of each filter (in its frequency range)

% The filter bank matrix output by this function has
% num_filters rows and Nfft/2 columns (since the right half
% of the DFT magnitude spectrum is assumed to be discarded)

% Code modified from http://www.cs.tut.fi/kurssit/SGN-4010/Ex2013/Create_MelFrequencyFilterBank.m
% Author: Max Fisher

function mel_filterbank = mel_filter(fs, Nfft, num_filters)
Hz_to_mel = @(f) 2595*log10(1+f/700);
mel_to_Hz = @(m) 700*(10.^(m/2595)-1);
% filter centre frequencies, linearly spaced in the mel scale from 0 to
% fs/2
filter_centre_mels = (0:(num_filters+1))/(num_filters+1) * Hz_to_mel(fs/2);
% filter centre frequencies, mapped to Hz scale
filter_centre_freqs = mel_to_Hz(filter_centre_mels);

% quantize mel centre freqencies (in Hz) to DFT bin indexes
filter_centre_bins = round(filter_centre_freqs/(fs/2)*(Nfft/2));

% Just correcting the first bin not to be 0.
filter_centre_bins(1) = 1;

if filter_centre_bins(2) == 0,
	error('Not enough DFT bins for this number of filters!')
end

%disp(filter_centre_bins);

% Each row of 'mel_filterbank' corresponds to one triangular filter.
mel_filterbank = zeros(num_filters, Nfft/2);

for filter_index = 1:num_filters
	% bins over which this filter has nonzero magnitude response
	low_bin = filter_centre_bins(filter_index);
	central_bin = filter_centre_bins(filter_index+1);
	high_bin = filter_centre_bins(filter_index+2);

	% set the magnitude response for the two edges of the triangle
	% note increase/decrease is linear in the Hz scale
	mel_filterbank(filter_index, low_bin:central_bin) = ...
		(0:(central_bin - low_bin))/(central_bin - low_bin);
	mel_filterbank(filter_index, central_bin:high_bin) = ...
		((high_bin - central_bin):-1:0)/(high_bin - central_bin);
end

% example plotting code (from source)
if ~nargout
figure('Name','Mel Frequency Filter Bank')
problem4fig = gcf;
set(problem4fig,'position',[1 700 560 420]);
plot(mel_filterbank')
axis([0, Nfft/2,  0,  1.1])
title('Mel-filterbank')
end
