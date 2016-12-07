% beat_alignment_function.m
% follows method in 'Context Dependent Beat tracking' to find the location of beats
% in a frame of an onset detection function.

% tempo spacing is given in samples.

% frame should be a column vector. Returns a column vector

% Author: Max Fisher
function baf = beat_alignment_function(frame, tempo_spacing)
	frame_length = size(frame, 1);
	if frame_length ~= length(frame)
		warning('Frame is likely not a column vector!');
	end

	% make mean zero but keep variance?
	% have to think about this
	%frame = frame - mean(frame);

	% we reverse frame to examine most recent onsets first
	frame = flipud(frame);
	SPIKE_WIDTH = 2;

	impulse_train = zeros(size(frame));
	% tempo is in lag_units (seconds)
	for w = 1:SPIKE_WIDTH
	    impulse_train(w:tempo_spacing:end-tempo_spacing) = 1;
	end
	% make impulses of decreasing height to weight more recent onsets

	% exponential decay of weighting:
	% w(n) = a^(-a*n/N)
	%choose a = 2
	a = 2;
	impulse_weighting = a.^(-a/frame_length*(1:frame_length)');
	impulse_train = impulse_train.*impulse_weighting;

	% normalise impulse train to have sum 1,
	% so it has the effect of a weighted sum
	% of odf samples
	%impulse_train = impulse_train/sum(impulse_train);

%  	figure; plot(impulse_train); hold on; plot(frame);
%  	title(sprintf('Reversed feature frame and impulses for comb width %d', tempo_spacing));
%  	xlabel('Samples');

	% only makes sense to consider shifts up to tempo_spacing, since
	% beyond that point the comb 'teeth' would fall off the end

	baf = zeros(tempo_spacing, 1);
	for shift = 0:tempo_spacing-1
		shifted_frame = circshift(frame, -shift);
		baf(shift+1) = sum(shifted_frame.*impulse_train);
	end;

	%normalise so that all possible beat alignments sum to 1: probability
	baf = baf/sum(baf);

% 	figure; plot(baf);
% 	title(sprintf('Beat alignment function - tempo spacing = %d', tempo_spacing));
% 	xlabel('Comb shift (samples)');
% 	ylabel('Strength');

end
