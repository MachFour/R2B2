% beat_alignment_function.m
% follows method in 'Context Dependent Beat tracking' to find the location of beats
% in a frame of an onset detection function.

% tempo spacing is given in samples.

% frame should be a column vector. Returns a column vector

% Author: Max Fisher
function baf = beat_alignment_function_3impulse(frame, tempo_spacing)
	SPIKE_WIDTH = 2;
	NUM_SPIKES = 3;

	frame_length = size(frame, 1);
	if frame_length ~= length(frame)
		warning('Frame is likely not a column vector!');
	end

	%normalise frame to power 1
	frame_power = sqrt(sum(frame.^2)/frame_length);
	frame = frame./frame_power;

	% we reverse frame to examine most recent onsets first
	frame = flipud(frame);

	impulse_train = zeros(size(frame));
	% tempo is in lag_units (seconds)
	for spike_idx = 1:NUM_SPIKES
		% double() for matlab type compatibility
		indices = spike_idx*double(tempo_spacing) + (1:SPIKE_WIDTH);
		if max(indices) > frame_length
			%warning('spikes do not fit inside frame');
		else
			impulse_train(indices) = 1;
		end

	end

	% normalise impulse train to have sum 1,
	% so it has the effect of a weighted sum
	% of odf samples
	% note that if spikes fall off the end, this won't be accounted for
	impulse_train = impulse_train/(SPIKE_WIDTH*NUM_SPIKES);

%  	figure; plot(impulse_train); hold on; plot(frame);
%  	title(sprintf('Reversed feature frame and impulses for comb width %d', tempo_spacing));
%  	xlabel('Samples');

	% only makes sense to consider shifts up to tempo_spacing, since
	% beyond that point the comb 'teeth' would fall off the end

	baf = zeros(tempo_spacing, 1);
	for shift = 0:tempo_spacing-1
		shifted_frame = circshift(frame, -shift);
		inner_product = shifted_frame.*impulse_train;
		baf(shift+1) = sum(inner_product);
	end

	%normalise so that all possible beat alignments sum to 1: probability
	%baf = baf/sum(baf);

% 	figure; plot(baf);
% 	title(sprintf('Beat alignment function - tempo spacing = %d', tempo_spacing));
% 	xlabel('Comb shift (samples)');
% 	ylabel('Strength');

end
