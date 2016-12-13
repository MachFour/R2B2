% autocorrelation.m
% calculates the autocorrelation of a feature frame, optionally using an
% equal-length frame of previous samples to slide it back against.
% RMS power is normalised to 1 (but mean is not subtracted)
% Both input should be column vectors

% Author: Max Fisher

function acf = autocorrelation(data_frame, prev_samples)
	data_len = length(data_frame);
	if data_len ~= size(data_frame, 1)
		warning('data frame is probably not a column vector!');
	end

	use_prev_samples = 0;

	if nargin == 2
		len_prev = length(prev_samples);
		if data_len == len_prev
			use_prev_samples = 1;
		else
			warning('previous samples are not of same length as data; ignoring.');
		end
	end

	% since this is constant, we can normalise its power outside the loop
	data_frame_power = sqrt(sum(data_frame.^2)/data_len);
	normalised_data_frame = data_frame/data_frame_power;

	acf = zeros(data_len, 1);

	% compound_frame is what we'll take a shifting window of, to compare with the 
	% data frame. When the shift is zero, the data frame is just prepended with zeros.

	if use_prev_samples
		compound_frame = [prev_samples; data_frame];
	else
		% replace prev samples with zeros
		compound_frame = [zeros(data_len, 1); data_frame];
	end

	for shift = 0:data_len-1
		shifted_frame = compound_frame(data_len - shift + 1: end - shift);
		shifted_frame_power = sqrt(sum(shifted_frame.^2)/data_len);
		normalised_shifted_frame = shifted_frame/shifted_frame_power;

		acf(shift+1) = sum(normalised_data_frame.*normalised_shifted_frame)/data_len;
	end

	% normalise so that first coefficient is 1, that way we can compare
	% between different autocorrelations
	% -> the RMS power normalisation in the previous step makes this
	% unnecessary.
end
