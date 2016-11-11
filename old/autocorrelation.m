% autocorrelation.m
% calculates the normalised autocorrelation of a feature frame
% 

% Author: Max Fisher

function ac = autocorrelation(feature_frame)
	frame_length = length(feature_frame);

	% normalise to mean zero and power 1
	mean = sum(feature_frame)/frame_length;
	power = sum(feature_frame.^2)/frame_length;
	feature_frame = (feature_frame - mean)/sqrt(power);

	ac = xcorr(feature_frame, 'none');
	% shift so that zero lag is at the first index
	% make sure to shift the right dimension!!
	% -> shift the acf along the dimension of the feature
	% frame that isn't 1
	% since the frame is half the length of the acf,
	% positive or negative shifts don't matter
	shift_array = size(feature_frame);
	shift_array(shift_array == 1) = 0;
	ac = circshift(ac, shift_array);
	% discard redundant second half of samples and correct for
	% bias caused by shifting (i.e. unbiased autocorrelation)
	ac = ac(1:frame_length)./(frame_length:-1:1)';

	% normalise so that first coefficient is 1, that way we can compare
	% between different autocorrelations
	ac = ac/ac(1);
end
