% autocorrelation.m
% calculates the autocorrelation of a feature frame
% by splitting it in half, sliding the second half back across the first
% half, and computing the inner product of the two.
% The frame is first normalised to have unit power.

% Author: Max Fisher

function acf = autocorrelation(feature_frame)
	len = length(feature_frame);
	if mod(len, 2) ~= 0
		warning('Autocorrelation frame is an odd length, can''t split in half');
	end

	half_len = len/2;

	% normalise RMS power to 1 (not mean zero)
	% -> but do normalisation inside the loop, since every half_len slice of the
	% feature frame has a different power

	% this is what we'll shift back and compare against the whole frame.
	% when the shift is zero, the second half of the frame is compared against itself.
	half_frame = feature_frame(half_len + 1: end);
	% since this is constant, we can normalise its power outside the loop
	half_frame_power = sqrt(sum(half_frame.^2)/half_len);
	normalised_half_frame = half_frame/half_frame_power;

	acf = zeros(size(half_frame));

	for shift = 0:half_len-1
		shifted_frame = feature_frame(half_len - shift + 1: end - shift);
		shifted_frame_power = sqrt(sum(shifted_frame.^2)/half_len);
		normalised_shifted_frame = shifted_frame/shifted_frame_power;

		acf(shift+1) = sum(normalised_half_frame.*normalised_shifted_frame)/half_len;
	end

	% normalise so that first coefficient is 1, that way we can compare
	% between different autocorrelations
	% -> the RMS power normalisation in the previous step makes this unnecessary
	%acf = acf/acf(1);
end
