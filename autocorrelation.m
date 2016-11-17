% autocorrelation.m
% calculates the autocorrelation of a feature frame
% by splitting it in half, sliding the second half back across the first
% half, and computing the inner product of the two.
% The frame is first normalised to have unit power.

% Author: Max Fisher

function acf = autocorrelation(feature_frame)
	len = length(feature_frame);
	if mod(len, 2) ~= 0
		warning('Autocorrelation frame cannot be split in half evenly');
	end

	% normalise RMS power to 1 (not mean zero)
	power = sqrt(sum(feature_frame.^2)/len);
	feature_frame = feature_frame./power;

	second_half = feature_frame(len/2+1:end);
	size(second_half)
	disp(len/2);

	acf = zeros(size(second_half));

	for shift = 0:len/2-1
		shifted_frame = feature_frame((len/2 + 1 - shift):end - shift);
		acf(shift+1) = sum(second_half.*shifted_frame)/(len/2);
	end

	% normalise so that first coefficient is 1, that way we can compare
	% between different autocorrelations
	% ac = ac/ac(1);
	% -> maybe actually normalise to area 1?

end
