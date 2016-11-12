% autocorrelation.m
% calculates the autocorrelation of a feature frame
% by sliding it back across itself and a previous frame
% this way problems of bias are avoided

% frames should be column vectors and have the same length
% make sure the frames have the same length

% Author: Max Fisher

function ac = autocorrelation(feature_frame, prev_frame)
	if size(feature_frame) ~= size(feature_frame)
		error('Frames are different sizes');
	end
	l = length(feature_frame);
	combined_frame = [prev_frame; feature_frame];

	% normalise to power 1 (not mean zero)
	power = sum(combined_frame.^2)/length(combined_frame);
	combined_frame = combined_frame/sqrt(power);
	feature_frame = feature_frame/sqrt(power);

	ac = zeros(size(feature_frame));

	for shift = 0:l-1
		shifted_frame = combined_frame(l-shift+1:end - shift);
		ac(shift+1) = sum(shifted_frame.*feature_frame);
	end

	% normalise so that first coefficient is 1, that way we can compare
	% between different autocorrelations
	ac = ac/ac(1);
end
