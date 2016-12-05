% beat_predictor.m
% Class that takes in multiple sequences of tempo and phase estimates,
% (for the same piece of music) and outputs predicted beat times
% over the course of the excerpt

% Author: Max Fisher

classdef (Abstract) beat_predictor < handle

properties (Constant)
	% output the predicted beat times to a file with this suffix
	DATA_OUTPUT_SUFFIX = '-beat-times.txt'

end % properties (Constant)

properties

end % Properties

methods
	function initialise(this, params, predictor_name, num_features)
	end


end % methods
end
