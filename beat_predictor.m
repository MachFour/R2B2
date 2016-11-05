% beat_predictor.m
% Class that takes in multiple sequences of tempo and phase estimates,
% (for the same piece of music) and outputs predicted beat times
% over the course of the excerpt

classdef (Abstract) beat_predictor < handle

properties (Constant)
	% put any constant properties here
	% (for example)

	% ADAPTIVENESS_VS_CONTINUITY_RATIO = 0.5

	% how many columns to give the estimate_data cell matrix on construction
	% each column corresponds to a set of tempo phase estimates made by a different feature
	NUM_FEATURES_INITIALLY = 4;

	% output the predicted beat times to a file with this suffix
	DATA_OUTPUT_SUFFIX = '-beat-times.txt'

end % properties (Constant)

properties
	predictor_name;

	% cell array of tempo and phase estimates produced by
	% tempo_phase_estimators corresponding to several features
	% rows represent time in the excerpt (i.e. different estimate times)
	% and columns index different features (independent estimates)
	estimate_data;

	% all input tempo and phase estimates must represent
	% the same audio, and must be made the same times/rate,
	% for results to make any sense!
	first_estimate_time;
	estimate_frequency;

	% how many features' data we have stored in the cell array
	num_features;

end % Properties

properties (Dependent)
	% inverse of estimate_frequency
	time_between_estimates;

end % properties (Dependent)


methods
	function initialise(this, first_estimate_time, estimate_frequency, predictor_name)
		this.predictor_name = predictor_name;
		this.first_estimate_time = first_estimate_time;
		this.estimate_frequency = estimate_frequency;

		% initialise estimate data cell array, making room for 4 features initially
		estimate_data = cell(num_estimates, NUM_FEATURES_INITIALLY);
		num_features = 0;
	end

	% add data from a tempo_phase estimator to this predictor's set of estimates
	% make sure that the data matches the initialised first_estimate_time
	% and estimate_frequency!
	function add_feature_estimates(tempo_phase_estimates)
		num_estimate_times = size(tempo_phase_estimates, 1);
		if size(tempo_phase_estimates, 1) ~= num_estimates
			error('Wrong number of estimates for this beat_predictor!');
		else
			num_features = num_features + 1;
			estimate_data{1:, num_features} = tempo_phase_estimates;
		end
	end


	% getters for dependent variables
	function t = get.time_between_estimates(this)
		t = 1/this.estimate_frequency;
	end

end % methods

methods (Abstract)

	% for a given time in seconds, corresponding to real audio time, returns
	% the time of the next (predicted) beat
	% Note that this method is inherently non-realtime, but it serves to demonstrate
	% the functionality. Predictions for a given time should still be made only
	% using past estimate data
	predict_next_beat_time(this, audio_time)

end % methods (Abstract)

end
