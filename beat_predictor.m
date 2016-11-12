% beat_predictor.m
% Class that takes in multiple sequences of tempo and phase estimates,
% (for the same piece of music) and outputs predicted beat times
% over the course of the excerpt

% Author: Max Fisher

classdef (Abstract) beat_predictor < handle

properties (Constant)
	% put any constant properties here
	% (for example)

	% ADAPTIVENESS_VS_CONTINUITY_RATIO = 0.5

	% output the predicted beat times to a file with this suffix
	DATA_OUTPUT_SUFFIX = '-beat-times.txt'

end % properties (Constant)

properties
	predictor_name;

	% cell matrix of tempo and phase estimates produced by
	% tempo_phase_estimators corresponding to several features.

	% Tempo and phase estimates should be in seconds!!!
	% Otherwise this class has to know about the feature sampling
	% rate, which it doesn't need to

	% Rows represent time in the excerpt (i.e. estimates that are in the
	% same row - same feature or not - should correspond to the same
	% time in audio terms. Columns index different features
	estimate_data;

	% list of winning estimate tuples, one for each row in the
	% estimate data matrix. Note that the tuples don't have
	% to be actual estimates, they can be averages (for example)
	winning_estimates;

	% all input tempo and phase estimates must represent
	% the same audio, and must be made the same times/rate,
	% for results to make any sense!
	first_estimate_time;
	estimate_frequency;

	% array of output beat times, computed after the winning tempo and
	% phase is known for each row of the estimate data
	beat_times;

end % Properties

properties (Dependent)
	% inverse of estimate_frequency
	time_between_estimates;

	% columns in cell matrix
	num_features;

end % properties (Dependent)


methods
	% Give this class the estimate data (which should not use samples but time for tempo and phase
	% estimates
	% first_estimate_time is the audio time to which the first row of estimate groups corresponds to
	% (i.e. the end of the first feature frame)
	% estimate_frequency is 1 divided by the time between the ends of successive feature frames
	function initialise(this, estimate_data, first_estimate_time, estimate_frequency, predictor_name)
		this.predictor_name = predictor_name;
		this.first_estimate_time = first_estimate_time;
		this.estimate_frequency = estimate_frequency;

		% initialise estimate data cell array, making room for 4 features initially
		this.estimate_data = estimate_data;

		% sequence of winning tempo/phase estimate tuples
		this.winning_estimates = zeros(size(estimate_data, 1), 4);
	end

	% return the audio time (in seconds), corresponding to when the
	% kth row of estimates in the estimate_data matrix (or the kth
	% winning estimate) were made
	function t = get_estimate_time(this, k)
		% this is equal to feature_win_time of a tempo_phase_estimator
		estimate_t0 = this.get_first_estimate_time;
		% this is equal to 1/estimate_update_rate of a tempo_phase_estimator
		estimate_dt = this.get_time_between_estimates;
		t = estimate_t0 + (k-1)*estimate_dt;
	end

	% getters for dependent variables
	function t = get.time_between_estimates(this)
		t = 1/this.estimate_frequency;
	end

	function n = get.num_features(this)
		n = size(this.estimate_data, 2);
	end

	% output the beat times for this predictor, based on the winning sequence
	% of tempo and phase estimates
	% Note that this method is inherently non-realtime, but it serves to demonstrate
	% the functionality. Predictions for a given time should still be made only
	% using past estimate data
	function compute_beat_times(this)
		% assume 2 beats predicted per row of estimate data,
		% just to be able to preallocate something
		this.beat_times = zeros(2*size(this.winning_estimates, 1));
		beat_index = 1;

		for k = 1:size(this.estimate_data, 1)
			% output beat predictions using the kth winning tempo/phase
			% estimate, up to the estimate time of the next frame
			estimate_time = this.get_estimate_time(k);
			% the estimates for the final frame may extend slightly over the
			% end of the audio
			next_estimate_time = this.get_estimate_time(k+1);

			% these should be in seconds. Phase should be negative
			kth_tempo = this.winning_estimates(k, 1);
			kth_phase = this.winning_estimates(k, 3);

			latest_observed_beat_time = estimate_time + kth_phase;

			predicted_beat_time = latest_observed_beat_time + kth_tempo;
			% do we need to allow for latency?
			overlap_allowed = 0.01;
			while predicted_beat_time <= next_estimate_time + overlap_allowed
				this.beat_times(beat_index) = predicted_beat_time;
				beat_index = beat_index + 1;
				% try to predict another beat time
				predicted_beat_time = predicted_beat_time + kth_tempo;
			end
		end
		% trim beat times array to remove zeros?
	end

	function output_beat_times(this, data_directory)
		filename = strcat(data_directory, '/', this.predictor_name, ...
			this.DATA_OUTPUT_SUFFIX);
		outfile = fopen(filename, 'w+');

		for beat_index = 1:length(this.beat_times)
			beat_time = this.beat_times(beat_index);
			if beat_index > 1 && beat_time ~= 0
				fprintf(outfile, '%.3f\n', this.beat_times(beat_index));
			end
		end
	end

end % methods

methods (Abstract)
	% basically populate the winning_estimates array
	compute_winning_tp_estimates(this)
end % methods (Abstract)

end
