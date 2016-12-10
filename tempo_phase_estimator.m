% Class that takes a vector of features (single dimension)
% and outputs estimates of likely tempo and alignment of beats,
% under the assumption that the feature has peaks when there are
% significant events in the audio.
% A 'perfect' feature would consist of a sequence of impulses
% exactly where the beats lie, but obviously this is unattainable
% in reality

% Author: Max Fisher

classdef (Abstract) tempo_phase_estimator < handle

properties
	% name of this estimator, to use when writing out data files
	estimator_name;

	% processing parameters object
	params;

	% the complete matrix of features for the entire audio length.
	% Rows index feature samples over time, columns index different features
	feature_matrix;

	num_feature_samples;
	num_features;

	% output tempo and beat alignment estimates to a file with the following suffix
	data_output_suffix = '-tp-estimates.txt';

	% tempo and beat phase estimates for each feature frame analysed
	% these are cell arrays: the n'th index of each contains a list of tuples
	% of the form
	% (tempo, tempo confidence, beat location, beat location confidence)
	% i.e. a set of estimates of possible tempo and beat phases produced by
	% the n'th frame of the input feature. (which corresponds to a section of
	% audio that is feature_win_time seconds long, and is updated at a frequency
	% of estimate_update_rate Hz.

	% the first and third values are in samples, and must be multiplied by the feature
	% rate in order to get a value in seconds.
	% the second and fourth values are unitless real numbers
	tp_estimates;
end

properties (Dependent)
	num_feature_frames;
end


methods
	function n = get.num_feature_frames(this)
		feature_win_length = this.params.feature_win_length;
		feature_hop_size = this.params.feature_hop_size;
		% number of feature frames k, is the smallest number such that
		%feature_win_length + k*feature_hop_size > num_feature_samples;

		n = ceil((this.num_feature_samples - feature_win_length)/feature_hop_size);

	end

	function initialise(this, params, feature_matrix, estimator_name)
		this.params = params;
		this.estimator_name = estimator_name;
		this.feature_matrix = feature_matrix;

		this.num_feature_samples = size(this.feature_matrix, 1);
		this.num_features = size(this.feature_matrix, 2);

		this.tp_estimates = cell(this.num_feature_frames, this.num_features);
	end

	% returns the index of the first sample in the kth feature frame, to index into
	% rows of the featurematrix
	function i = frame_first_sample_row(this, k)
		i = (k - 1)*this.params.feature_hop_size + 1;
	end

	function i = frame_last_sample_row(this, k)
		i = (k - 1)*this.params.feature_hop_size + this.params.feature_win_length;
	end

	% returns a frame of feature_win_length samples
	% padded with zeros if necessary. 1-indexed.
	function f = get_feature_frame(this, k)
		f = zeros(this.params.feature_win_length, this.num_features);

		frame_data = this.feature_matrix(this.frame_first_sample_row(k): ...
			min(this.frame_last_sample_row(k), end), :);

		% if there is not enough feature data left to finish the frame, the rest is
		% just filled with zeros.
		f(1:size(frame_data, 1), :) = frame_data;
	end

	% exports the computed data to files, one for each features
	function output_tempo_phase_data(this, data_directory)
		outfile = cell(this.num_features, 1);
		for n = 1:this.num_features;
			filename = strcat(data_directory, '/', this.estimator_name, ...
				sprintf('-feature%d', n), this.data_output_suffix);
			outfile{n} = fopen(filename, 'w+');
		end

		output_format_string = strcat(...
			'frame=%d\t', ...
			'time=%f\t', ...
			'tempo=%f\t', ...
			'tempo_confidence=%f\t', ...
			'phase=%f\t', ...
			'phase_confidence=%f\n' ...
		);

		feature_sample_rate = this.params.feature_sample_rate;

		for k = 1:this.num_feature_frames
			% set time for estimate estimates from frame k
			% to be at the end of that frame
			% so first frame's estimates correspond to feature_win_time
			frame_end_time = this.params.feature_frame_end_time(k);

			for n = 1:this.num_features
				curr_tp_estimates = this.tp_estimates{k, n};

				for estimate_index = 1:size(curr_tp_estimates, 1)
					% tempo and alignment are expresseed in samples,
					% convert to time by dividing by feature rate
					est_tempo = curr_tp_estimates(estimate_index, 1)/feature_sample_rate;
					est_tempo_confidence = curr_tp_estimates(estimate_index, 2);
					est_phase = curr_tp_estimates(estimate_index, 3)/feature_sample_rate;
					est_phase_confidence = curr_tp_estimates(estimate_index, 4);

					fprintf(outfile{n}, output_format_string, ...
						k, ...
						frame_end_time, ...
						est_tempo, ...
						est_tempo_confidence, ...
						est_phase, ...
						est_phase_confidence ...
					);
				end
			end
		end
	end

end % methods

methods (Abstract)
	% populate the tempo_phase_estimate cell array with estimates from the
	% supplied feature
	compute_tempo_phase_estimates(this, frame_idx)

	% plots relevant intermediate processing data for the given sample frames
	% e.g. the graph which is peak picked to choose a tempo estimate.
	plot_sample_intermediate_data(this, sample_frames, feature_idx)

end % methods (Abstract)

end % classdef

