% Class that takes a vector of features (single dimension)
% and outputs estimates of likely tempo and alignment of beats,
% under the assumption that the feature has peaks when there are
% significant events in the audio.
% A 'perfect' feature would consist of a sequence of impulses
% exactly where the beats lie, but obviously this is unattainable
% in reality

% Author: Max Fisher

classdef (Abstract) tempo_alignment_estimator < handle

properties
	% to use when writing out data files
	name;

	% processing parameters object
	params;

	% the complete matrix of features for the entire audio length.
	% Rows index feature samples over time, columns index different features
	feature_matrix;

	num_feature_samples;
	num_features;

	% output tempo and beat alignment estimates to a file with the following suffix
	data_output_suffix = '-ta-estimates.txt';

	% tempo and beat estimates estimates for each feature frame analysed
	% these are cell arrays: the n'th index of each contains a list of tuples
	% of the form
	% (tempo, tempo confidence, beat alignment, beat alignment confidence)
	% i.e. a set of estimates of possible tempos and beat alignments produced by
	% the n'th frame of the input feature. (which corresponds to a section of
	% audio that is feature_win_time seconds long, and is updated at a frequency
	% of estimate_update_rate Hz.

	% the first and third values are in samples, and must be divided by the feature
	% rate in order to get a value in seconds.
	% the second and fourth values are unitless real numbers
	tempo_alignment_estimates;
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

	function estimator = tempo_alignment_estimator(params, feature_matrix, name)
		if nargin == 0
			warning('estimator initialised with default values');
			params = {};
			feature_matrix = [];
			name = 'tempo_alignment_estimator';
		end

		estimator.params = params;
		estimator.name = name;
		estimator.feature_matrix = feature_matrix;

		estimator.num_feature_samples = size(estimator.feature_matrix, 1);
		estimator.num_features = size(estimator.feature_matrix, 2);

		estimator.tempo_alignment_estimates = ...
			cell(estimator.num_feature_frames, estimator.num_features);
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

	% gets the feature frame that ends right before the one returned by
	% this.get_feature_frame(k) begins. This is useful for doing the
	% autocorrelation, which need to slide back over previous data.
	% For early frames, there won't be enough samples to get a previous
	% nonoverlapping frame. In this case, the unavailable samples will become zeros
	function f = get_prev_nonoverlapping_frame(this, k)
		feature_win_length = this.params.feature_win_length;

		prev_frame_end_idx = this.frame_first_sample_row(k) - 1;
		prev_frame_start_idx = prev_frame_end_idx - feature_win_length;

		if prev_frame_start_idx <= 0
			% not enough previous samples; start with zeros and fill the rows that
			% can be filled
			num_zero_rows = feature_win_length - prev_frame_end_idx;
			f = zeros(feature_win_length, this.num_features);
			f((num_zero_rows + 1):end, :) = this.feature_matrix(1:prev_frame_end_idx, :);
		else
			prev_frame_rows = prev_frame_start_idx + (1:feature_win_length);
			f = this.feature_matrix(prev_frame_rows, :);
		end
	end


	% exports the computed data to files, one for each features
	function output_estimate_data(this, data_directory)
		outfile = cell(this.num_features, 1);
		for n = 1:this.num_features;
			filename = strcat(data_directory, '/', this.name, ...
				sprintf('-feature%d', n), this.data_output_suffix);
			outfile{n} = fopen(filename, 'w+');
		end

		output_format_string = strcat(...
			'frame=%d\t', ...
			'time=%.2f\t', ...
			'lag=%d\t', ...
			'bpm=%.1f\t', ...
			'tempo_conf=%.4f\t', ...
			'beat_align=%d\t', ...
			'beat_loc=%.2f\t', ...
			'alignment_conf=%.4f\n' ...
		);

		feature_sample_rate = this.params.feature_sample_rate;

		for k = 1:this.num_feature_frames
			% set time for estimate estimates from frame k
			% to be at the end of that frame
			% so first frame's estimates correspond to feature_win_time
			frame_end_time = this.params.feature_frame_end_time(k);

			for n = 1:this.num_features
				feature_n_estimates = this.tempo_alignment_estimates{k, n};

				for estimate_index = 1:size(feature_n_estimates, 1)
					% tempo and alignment are expresseed in samples,
					% convert to time by dividing by feature rate
					lag = feature_n_estimates(estimate_index, 1);
					bpm = 60*feature_sample_rate/lag;
					tempo_confidence = feature_n_estimates(estimate_index, 2);
					alignment = feature_n_estimates(estimate_index, 3);
					beat_location = frame_end_time + alignment/feature_sample_rate;
					alignment_confidence = feature_n_estimates(estimate_index, 4);

					fprintf(outfile{n}, output_format_string, ...
						k, ...
						frame_end_time, ...
						lag, ...
						bpm, ...
						tempo_confidence, ...
						alignment, ...
						beat_location, ...
						alignment_confidence ...
					);
				end
			end
		end

		for n = 1:this.num_features
			 fclose(outfile{n});
		end
	end

end % methods

methods (Abstract)
	% produce a set of tempo and beat alignment estimates from each feature, for the
	% given frame number. Also save internally in the tempo_alignment_estimates cell
	% array.
	pick_tempo_and_alignment_estimates(this, frame_number)

	% plots relevant intermediate processing data for the given sample frame
	% e.g. the graph which is peak picked to choose a tempo estimate.
	plot_frame_data(this, frame_number, feature_number)

end % methods (Abstract)

end % classdef

