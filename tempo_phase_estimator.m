% Class that takes a vector of features (single dimension)
% and outputs estimates of likely tempo and alignment of beats,
% under the assumption that the feature has peaks when there are
% significant events in the audio.
% A 'perfect' feature would consist of a sequence of impulses
% exactly where the beats lie, but obviously this is unattainable
% in reality

classdef (Abstract) tempo_phase_estimator < handle

properties (Constant)
	% minimum length of audio to use for tempo estimation at each time point (seconds)
	% actual length of feature frames will be power of 2 closest to this.
	FEATURE_WIN_TIME_APPROX = 5; %seconds

	% This can equivalently be specified using the window overlap
	%TEMPO_UPDATE_FREQUENCY_APPROX = 2; % Hz

	FEATURE_WIN_OVERLAP_PERCENT = 87.5;

	% or maybe use a window that weights recent samples more than older (by a few seconds) samples
	FEATURE_WIN_TYPE = 'rect';

	% output the tempo and beat alignment values to a text file with the following suffix
	DATA_OUTPUT_SUFFIX = '-tp-estimates.txt';

end


properties
	% name of this estimator, to use when writing out data files
	estimator_name;

	% the (single dimensional) feature vector
	feature_data;
	% effective sample rate of feature (in Hz)
	feature_sample_rate;

	% actual window length, after rounding to the nearest
	% power of 2
	feature_win_length;

	% BPM ranges to allow when detecting tempo
	min_bpm = 40;
	max_bpm = 240;

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
	tempo_phase_estimates;
end

properties (Dependent)
	% BPM ranges translated to autocorrelation time lag (seconds)
	% max lag needs to also be smaller than the frame time divided by 4, for
	% the tempo strength function to work properly
	max_lag;
	min_lag;

	% corresponding time of feature window
	feature_win_time;
	% how many samples each frame of features will advance
	feature_hop_size;

	num_feature_frames;


	% How often tempo and phase is (re)estimated on
	% a new frame of features. Given a fixed feature_sample_rate,
	% this is inversely proportional to feature_hop_size
	estimate_update_rate;

	% useful for translating between samples and time
	feature_time_axis;
end

methods
	function initialise(this, feature_data, feature_sample_rate, estimator_name)
		this.estimator_name = estimator_name;
		this.feature_data = feature_data;
		this.feature_sample_rate = feature_sample_rate;
		% make the actual window length a power of 2
		this.feature_win_length = ...
            2^round(log2(this.feature_sample_rate*this.FEATURE_WIN_TIME_APPROX));

		this.tempo_phase_estimates = cell(this.num_feature_frames, 1);
	end

	% getters for dependent properties
	function l = get.max_lag(this)
		l = 60/this.min_bpm;
	end
	function l = get.min_lag(this)
		l = 60/this.max_bpm;
	end
	function t = get.feature_win_time(this)
		t = this.feature_win_length/this.feature_sample_rate;
	end
	function h = get.feature_hop_size(this)
		h = this.feature_win_length*(1 - this.FEATURE_WIN_OVERLAP_PERCENT/100);
	end
	function f = get.estimate_update_rate(this)
		f = this.feature_sample_rate/this.feature_hop_size;
	end
	function n = get.num_feature_frames(this)
		n = ceil(length(this.feature_data)/this.feature_hop_size);
	end
	function a = get.feature_time_axis(this)
		a = (1:this.feature_win_length)/this.feature_sample_rate;
	end

	% returns a frame of feature_win_length samples
	% padded with zeros if necessary. 1-indexed.
	function feature_frame = get_feature_frame(this, k)
		feature_frame = this.feature_data((k - 1)*this.feature_hop_size + 1: ...
			min(((k - 1)*this.feature_hop_size + this.feature_win_length), end));
		if length(feature_frame) < this.feature_win_length
			feature_frame(length(feature_frame+1):this.feature_win_length) = 0;
		end
	end

	function print_properties(this)
		fprintf('Properties for %s:\n', this.estimator_name);
		fprintf('Feature window time: \t %f s \n', this.feature_win_time);
		fprintf('Feature window length: \t %d samples\n', this.feature_win_length);
		fprintf('Feature hop size: %d samples\n', this.feature_hop_size);
		fprintf('Number of feature frames: \t %d\n', this.num_feature_frames);
		fprintf('Tempo Update Frequency \t: %3.3f Hz\n', this.estimate_update_rate);
	end

	% exports the computed data to files
	function output_tempo_phase_data(this, data_directory)
		filename = strcat(data_directory, '/', this.estimator_name, ...
			this.DATA_OUTPUT_SUFFIX);
		outfile = fopen(filename, 'w+');
		
		output_format_string = strcat(...
			'frame=%d\t', ...
			'time=%f\t', ...
			'tempo=%f\t', ...
			'tempo_confidence=%f\t', ...
			'phase=%f\t', ...
			'phase_confidence=%f\n' ...
		);

		for k = 1:this.num_feature_frames
			% set time for estimate estimates from frame k
			% to be at the end of that frame
			% so first frame's estimates correspond to FEATURE_TIME
			frame_time = this.feature_win_time + (k-1)/this.estimate_update_rate;
			curr_tp_estimates = this.tempo_phase_estimates{k};
			for estimate_index = 1:length(curr_tp_estimates)
				% tempo and alignment are expresseed in samples,
				% convert to time by dividing by feature rate
				est_tempo = curr_tp_estimates(estimate_index, 1);
				est_tempo_confidence = curr_tp_estimates(estimate_index, 2);
				est_phase = curr_tp_estimates(estimate_index, 3);
				est_phase_confidence = curr_tp_estimates(estimate_index, 4);

				fprintf(outfile, output_format_string, ...
					k, ...
					frame_time, ...
					est_tempo/this.feature_sample_rate, ...
					est_tempo_confidence, ...
					est_phase/this.feature_sample_rate, ...
					est_phase_confidence ...
				);
			end
		end
	end
	
end % methods



methods (Abstract)
	% populate the tempo_phase_estimate cell array with estimates from the
	% supplied feature
	compute_tempo_phase_estimates(this)

	% plots relevant intermediate processing data for the given sample frames
	% e.g. the graph which is peak picked to choose a tempo estimate.
	plot_sample_intermediate_data(this, sample_frames)

end % methods (Abstract)

end % classdef

