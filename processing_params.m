
% holds all the important global parameters for the algorithm

classdef processing_params < handle

properties
	audio_sample_rate;

	% which type of window function to use for audio windowing
	audio_win_type;

	% window type used for windowing feature data into feature frames
	% this variable is currently unused
	feature_win_type;

	% some features (Klapuri's) are upsampled to provide a better
	% resolution when calculating autocorrelations.
	% upsample factor is stored here
	feature_upsample_factor;


	% actual window lengths, after rounding to the nearest power of 2
	audio_win_length;
	feature_win_length;

	% number of samples advanced by each audio window
	audio_hop_size;

	% how many samples each frame of features will advance
	feature_hop_size;

	% effective sampling rate of the feature measurement, in Hz
	feature_sample_rate;

	% BPM ranges translated to autocorrelation time lag (samples)
	% derived from {min,max}_bpm_approx respectively, and the feature
	% sampling rate. Rounded to the nearest integer, of course.
	max_lag_samples;
	min_lag_samples;

	% how many peaks to pick in tpe_autocorrelation and friends
	% the maximum total number of tempo/beat alignment estimates for each feature
	% frame is given by MAX_TEMPO_PEAKS*MAX_PHASE_PEAKS
	max_tempo_peaks;
	max_alignment_peaks;



	%
	% The following properties depend entirely on the previous ones, but since
	% the parameters do not change after initalisation, they are precomputed and
	% stored as normal properties (not dependent ones)
	%

	% actual duration of each audio window, after window length has been fixed
	audio_win_time;

	% BPM ranges translated to autocorrelation time lag (seconds)
	% these are equivalent to {max,min}_lag*feature_sample_rate,
	max_lag;
	min_lag;

	% recomputed min and max bpm after rounding to nearest sample
	min_bpm;
	max_bpm;

	% range of allowable tempos (in samples)
	% this is simply min_lag_samples:max_lag_samples
	% used to e.g. limit range of peak picking in autocorrelation
	tempo_lag_range;

	% corresponding time of feature window
	feature_win_time;

	% How often tempo and beat alignment is (re)estimated on a new frame of features.
	% Given a fixed feature_sample_rate, this is proportional to feature_hop_size
	time_between_estimates;

	% how long until the first estimate is output by the algorithm, in audio terms
	first_estimate_time;


end % properties


methods

	% PARAMETERS:
	% audio_sample_rate
	% 	is the sampling rate of the input audio signal, in Hz.
	% audio_win_time
	% 	is the approximate time of audio to use for windowing.
	% 	 Window size in samples will be rounded to the nearest power of 2.
	% audio_win_overlap_proportion
	% 	is how much each window will overlap with each other window;
	% 	each audio sample might be contained in more than one frame
	% audio_win_type
	% 	is which type of window function to use for audio windowing
	% min_bpm, max_bpm
	% 	are self explanatory. Actual values are rounded to nearest sample
	% 	equivalents.
	% feature_win_time
	% 	is the approximate length of audio to use for tempo estimation at each time
	% 	point (in seconds). The actual length of feature frames will be the closest
	% 	power of 2 to this number.
	% feature_win_overlap_proportion
	% 	defines how much feature frames overlap
	% feature_win_type
	% 	is the windowing method used to create feature frames from the continuous
	% 	onset detection function samples
	% feature_upsample_factor
	% 	is the factor by which features are upsampled before further processing.
	% max_tempo_peaks
	% 	is how many tempo peaks to pick, per feature, which computing tempo estimates
	function params = processing_params( ...
			audio_sample_rate, ...
			audio_win_time, audio_win_overlap_proportion, audio_win_type, ...
			feature_win_time, feature_win_overlap_proportion, feature_win_type, ...
			feature_upsample_factor, ...
			min_bpm, max_bpm, max_tempo_peaks, max_alignment_peaks ...
		)

		% default values
		if nargin == 0
			warning('Setting default sample rate of 44.1 kHz');
			audio_sample_rate = 44100;
			% (the following three values are commonly used in speech processing)
			audio_win_time = 20/1000; %seconds,
			audio_win_overlap_proportion = 0.5;
			audio_win_type = @hann;

			min_bpm = 35;
			max_bpm = 320;

			max_tempo_peaks = 8;
			max_alignment_peaks = 4;

			feature_win_time = 3; %seconds
			feature_win_overlap_proportion = 0.75;
			% or maybe use a window that weights recent samples more than older (by a
			% few seconds) samples
			feature_win_type = 'rect';

			feature_upsample_factor = 2;
		end

		% stuff needed to calculate other properties
		audio_win_length = 2^round(log2(audio_sample_rate*audio_win_time));
		audio_hop_size = round(audio_win_length*(1 - audio_win_overlap_proportion));

		feature_sample_rate = audio_sample_rate*feature_upsample_factor/audio_hop_size;
		feature_win_length = 2^round(log2(feature_sample_rate*feature_win_time));
		feature_hop_size = round(feature_win_length*(1-feature_win_overlap_proportion));

		% rounding of supplied approximate values
		params.min_lag_samples = round(feature_sample_rate*60/max_bpm);
		params.max_lag_samples = round(feature_sample_rate*60/min_bpm);

		% fill in other values
		params.audio_sample_rate = audio_sample_rate;
		params.audio_win_length = audio_win_length;
		params.audio_hop_size = audio_hop_size;

		params.feature_sample_rate = feature_sample_rate;
		params.feature_win_length = feature_win_length;
		params.feature_hop_size = feature_hop_size;
		params.feature_upsample_factor = feature_upsample_factor;

		params.audio_win_type = audio_win_type;
		params.feature_win_type = feature_win_type;

		params.max_tempo_peaks = max_tempo_peaks;
		params.max_alignment_peaks = max_alignment_peaks;


		% precompute dependent features

		params.set_dependent_features;

	end

	% return the audio time (in seconds) that the estimates for feature frame k
	% will be made. Equivalently, it is the end time of the kth feature frame.
	function t = estimate_time(this, k)
		t = this.feature_frame_end_time(k);
	end

	% returns the time when frame k ends, in seconds
	% first frame ends at this.feature_win_time
	function t = feature_frame_end_time(this, k)
		t = this.feature_win_time + (k-1)*this.time_between_estimates;
	end

	function t = feature_frame_start_time(this, k)
		t = (k-1)*this.time_between_estimates;
	end

	function a = feature_time_axis(this)
		a = (1:this.feature_win_length)/this.feature_sample_rate;
	end

	function set_dependent_features(this)

		this.audio_win_time = this.audio_win_length/this.audio_sample_rate;

		this.max_lag = this.max_lag_samples/this.feature_sample_rate;

		this.min_bpm = 60/this.max_lag;

		this.min_lag = this.min_lag_samples/this.feature_sample_rate;

		this.max_bpm = 60/this.min_lag;

		this.tempo_lag_range = (this.min_lag_samples:this.max_lag_samples);

		this.feature_win_time = this.feature_win_length/this.feature_sample_rate;

		this.first_estimate_time = this.feature_win_time;

		this.time_between_estimates = this.feature_hop_size/this.feature_sample_rate;
	end


end % methods

end % classdef


