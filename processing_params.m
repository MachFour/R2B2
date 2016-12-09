
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
	max_tempo_peaks;
	max_phase_peaks;

	
	
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

	% How often tempo and phase is (re)estimated on
	% a new frame of features. Given a fixed feature_sample_rate,
	% this is inversely proportional to feature_hop_size
	prediction_rate;

	% inverse of prediction rate
	time_between_predictions;

	% how long until the first prediction is output by the algorithm, in audio terms
	first_prediction_time;


end % properties


methods

	% PARAMETERS:
	% audio_sample_rate
	% 	is the sampling rate of the input audio signal, in Hz.
	% audio_win_time
	% 	is the approximate time of audio to use for windowing.
	% 	 Window size in samples will be rounded to the nearest power of 2.
	% audio_win_overlap_percent
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
	% feature_win_overlap_percent
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
			audio_win_time, audio_win_overlap_percent, audio_win_type, ...
			feature_win_time, feature_win_overlap_percent, feature_win_type, ...
			feature_upsample_factor, ...
			min_bpm, max_bpm, max_tempo_peaks ...
		)

		% default values
		if nargin == 0
			warning('Setting default sample rate of 44.1 kHz');
			audio_sample_rate = 44100;
			% (the following three values are commonly used in speech processing)
			audio_win_time = 20/1000; %seconds,
			audio_win_overlap_percent = 50;
			audio_win_type = @hann;

			min_bpm = 35;
			max_bpm = 320;

			max_tempo_peaks = 8;
			max_phase_peaks = 4;

			feature_win_time = 5; %seconds
			feature_win_overlap_percent = 87.5;
			% or maybe use a window that weights recent samples more than older (by a
			% few seconds) samples
			feature_win_type = 'rect';

			feature_upsample_factor = 2;
		end

		% stuff needed to calculate other properties
		audio_win_length = 2^round(log2(audio_sample_rate*audio_win_time));
		audio_hop_size = round(audio_win_length*(1 - audio_win_overlap_percent/100));

		feature_sample_rate = audio_sample_rate*feature_upsample_factor/audio_hop_size;
		feature_win_length = 2^round(log2(feature_sample_rate*feature_win_time));
		feature_hop_size = round( ...
			feature_win_length*(1 - feature_win_overlap_percent/100));

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
		params.max_phase_peaks = max_phase_peaks;
		
		
		% precompute dependent features
		
		params.set_dependent_features;
		
	end

	% return the audio time (in seconds) that the predictions for feature frame k
	% will be made. Equivalently, it is the end time of the kth feature frame.
	function t = prediction_time(this, k)
		estimate_t0 = this.first_prediction_time;
		estimate_dt = this.time_between_predictions;
		t = estimate_t0 + (k-1)*estimate_dt;
	end

	function set_dependent_features(this)
	
		this.audio_win_time = this.audio_win_length/this.audio_sample_rate;

		this.max_lag = this.max_lag_samples/this.feature_sample_rate;
	
		this.min_bpm = 60/this.max_lag;

		this.min_lag = this.min_lag_samples/this.feature_sample_rate;
		
		this.max_bpm = 60/this.min_lag;

		this.tempo_lag_range = (this.min_lag_samples:this.max_lag_samples)';

		this.feature_win_time = this.feature_win_length/this.feature_sample_rate;

		this.prediction_rate = this.feature_sample_rate/this.feature_hop_size;

		this.first_prediction_time = this.feature_win_time;
	
		this.time_between_predictions = this.feature_hop_size/this.feature_sample_rate;
	end
	

end % methods

end % classdef


