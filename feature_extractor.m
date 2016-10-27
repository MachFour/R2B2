
% Class that takes audio data and produces some sort of feature which is used for
% tempo analysis by a member of the tempo_analyser class

% Handles windowing and signal processing of data, and so the relevant parameters
% are properties of this class

% Most likely all the features we will use are onset detection functions, i.e.
% they have a peak where a 'significant event' is in the music.
% Three features are planned for implementation:
% - Klapuri's multiband feature measuring change in energy across psychoacoustic frequency bands
% - Weighted phase deviation
% - The 13 lowest MFCCs excluding the zeroth/first, which measure the shape of the power spectrum on a logarithmic scale, often used in speech processing

% Abstract class because subclasses need to implement the compute feature
% function
% subclasses handle so that functions involving instances of the class
% act on references to the instance/object, not copies
classdef (Abstract) feature_extractor < handle

%
% These variables are independent of the audio length or sample rate
% For flexibility, we do not (yet) assume a fixed sample rate
% later on, we might fix it to 44100Hz (or 22050Hz, if the 2x speedup
% is necessary)
%
properties (Constant)
	% approximate time of audio to use for windowing.
	% Window size in samples will be rounded to the nearest power of 2.
	AUDIO_WIN_TIME_APPROX = 20/1000; %seconds,
	% How much each window will overlap with each other window;
	% each audio sample might be contained in more than one frame
	AUDIO_WIN_OVERLAP_PERCENT = 50;
	% (both of these values are commonly used in speech processing)

	% which type of window function to use for audio windowing
	AUDIO_WIN_TYPE = @hann

	% output the feature values to a time-value CSV file with the following suffix
	DATA_OUTPUT_SUFFIX = '-feature-data.txt';
end

% These variables are initialised when the audio and its sample rate
% are provided, using the above parameters.
% see the function initialise_audio_and_parameters
properties
	% name of this feature; to use when writing files
	feature_name;

	% the raw/uncompressed audio data
	audio_data;
	% sampled at... (in Hz)
	audio_sample_rate;

	% actual window length, after rounding to the nearest power of 2
	% at 44100Hz sample rate, this will probably be equal to 1024,
	% if AUDIO_WIN_TIME_APPROX = 20ms
	audio_win_length;

	% The matrix holding the actual feature values, after calculation. Should have num_audio_frames rows and  num_feature_channels columns
	feature_matrix;

	% actual samples of the window function, calculated when window length is known
	window_function;
end

%
% These properties are derived entirely from those above, and so cannot 
% be directly set
%
properties (Dependent)
	% actual duration of each audio window, after window length has been fixed
	audio_win_time;
	% number of samples advanced by each window
	audio_hop_size;

	num_audio_frames;

	% effective sampling rate of the feature measurement, in Hz
	feature_sample_rate;
	
	% how many feature samples are calculated
	feature_data_length;

end;

properties (Abstract)
	% dimensionality of the computed feature, i.e.
	% how many columns in the feature matrix
	num_feature_channels;
	
	
	% some features (Klapuri's) are upsampled to provide a better
	% resolution when calculating autocorrelations.
	% upsample factor is stored here
	feature_upsample_factor;
end;

methods
	% gives the instance its audio. Audio windowing could potentially
	% be separated into another layer, but some features require two
	% successive frames in order to calculate a difference value,
	% so it's not super simple. Maybe it's enough to supply a current and previous
	% frame of audio, however...?
	function initialise(this, feature_name, audio_data, audio_sample_rate)
		this.feature_name = feature_name;
		this.audio_data = audio_data;
		this.audio_sample_rate = audio_sample_rate;
		% make the actual window length a power of 2
		this.audio_win_length = ...
            2^round(log2(audio_sample_rate*this.AUDIO_WIN_TIME_APPROX));

		this.window_function = window(this.AUDIO_WIN_TYPE, this.audio_win_length);

		this.feature_matrix = zeros(this.feature_data_length, ...
			this.num_feature_channels);
	end

	% Getters for dependent properties
	function n = get.num_audio_frames(this)
		n = ceil(length(this.audio_data)/this.audio_hop_size);
	end
	function a = get.audio_win_time(this)
		a = this.audio_win_length/this.audio_sample_rate;
	end
	function r = get.feature_sample_rate(this)
		r = this.audio_sample_rate*this.feature_upsample_factor/this.audio_hop_size;
	end
	function h = get.audio_hop_size(this)
		h = this.audio_win_length*(1 - this.AUDIO_WIN_OVERLAP_PERCENT/100);
	end
	function l = get.feature_data_length(this)
		l = this.num_audio_frames*this.feature_upsample_factor;
	end

	% extracts and windows the kth audio frame from the audio data, 1-indexed
	function data = get_windowed_audio_frame(this, k)
		% if audio data has multiple colums (e.g. stereo), just return the first one
		% if we're at the end of the data, so there aren't enough samples for a full
		% frame, we just pad it out with zeros so that the returned data is always
		% audio_win_length samples in length.
		data = this.audio_data(((k-1)*this.audio_hop_size + 1): ...
			min(((k-1)*this.audio_hop_size + this.audio_win_length), end), 1);

		if length(data) < this.audio_win_length
			data(length(data)+1:this.audio_win_length) = 0;
		end

		data = data.*this.window_function;
	end

	function print_properties(this)
		fprintf('Feature extractor properties: %s\n', this.feature_name);
		fprintf('Audio sample rate: \t %3.2f kHz\n', this.audio_sample_rate/1000);
		fprintf('Audio window time: \t %3.2f ms\n', this.audio_win_time*1000);
		fprintf('Audio window length: \t %d samples\n', this.audio_win_length);
		fprintf('Audio hop size: \t %d samples\n', this.audio_hop_size);
		fprintf('Number of audio frames: \t %d\n', this.num_audio_frames);
		fprintf('Feature sample rate: \t %3.2f Hz\n', this.feature_sample_rate);
		fprintf('Dimensionality of feature: \t %d\n', this.num_feature_channels');
	end

	% export feature data to time instants to import into Sonic Visualiser
	% feature_matrix[k] represents the difference between frame k and frame k-1.
	% However the frames overlap - at what point should this measurement be
	% attributed to?
	% in terms of time, we say that this value corresponds to the end of frame k-1.
	% Then we assign feature_matrix(k) to be at (k-1)/feature_sample_rate seconds

	function output_feature_data(this, data_directory)
		outfile = cell(this.num_feature_channels, 1);
		for channel = 1:this.num_feature_channels
			filename_for_channel = strcat(data_directory, '/', ...
				this.feature_name, sprintf('-ch%d', channel), ...
				this.DATA_OUTPUT_SUFFIX);
			outfile{channel} = fopen(filename_for_channel, 'w+');
		end

		for k = 1:this.feature_data_length
			for channel = 1:this.num_feature_channels
				fprintf(outfile{channel}, '%f\t%f\n', (k-1)/this.feature_sample_rate, ...
					this.feature_matrix(k, channel));
			end
		end
	end
end

methods (Abstract)
	% returns a matrix of feature values over time, with rows indexing time/frames,
	% with period between samples corresponding to 1/feature_rate seconds

	% For features that are multi dimensional, the columns of the returned matrix
	% will index the different 'channels' of the feature. See the subclass for details.
	compute_feature(this);

	% plots some intermediate processing data as computed for each frame in sample_frames
	% e.g the power spectrum of the sample frame
	plot_sample_intermediate_data(this, sample_frames);
end

end
