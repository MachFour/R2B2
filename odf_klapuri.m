% odf_klapuri
% takes in the entire sample of audio, performs windowing
% and calculates the energy in critical bands over time
% As done in Klapuri et al. (2006), "Analysis of the meter of acoustic musical signals"

classdef odf_klapuri < feature_extractor

properties (Constant)
    % low pass filter for interpolating the log magnitude spectrum
	LP_CUTOFF_FREQ = 15; %Hz
	LP_ORDER = 6;

	NUM_MEL_FILTERS = 36;

	% parameter for log-compression
	MU = 100;

	% parameter for weighted sum of differenced feature with non-differenced
	LAMBDA = 0.8;

end

properties
	num_feature_channels = 4;

    % following Klapuri, we double the feature rate by interpolating
	% the values calculated for each frame by 2
	feature_upsample_factor = 2;

	% stores previously computed values, needed to take a difference
	compressed_mel_band_energy;
	interpolated_mel_band_energy;

	% Mel filterbank: used for MFCC's and also for Klapuri's features
	% returns a NUM_MEL_FILTERS*AUDIO_WIN_LENGTH/2 size matrix,
	% perfect for left multiplying the DFT magnitude spectrum by
	mel_filterbank;
	% handle to get vector of which filters correspond to which channel
	mel_filters_in_channel;
	
	% filter for doing interpolation of feature data
	LP_num;
	LP_denom;
	% when saving the final feature value, we multiply by this ratio
	% to compensate for the small gain of the low pass filter
	LP_compensation_factor;

    % divide up filters evenly into NUM_FEATURE_CHANNELS channels
	filter_divisions;
end


methods
	% print feature specific properties
	function print_properties(this)
		print_properties@feature_extractor(this);
		fprintf('Upsampling factor: \t %d\n', this.feature_upsample_factor);
		fprintf('Number of filters used in Mel filterbank: \t %d\n', this.NUM_MEL_FILTERS);
	end

	% set up feature specific parameters
	function initialise(this, name, audio_data, audio_sample_rate)
		initialise@feature_extractor(this, name, audio_data, audio_sample_rate);

		% in units of pi/2 radians per sample
		LP_cutoff_freq_normalised = this.LP_CUTOFF_FREQ/(this.feature_sample_rate/2);

		[this.LP_num, this.LP_denom] = butter(this.LP_ORDER, LP_cutoff_freq_normalised);
		this.LP_compensation_factor = this.feature_sample_rate/this.LP_CUTOFF_FREQ;

		this.compressed_mel_band_energy = ...
			zeros(this.num_audio_frames, this.NUM_MEL_FILTERS);
		this.interpolated_mel_band_energy = ...
			zeros(this.feature_upsample_factor*this.num_audio_frames, this.NUM_MEL_FILTERS);

		this.mel_filterbank = mel_filter(this.audio_sample_rate, ...
			this.audio_win_length, this.NUM_MEL_FILTERS);

		this.filter_divisions = round(linspace(1,this.NUM_MEL_FILTERS, ...
			this.num_feature_channels+1));
		
		this.mel_filters_in_channel = @(c) this.filter_divisions(c):this.filter_divisions(c+1);

	end;


	function compute_feature(this)
		for k = 1:this.num_audio_frames
			curr_frame = this.get_windowed_audio_frame(k);

			mean = sum(curr_frame)/length(curr_frame);
			variance = sum(curr_frame.^2)/(length(curr_frame) - 1);
			% normalise to mean 0 variance 1
			% does doing this introduce noise?
            % since successive frames would have different means and variances...

			%curr_frame = (curr_frame - mean)/sqrt(variance);

			power_spectrum = abs(fft(curr_frame, this.audio_win_length));
			% throw away (symmetrical) right half of spectrum
			power_spectrum = power_spectrum(1:this.audio_win_length/2);

			mel_band_energy = this.mel_filterbank * power_spectrum;

			% use log-compression as Klapuri does, mu = 100
			this.compressed_mel_band_energy(k, :) = log(1 + this.MU*mel_band_energy)/log(1 + this.MU);
		end

		% here Klapuri performs a 2x interpolation of the compressed
		% mel band energy over time before calculating the difference between frames
		% have to find a way of doing this inside the loop
		this.interpolated_mel_band_energy = upsample(this.compressed_mel_band_energy, this.feature_upsample_factor);
		% do this with [Y, Zf] = filter(B, A, X, Zi, 1);
		this.interpolated_mel_band_energy = filter(this.LP_num, this.LP_denom, this.interpolated_mel_band_energy, [], 1);

		for k = 2:this.feature_data_length
			% take the half-wave rectified difference between this and prev frame's
			% log power spectrum, to create a feature that peaks at the onset time

			% Bock: second difference so that we take the difference with a 50%
			% overlapping window, but have double the feature rate

			unrectified_difference = this.interpolated_mel_band_energy(k, :) ...
				- this.interpolated_mel_band_energy(k - 1, :);

			% half-wave rectified difference: set all negative entries to 0
			HWR_difference = unrectified_difference;
			HWR_difference(HWR_difference < 0) = 0;

			for c = 1:this.num_feature_channels
				mel_filters_to_sum = this.mel_filters_in_channel(c);
				undifferenced_band_energy = ...
					sum(this.interpolated_mel_band_energy(k, mel_filters_to_sum), 2);
				differenced_band_energy = ...
					sum(HWR_difference(mel_filters_to_sum));
				
				% weighted sum the result with the undifferentiated log_power_spectrum
				% makes analysis better according to [4]
				% also multiply differenced filter output by ratio of
				% sampling frequencies of audio and low pass filter, (again following [4])
				this.feature_matrix(k, c) = (1 - this.LAMBDA)*undifferenced_band_energy + ...
					this.LAMBDA*this.LP_compensation_factor*differenced_band_energy;
			end
		end


	end


	function plot_sample_intermediate_data(this, sample_frames)
		for frame = sample_frames
			% recalculate feature differential; but if frame is 1 there's
			% nothing to take a differential with
			if frame == 1
				warning('Cannot calculate feature differential for frame 1');
				unrectified_difference = this.interpolated_mel_band_energy(frame, :);
			else
				unrectified_difference = this.interpolated_mel_band_energy(frame, :) ...
				- this.interpolated_mel_band_energy(frame-1, :);
			end

			HWR_difference = this.LP_compensation_factor*unrectified_difference;
			HWR_difference(HWR_difference < 0) = 0;

			figure;
			subplot(2, 1, 1);
			plot(this.interpolated_mel_band_energy(frame, :));
			title(sprintf('mel band energies, frame %d', frame));
			subplot(2, 1, 2);
			plot(HWR_difference);
			title(sprintf('LP compensated HWR difference, frame %d', frame));
		end

	end
	end

end
