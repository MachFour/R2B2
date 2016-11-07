% odf_klapuri
% takes in the entire sample of audio, performs windowing
% differentials of the lowest Mel-Frequency Cepstral Coefficients
% over time

% Author: Max Fisher

classdef odf_mfcc < feature_extractor

properties
    % low pass filter for interpolating the log magnitude spectrum
	lp_cutoff_freq = 10; %Hz
	lp_order = 6;

	num_mel_filters = 36;

	% parameter for log-compression
	mu = 100;

	% parameter for weighted sum of differenced feature with non-differenced
	lambda = 1;

	% number of MFCCs used
	num_feature_channels = 6;

    % following Klapuri, we double the feature rate by interpolating
	% the values calculated for each frame by 2
	feature_upsample_factor = 2;

	% stores previously computed values, needed to take a difference
	frame_mfccs;

	% Mel filterbank: used for MFCC's and also for Klapuri's features
	% returns a num_mel_filters*AUDIO_WIN_LENGTH/2 size matrix,
	% perfect for left multiplying the DFT magnitude spectrum by
	mel_filterbank;

	% filter for doing interpolation of feature data
	LP_num;
	LP_denom;
	% when saving the final feature value, we multiply by this ratio
	% to compensate for the small gain of the low pass filter
	LP_gain;

end % properties

properties (Dependent)
	% just this.num_feature_channels
	num_feature_mfccs;
end % properties (Dependent)

methods
	% print feature specific properties
	function print_properties(this)
		print_properties@feature_extractor(this);
		fprintf('Upsampling factor: \t %d\n', this.feature_upsample_factor);
		fprintf('Number of filters used in Mel filterbank: \t %d\n', this.num_mel_filters);
	end

	% set up feature specific parameters
	function initialise(this, audio_data, audio_sample_rate, name)
		initialise@feature_extractor(this, audio_data, audio_sample_rate, name);

		% in units of pi/2 radians per sample
		LP_cutoff_freq_normalised = this.lp_cutoff_freq/(this.feature_sample_rate/2);

		[this.LP_num, this.LP_denom] = butter(this.lp_order, LP_cutoff_freq_normalised);
		this.LP_gain = this.feature_sample_rate/this.lp_cutoff_freq;

		% stored the (compressed, upsampled and filtered)
		% mel band energy in successive frames
		this.frame_mfccs = ...
			zeros(this.feature_upsample_factor*this.num_audio_frames, this.num_mel_filters);

		this.mel_filterbank = mel_filter(this.audio_sample_rate, ...
			this.audio_win_length, this.num_mel_filters);
	end


	function compute_feature(this)
		% store filter data in between loop iterations
		filter_delays = zeros(this.lp_order, this.num_mel_filters);

		for n = 1:this.num_audio_frames
			curr_frame = this.get_windowed_audio_frame(n);

			%avg = sum(curr_frame)/length(curr_frame);
			%variance = sum((curr_frame - avg).^2)/(length(curr_frame) - 1);
			% normalise to mean 0 variance 1
			% does doing this introduce noise?
            % since successive frames would have different means and variances...

			%curr_frame = (curr_frame - avg)/sqrt(variance);

			power_spectrum = abs(fft(curr_frame, this.audio_win_length));
			% throw away (symmetrical) right half of spectrum
			power_spectrum = power_spectrum(1:this.audio_win_length/2);

			mel_band_energy = this.mel_filterbank * power_spectrum;

			% mel band energy is a vector of length this.num_mel_filters

			% use log-compression as Klapuri does, mu = 100
			compressed_band_energy = log(1 + this.mu*mel_band_energy)/log(1 + this.mu);

			% here Klapuri performs a 2x interpolation of the compressed
			% mel band energy over time before calculating the difference between frames

			% put feature data as a row and feature_upsample_factor - 1 rows of zeros into
			% a matrix to be low pass filtered

			if this.feature_upsample_factor > 1
				upsampled_band_energy = [
					compressed_band_energy'
					zeros(this.feature_upsample_factor - 1, this.num_mel_filters)
				];

				% filter along first dimension (frames), using the low pass filter, and store the
				% filter delays for the next loop iteration
				[interpolated_band_energy, filter_delays] = ...
					filter(this.LP_num, this.LP_denom, upsampled_band_energy, filter_delays, 1);
			else
				interpolated_band_energy = compressed_band_energy';
			end

			% correct slice of length this.feature_upsample_factor
			% e.g. for feature_upsample_factor = 2, this will be
			% 2*n-1:2*n
			save_indices = (this.feature_upsample_factor*(n-1)+1):this.feature_upsample_factor*n;

			% dct operates column wise
			mfccs = (dct(interpolated_band_energy'))';
			this.frame_mfccs(save_indices, :) = mfccs;

			% now go through the interpolated energy vector and take 1-sample
			% differences, adding positive differences from adjacent bands
			% together to create feature channels, where each one peaks at
			% energy rises occuring somewhere in their set of mel frequencies

			% for each audio frame we have this.num_feature_channels
			% rows of frame_mfccs to process. Do this matrix-wise

			difference_indices = save_indices - 1;
			if n == 1
				 % first frame will just be subtracted from itself
				difference_indices(1) = 1;
			end

			unrectified_difference = mfccs - this.frame_mfccs(difference_indices, :);

			% take absolute value of difference
			abs_diff_mfccs = abs(unrectified_difference);
			%abs_diff_mfccs(abs_diff_mfccs < 0) = 0;

			% weighted sum the result with the undifferentiated log_power_spectrum
			% makes analysis better according to [4]
			% also multiply differenced filter output by ratio of
			% sampling frequencies of audio and low pass filter, (again following [4])
			this.feature_matrix(save_indices, :) = ...
				(1 - this.lambda)*mfccs(:, 1:this.num_feature_mfccs) + ...
				this.lambda*this.LP_gain*abs_diff_mfccs(:, 1:this.num_feature_mfccs);
		end
	end
	
	function n = get.num_feature_mfccs(this)
		n = this.num_feature_channels;
	end

	function plot_sample_intermediate_data(this, sample_frames)
		for frame = sample_frames
			% recalculate feature differential; but if frame is 1 there's
			% nothing to take a differential with
			if frame == 1
				warning('Cannot calculate feature differential for frame 1');
				unrectified_difference = this.frame_mfccs(frame, :);
			else
				unrectified_difference = this.frame_mfccs(frame, :) ...
				- this.frame_mfccs(frame-1, :);
			end

			HWR_difference = this.LP_gain*unrectified_difference;
			HWR_difference(HWR_difference < 0) = 0;

			figure;
			subplot(2, 1, 1);
			plot(this.frame_mfccs(frame, :));
			title(sprintf('mel band energies, frame %d', frame));
			subplot(2, 1, 2);
			plot(HWR_difference);
			title(sprintf('Rectified difference, frame %d->%d', frame-1, frame));
		end
	end

end % methods

end % classdef
