classdef odf_wpd < feature_extractor
    
properties
    num_feature_channels = 1;
    
    num_feature_dimensions = 1;
    
    feature_upsample_factor = 2;
    
    frame_wpd;
    
    lp_cutoff_freq = 10; %Hz
	lp_order = 6;
    
    % filter for doing interpolation of feature data
	LP_num;
	LP_denom;
	% when saving the final feature value, we multiply by this ratio
	% to compensate for the small gain of the low pass filter
	LP_gain;
    
end

methods 
    
    function initialise(this, audio_data, audio_sample_rate, name)
        initialise@feature_extractor(this, audio_data, audio_sample_rate, name);
            
        this.frame_wpd = ...
            zeros(this.feature_upsample_factor*this.num_audio_frames,this.num_feature_dimensions);
        
        LP_cutoff_freq_normalised = this.lp_cutoff_freq/(this.feature_sample_rate/2);

		[this.LP_num, this.LP_denom] = butter(this.lp_order, LP_cutoff_freq_normalised);
		this.LP_gain = this.feature_sample_rate/this.lp_cutoff_freq;
    end
    
    function compute_feature(this)
        % Store delays between each iteration
        % Used for LP filter when interpolating
        filter_delays = zeros(this.lp_order, this.num_feature_dimensions);
        % x_ph_prev stores the phase info of the once time delayed audio
        % frame
        x_ph_prev = zeros(this.audio_win_length/2, this.num_feature_dimensions);
        % x_ph_first_diff_prev contains the first order difference between
        % the phase of the once time delayed and twice time delayed audio
        % frame (w.r.t current frame)
        x_ph_first_diff_prev = zeros(this.audio_win_length/2, this.num_feature_dimensions);
        % n frames, k frequency bins
        for n = 1:this.num_audio_frames
            
           % Fetch audio frame
           curr_frame = this.get_windowed_audio_frame(n);
           % Initialise wpd
           wpd = 0;
           % Spectrum of current frame
           x_fft = fft(curr_frame, this.audio_win_length);
           x_fft = x_fft(1:this.audio_win_length/2);
           % Phase of current frame
           x_ph_curr = angle(x_fft);
           % First order difference between current phase and once time
           % delayed phase
           x_ph_first_diff_cur = x_ph_curr - x_ph_prev;
           % Second order difference between current phase and once time
           % delayed phase
           x_ph_sec_diff_cur = x_ph_first_diff_cur - x_ph_first_diff_prev;
           % Cycle through the spectrum and compute wpd
           for k = 1 : (this.audio_win_length/2) - 1
              wpd = wpd + abs(x_fft(k)*x_ph_sec_diff_cur(k)); 
           end   
           wpd = (1/this.audio_win_length).*wpd;
           
           % Pushing zeros in for interpolation
           upsampled_wpd = [
              wpd
              zeros(this.feature_upsample_factor-1,1)
           ];
        
           % Filter
           [interpolated_wpd, filter_delays] = ...
				filter(this.LP_num, this.LP_denom, upsampled_wpd, filter_delays, 1);
           
           save_indices = ...
              (this.feature_upsample_factor*(n-1)+1):this.feature_upsample_factor*n;
          
           this.frame_wpd(save_indices,:) = interpolated_wpd;
           % Place frame zero crossing rate in feature matrix
           this.feature_matrix(save_indices, :) = interpolated_wpd;
           % Store current phase info
           x_ph_prev = x_ph_curr;
           x_ph_first_diff_prev = x_ph_first_diff_cur;
        end
    end
    
    function plot_sample_intermediate_data(this, sample_frames)
        for frame = sample_frames
            
            if frame == 1
                    warning('Cannot calculate feature differential for frame 1');
                    unrectified_difference = this.frame_wpd(frame, :);
                else
                    unrectified_difference = this.frame_wpd(frame, :) ...
                    - this.frame_wpd(frame-1, :);
            end
            
            HWR_difference = unrectified_difference;
            HWR_difference(HWR_difference < 0) = 0;
            
            diff(frame,1) = HWR_difference;
            
            
        end
        figure;
        subplot(2, 1, 1);
        plot(autocorrelation(diff));
        title(sprintf('weighted phase difference, frame %d', frame));

        subplot(2, 1, 2);
        plot(autocorrelation(this.frame_wpd(sample_frames,1)));
        title(sprintf('Difference, frame %d->%d', frame-1, frame)); 
    end
end

end
