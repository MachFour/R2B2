classdef odf_zero < feature_extractor
    
properties
      
    num_feature_channels = 1;
    
    num_feature_dimensions = 1;
    
    feature_upsample_factor = 2;

    frame_zero_crossing_rate;
   
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

        this.frame_zero_crossing_rate = ...
            zeros(this.feature_upsample_factor*this.num_audio_frames,this.num_feature_dimensions);
       
        LP_cutoff_freq_normalised = this.lp_cutoff_freq/(this.feature_sample_rate/2);

		[this.LP_num, this.LP_denom] = butter(this.lp_order, LP_cutoff_freq_normalised);
		this.LP_gain = this.feature_sample_rate/this.lp_cutoff_freq;

    end
    
    function compute_feature(this)
        % Store delays between each iteration
        % Used for LP filter when interpolating
        filter_delays = zeros(this.lp_order, this.num_feature_dimensions);

        for n = 1:this.num_audio_frames
            %Initialise zero crossing rate
            zero_crossing_rate = 0;
            % Fetch audio frame
            curr_frame = this.get_windowed_audio_frame(n);
            % Finding the amount of times the sign changes in each frame
            % ie zero crossing rate
            for i = 2:(this.audio_win_length - 1)
                if (curr_frame(i - 1) < 0 && curr_frame(i) > 0) || ...
                        (curr_frame(i - 1) > 0 && curr_frame(i) < 0)   
                    zero_crossing_rate = zero_crossing_rate + 1;
                end           
            end
           
            % Pushing zeros in for interpolation
            upsampled_zero_crossing = [
              zero_crossing_rate
              zeros(this.feature_upsample_factor-1,this.num_feature_dimensions)
            ];
            % Filter
            [interpolated_zero_crossing, filter_delays] = ...
				filter(this.LP_num, this.LP_denom, upsampled_zero_crossing, filter_delays, 1);

            save_indices = (this.feature_upsample_factor*(n-1)+1):this.feature_upsample_factor*n;
            % Is this necessary?
            this.frame_zero_crossing_rate(save_indices,:) = interpolated_zero_crossing;
            % Place frame zero crossing rate in feature matrix
            this.feature_matrix(save_indices, :) = interpolated_zero_crossing;
        end         
    end 
    
    function plot_sample_intermediate_data(this, sample_frames)
        for frame = sample_frames
            if frame == 1
                    warning('Cannot calculate feature differential for frame 1');
                    difference = this.frame_zero_crossing_rate(frame, :);
                else
                    difference = this.frame_zero_crossing_rate(frame, :) ...
                    - this.frame_zero_crossing_rate(frame-1, :);
            end
            % No need to rectify as rate will always be >= 0
            
            HWR_difference = difference;
            HWR_difference(HWR_difference < 0) = 0;
            
            diff(frame,1) = HWR_difference;
            
           
        end
        figure;
        subplot(2, 1, 1);
        plot(autocorrelation(diff));
        title(sprintf('zero crossing rate, frame %d', frame));

        subplot(2, 1, 2);
        plot(autocorrelation(this.frame_zero_crossing_rate(sample_frames,1)));
        title(sprintf('Difference, frame %d->%d', frame-1, frame));        
    end % plot_sample_intermediate_data
end
   
end