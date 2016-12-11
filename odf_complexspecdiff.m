classdef odf_complexspecdiff < feature_extractor
   properties
      % One channel 
      num_feature_channels = 1;
      % No upsampling performed
      feature_upsample_factor = 2;
      % Stores complex spectral diff of each audio frame
      frame_complex_spectral_diff;
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
       %Init
       function initialise(this, audio_data, audio_sample_rate, name)
           
            initialise@feature_extractor(this, audio_data, audio_sample_rate, name);
            
            this.frame_complex_spectral_diff = ...
            zeros(this.feature_upsample_factor*this.num_audio_frames,1);
        
            LP_cutoff_freq_normalised = this.lp_cutoff_freq/(this.feature_sample_rate/2);

            [this.LP_num, this.LP_denom] = butter(this.lp_order, LP_cutoff_freq_normalised);
            this.LP_gain = this.feature_sample_rate/this.lp_cutoff_freq;
        
            
       end
       %Compute
       function compute_feature(this)
           % Store delays between each iteration
           % Used for LP filter when interpolating
           filter_delays = zeros(this.lp_order, 1);
           % prev_x_fft stores one delayed frame fft 
           prev_x_fft = zeros(this.audio_win_length/2,1);
           % prev_prev_x_fft stores two delayed frame fft
           prev_prev_x_fft = zeros(this.audio_win_length/2,1);
           
           for n = 1:this.num_audio_frames
                % Get audio frame
                curr_frame = this.get_windowed_audio_frame(n);              
                % FFT of current frame
                curr_x_fft = (fft(curr_frame, this.audio_win_length));
                % Consider half of spectrum
                curr_x_fft = curr_x_fft(1:this.audio_win_length/2);
                % Estimate phase from two prev frames (might need to unwrap?)
                curr_x_ph_est = 2.*angle(prev_x_fft)-angle(prev_prev_x_fft);
                % Estimate magnitude of current frame fft from previous
                curr_x_fft_est = abs(prev_x_fft).*curr_x_ph_est;             
                % Compute euclidean distance between estimate and actual
                % Over all k
                csd = 0;
                for k = 1: (this.audio_win_length/2) -1 
                    csd = csd + abs(curr_x_fft(k)-curr_x_fft_est(k));
                end
                % Pushing zeros in for interpolation
                upsampled_csd = [
                    csd
                    zeros(this.feature_upsample_factor-1,1)
                ];
                % Filter
                [interpolated_csd, filter_delays] = ...
                    filter(this.LP_num, this.LP_denom, upsampled_csd, filter_delays, 1);
                save_indices = ...
                  (this.feature_upsample_factor*(n-1)+1):this.feature_upsample_factor*n;               
                % Store
                this.frame_complex_spectral_diff(save_indices, :) = interpolated_csd;
                this.feature_matrix(save_indices, :) = interpolated_csd;
                % Update values
                prev_prev_x_fft = prev_x_fft;
                prev_x_fft = curr_x_fft;
           end       
       end
       
       function plot_sample_intermediate_data(this, sample_frames)
            for frame = sample_frames
                if frame == 1
                    warning('Cannot calculate feature differential for frame 1');
                    difference = this.frame_complex_spectral_diff(frame, :);
                else
                    difference = this.frame_complex_spectral_diff(frame, :) ...
                        - this.frame_complex_spectral_diff(frame-1, :);
                end
                diff(frame) = difference;   
            end
        
            figure;
            subplot(2, 1, 1);
            plot(autocorrelation(diff));
            title(sprintf('complex diff, frame %d', frame));

            subplot(2, 1, 2);
            plot(xcorr(this.frame_complex_spectral_diff(sample_frames,1)));
            title(sprintf('Difference, frame %d->%d', frame-1, frame));        
        end
   end
   
end