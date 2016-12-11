classdef odf_spectralflux < feature_extractor
   properties
      num_feature_channels = 1;
      feature_upsample_factor = 2;
      
      frame_spectral_flux;
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
        
            this.frame_spectral_flux = ...
            zeros(this.feature_upsample_factor*this.num_audio_frames,1);
        
            LP_cutoff_freq_normalised = this.lp_cutoff_freq/(this.feature_sample_rate/2);

            [this.LP_num, this.LP_denom] = butter(this.lp_order, LP_cutoff_freq_normalised);
            this.LP_gain = this.feature_sample_rate/this.lp_cutoff_freq;
       end
       
       function compute_feature(this)
           % Store delays between each iteration
           % Used for LP filter when interpolating
           filter_delays = zeros(this.lp_order, 1);
           % prev_x_fft_mag stores magnitude of once delayed delayed fft
           prev_x_fft_mag = zeros(this.audio_win_length/2,1);
           for n = 1:this.num_audio_frames
                % Fetch current frame
                curr_frame = this.get_windowed_audio_frame(n);
                % Magnitude of current frame
                curr_x_fft_mag = abs(fft(curr_frame, this.audio_win_length));
                curr_x_fft_mag = curr_x_fft_mag(1:this.audio_win_length/2);
                % Take the difference between cur and prev fft over all k
                % and apply HWR function: HWR = ( x - |x| ) / 2
                spec_diff = 0;
                for k = 1: (this.audio_win_length/2) -1   
                    spec_diff = spec_diff + ( (curr_x_fft_mag(k)-prev_x_fft_mag(k)) ...
                        - abs((curr_x_fft_mag(k)-prev_x_fft_mag(k))) )/2;        
                end
                
                % Pushing zeros in for interpolation
                upsampled_spec_diff = [
                  spec_diff
                  zeros(this.feature_upsample_factor-1,1)
                ];
        
                % Filter
                [interpolated_spec_diff, filter_delays] = ...
                    filter(this.LP_num, this.LP_denom, upsampled_spec_diff, filter_delays, 1);

                save_indices = ...
                  (this.feature_upsample_factor*(n-1)+1):this.feature_upsample_factor*n;
                % Store
                this.frame_spectral_flux(save_indices, :) = interpolated_spec_diff;
                this.feature_matrix(save_indices, :) = interpolated_spec_diff;
                % Save cur fft
                prev_x_fft_mag = curr_x_fft_mag;        
           end
       end
           
        function plot_sample_intermediate_data(this, sample_frames)
            for frame = sample_frames
                if frame == 1
                    warning('Cannot calculate feature differential for frame 1');
                    difference = this.frame_spectral_flux(frame, :);
                else
                    difference = this.frame_spectral_flux(frame, :) ...
                        - this.frame_spectral_flux(frame-1, :);
                end
                diff(frame) = difference;   
            end
        
            figure;
            subplot(2, 1, 1);
            plot(autocorrelation(diff));
            title(sprintf('spectral flux, frame %d', frame));

            subplot(2, 1, 2);
            plot(xcorr(this.frame_spectral_flux(sample_frames,1)));
            title(sprintf('Difference, frame %d->%d', frame-1, frame));        
        end
       
   end% Methods    
end