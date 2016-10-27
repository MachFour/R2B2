
% Class that performs tempo and phase estimation from a single dimensional
% feature, by using autocorrelation to detect likely periodicities (i.e.
% possible tempos), and then for each of those, using a sequence of impulses
% spaced at that tempo period to locate where the periodic peaks of the
% feature actually are (i.e. the beat phase)

classdef tpe_autocorrelation < tempo_phase_estimator

properties (Constant)
	% how many peaks to pick from the autocorrelation function, per frame
	MAX_TEMPO_PEAKS = 4;
	% how many peaks to pick from the beat alignment (phase) function, per frame
	MAX_PHASE_PEAKS = 4;

	% the maximum total number of tempo/phase estimates for each feature
	% frame is given by MAX_TEMPO_PEAKS*MAX_PHASE_PEAKS
end

properties
	autocorrelation_data;

	% store a measure of the tempo strength, for each candidate tempo,
	% as an alternative to peak picking of the raw autocorrelation function.
	tempo_strength_data;
end


methods
	function initialise(this, estimator_name, feature_data, feature_sample_rate)
		initialise@tempo_phase_estimator(this, estimator_name, feature_data, feature_sample_rate);
		% there needs to be a margin between 4*feature_win_time and max_lag,
		% since at the end of the unbiased autocorrelation function, the amplitude goes
		% up artificially. When using 4 equally spaced impulses at a given lag
		% to measure the periodicity (autocorrelation comb), at higher lags the 4th
		% impulse will pass through this artifically high zone. So we limit max_lag.
		if this.max_lag*this.feature_sample_rate > this.feature_win_length/4 - 30
			warning('Reduced maximum lag to be 30 samples less than feature_win_length/4');
			new_max_lag_samples = this.feature_win_length/4 - 30;
			new_max_lag = new_max_lag_samples/this.feature_sample_rate;
			this.min_bpm = 60/new_max_lag;
			fprintf('New min_bpm: %3.2f BPM\n', this.min_bpm);
			fprintf('New max_lag: %3.2f s\n', this.max_lag);

		end

		this.autocorrelation_data = zeros(this.num_feature_frames, this.feature_win_length);
		% initialise tempo and phase estimate array
		for k = 1:this.num_feature_frames
			% each row stores a different tempo/phase estimates,
			% which is a tuple of the form
			% (tempo_estimate, tempo_confidence, phase estimate, phase confidence)
			% the first and third values are in samples, and must be multiplied by the feature
			% rate in order to get a value in seconds.
			% the second and fourth values are unitless real numbers
			this.tempo_phase_estimates{k} = zeros(this.MAX_TEMPO_PEAKS*this.MAX_PHASE_PEAKS, 4);
		end


	end

	function print_properties(this)
		print_properties@tempo_phase_estimator(this);
		fprintf('Maximum number of tempo estimates per feature frame: \t %d\n', ...
			this.MAX_TEMPO_PEAKS);
		fprintf(strcat('Maximum number of phase estimates per tempo estimate, ', ...
			'per feature frame: \t %d\n'), this.MAX_TEMPO_PEAKS);
		end

	function compute_tempo_phase_estimates(this)
		for k = 1:this.num_feature_frames

			% go through and calculate autocorrelations for each slice of
			% FEATURE_LENGTH detection function features
			% no windowing so far
			curr_feature_frame = this.get_feature_frame(k); % no windowing so far

			this.autocorrelation_data(k, :) = autocorrelation(curr_feature_frame);

			% range of lags to consider (in samples)
			candidate_tempo_range = round(this.feature_sample_rate*this.min_lag) : ...
				round(this.feature_sample_rate*this.max_lag);

			% now pick peaks!
			% note that there may be peaks at smaller lag than MAX_BPM that reflect the tatum
			% (semiquaver/quaver) pulse
			% find peaks until confidence is within a certain ratio of the highest  
			% MinPeakProminence: necessary vertical descent on both sides in order to count as
			% a peak
			% MinPeakDistance: horizonal separation between peaks.
			[tempo_confidences, tempo_estimates]  = ...
				findpeaks(this.autocorrelation_data(k, candidate_tempo_range),  ...
					candidate_tempo_range, ...
					'MinPeakProminence', 0.05, ...
					'MinPeakDistance', this.min_lag/2, ...
					'MinPeakHeight', 0, ...
					'NPeaks', this.MAX_TEMPO_PEAKS, ...
					'SortStr', 'descend');
			% remove peak at 0, and peaks with confidence below 0.

			%%%%%%%%%%%%%%%
			% alternative: 'SHIFT-INVARIANT-COMB-FILTERBANK'
			% multiply with periodic comb for each possible tempo in range
			% range of comb separations to try (in samples)
			%%%%%%%%%%%%%%
			tempo_strength = zeros(length(candidate_tempo_range));
			for m = 1:length(tempo_estimates)
				acf_comb = autocorrelation_comb(this.feature_win_length, tempo_estimates(m));
				%%% CHECK THIS
				tempo_strength(m) = this.autocorrelation_data(k, :)*acf_comb;
			end;

			% if there were no significant peaks in the autocorrelation, just skip
			% calculating phase alignments
			if isempty(tempo_estimates)
				continue;
			end;

			%%%%%%%%%%%%%%%%%%%%
			% Compute phase estimates for each tempo hypothesis
			% following procedure in 'Context-Dependent Beat tracking'
			%%%%%%%%%%%%%%%%%%%%%

			for tempo_index = 1:length(tempo_estimates)
				% in samples
				curr_tempo_estimate = tempo_estimates(tempo_index);
				curr_tempo_confidence = tempo_confidences(tempo_index);

				% make an impulse train for each tempo hypothesis,
				% then slide it along the detection function for one tempo period,
				% to find where it lines up the best
				% give feature as row vector, as it is flipped left-to-right
				alignment_function = ...
					beat_alignment_function(curr_feature_frame', curr_tempo_estimate);

				% max distance between peaks depends on the tempo
				% assume that 'off by half a semiquaver' is the same beat
				% location
				[phase_confidences, phase_estimates]  = ...
					findpeaks(alignment_function, 1:curr_tempo_estimate, ...
						'MinPeakHeight', 0.01, ...
						'MinPeakDistance', curr_tempo_estimate/8, ...
						'NPeaks', this.MAX_PHASE_PEAKS, ...
						'SortStr', 'descend');

				% make sign negative to reflect that the offsets are calculated
				% backwards from the end of the frame

				phase_estimates = -1*phase_estimates;
				
				% could calculate and store beat locations here?
				%beat_locations = zeros(this.feature_win_length, this.MAX_PHASE_PEAKS);
				num_phase_peaks = length(phase_confidences);

				% append tempo_phase_estimate tuples to the kth array of the
				% tempo_phase_estimates cell array
				for phase_index = 1:num_phase_peaks
					tempo_phase_estimate_tuple = [
						curr_tempo_estimate;
						curr_tempo_confidence;
						phase_estimates(phase_index);
						phase_confidences(phase_index);
					];

					this.tempo_phase_estimates{k}(1 + phase_index ...
						+ (tempo_index - 1)*this.MAX_PHASE_PEAKS, :) ...
						= tempo_phase_estimate_tuple;
				end
			end

			% remove null estimates (all entries are zero) from the
			% array of estimates for the current frame
		end
	end

	% for each sample frame, plot the detection function,
	% ACF with peaks indicated,
	% and then a plot of superimposed tempo and phase estimates
	
	% also do a plot of tempo vs phase; plot different tempi on different
	% subplots
	function plot_sample_intermediate_data(this, sample_frames)
		time_axis = this.feature_time_axis;

		for k = sample_frames
			figure; 
			subplot(3, 1, 1);
			curr_frame = this.get_feature_frame(k);
			plot(time_axis, curr_frame);

			title(sprintf('Detection function: t=%3.2f s', ...
				(k-1)*this.feature_hop_size/this.feature_sample_rate));
			xlabel('Time (seconds)');

			subplot(3, 1, 2);
			plot(time_axis, this.autocorrelation_data(k, :));
			hold on;
			title(sprintf('ACF: t=%3.2f s', ...
				(k-1)*this.feature_hop_size/this.feature_sample_rate));
			xlabel('Lag (seconds)');
			ylabel('Similarity');
			stem(this.min_lag, 1);
			stem(this.max_lag, 1);

			curr_tp_estimates = this.tempo_phase_estimates{k};

			% add indications of peaks picked
			for estimate_idx = 1:length(curr_tp_estimates)
				tempo_estimate = curr_tp_estimates(estimate_idx, 1);
				tempo_confidence = curr_tp_estimates(estimate_idx, 2);
				if tempo_estimate ~= 0
					stem(time_axis(tempo_estimate), tempo_confidence);
				end
			end

			% plot autocorrelation comb function
			%plot((candidate_tempos/this.feature_sample_rate), tempo_strength(:, c));
			%title(sprintf('Tempo strength'));
			%xlabel('Tempo (BPM)');

			subplot(3, 1, 3);

			% plot esimated beat locations as impulses with height equal to
			% the alignment confidence, which can be overlaid on top of the feature
			% vector

			estimated_beat_locs = zeros(length(time_axis), size(curr_tp_estimates, 1));

			for estimate_idx = 1:length(curr_tp_estimates)
				tempo_estimate = curr_tp_estimates(estimate_idx, 1);
				tempo_confidence = curr_tp_estimates(estimate_idx, 2);
				phase_estimate = curr_tp_estimates(estimate_idx, 3);
				phase_confidence = curr_tp_estimates(estimate_idx, 4);

				% this sort of makes sense, since the phase estimates are
				% conditional on the tempo estimates
				combined_confidence = tempo_confidence*(1 + phase_confidence);
				if tempo_estimate ~= 0 && phase_confidence > 0
					% calculate beat locations by adding multiples of the tempo
					% hypothesis to the beat alignment estimate
					estimated_beat_locs(end+phase_estimate:-tempo_estimate:1, ...
						estimate_idx) = combined_confidence;
				end
			end

			plot(time_axis, estimated_beat_locs);
			%title(sprintf('tempo = %2.3f, t = %3.2f s, band=%d', ...
			%	tempo_hypothesis/FEATURE_RATE, ...
			%	(k-1)*FEATURE_HOP_SIZE/FEATURE_RATE, ...
			%	c)); 
			%xlabel('Offset from end of frame (seconds)');
			%ylabel('Strength');
		end
	end

end
end
