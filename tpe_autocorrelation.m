
% Class that performs tempo and phase estimation from a single dimensional
% feature, by using autocorrelation to detect likely periodicities (i.e.
% possible tempos), and then for each of those, using a sequence of impulses
% spaced at that tempo period to locate where the periodic peaks of the
% feature actually are (i.e. the beat phase)

% Author: Max Fisher

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

	% intermediate storage of tempo estimates
	% this is a cell array in the same format as tempo_phase_estimates
	% but each element of this cell array is a list of tuples of the form
	% (tempo, tempo_confidence)
	tempo_estimates;

	% store a measure of the tempo strength, for each candidate tempo,
	% as an alternative to peak picking of the raw autocorrelation function.
	tempo_strength_data;
end


methods
	function initialise(this, feature_data, feature_sample_rate, estimator_name)
		initialise@tempo_phase_estimator(this, feature_data, feature_sample_rate, estimator_name);
		% there needs to be a margin between 4*feature_win_time and max_lag,
		% since at the end of the unbiased autocorrelation function, the amplitude goes
		% up artificially. When using 4 equally spaced impulses at a given lag
		% to measure the periodicity (autocorrelation comb), at higher lags the 4th
		% impulse will pass through this artifically high zone. So we limit max_lag.
		if this.max_lag_samples > this.feature_win_length/4 - 30
			warning('Reduced maximum lag to be 30 samples less than feature_win_length/4');
			this.max_lag_samples = this.feature_win_length/4 - 30;
			fprintf('New min_bpm: %3.2f BPM\n', this.min_bpm);
			fprintf('New max_lag: %3.2f s\n', this.max_lag);

		end

		this.autocorrelation_data = zeros(this.num_feature_frames, this.feature_win_length);
		this.tempo_strength_data = zeros(this.num_feature_frames, length(this.tempo_lag_range));
	end

	function print_properties(this)
		print_properties@tempo_phase_estimator(this);
		fprintf('Maximum number of tempo estimates per feature frame: \t %d\n', ...
			this.MAX_TEMPO_PEAKS);
		fprintf(strcat('Maximum number of phase estimates per tempo estimate, ', ...
			'per feature frame: \t %d\n'), this.MAX_TEMPO_PEAKS);
		end

	function compute_tempo_phase_estimates(this)
		this.compute_tempo_estimates
		this.compute_phase_estimates
	end

	% Performs an autocorrelation to detect the most likely periodicities in each frame,
	% and stores them in the intermediate cell array, tempo_estimates
	function compute_tempo_estimates(this)
		for k = 1:this.num_feature_frames

			% go through and calculate autocorrelations for each slice of
			% FEATURE_LENGTH detection function features
			% no windowing so far

			% it's a column vector
			curr_feature_frame = this.get_feature_frame(k);
			this.autocorrelation_data(k, :) = autocorrelation(curr_feature_frame);

			% now pick peaks!
			% make sure to correct for 1-indexing in arrays
			acf_data_in_tempo_range = this.autocorrelation_data(k, this.tempo_lag_range+1);

			% note that there may be peaks at smaller lag than MAX_BPM that reflect the tatum
			% (semiquaver/quaver) pulse
			% find peaks until confidence is within a certain ratio of the highest?
			% MinPeakProminence: necessary vertical descent on both sides in order to count as
			% a peak
			% MinPeakDistance: horizonal separation between peaks.

			[curr_tempo_confidences, curr_tempo_estimates]  = ...
				findpeaks(acf_data_in_tempo_range,  ...
					this.tempo_lag_range, ...
					'MinPeakProminence', 0.05, ...
					'MinPeakDistance', this.min_lag_samples/2, ...
					'MinPeakHeight', 0, ...
					'NPeaks', this.MAX_TEMPO_PEAKS, ...
					'SortStr', 'descend');
			% remove peak at 0, and peaks with confidence below 0.


			%%%%%%%%%%%%%%%
			% alternative: 'SHIFT-INVARIANT-COMB-FILTERBANK'
			% multiply with periodic comb for each possible tempo in range
			% range of comb separations to try (in samples)
			%%%%%%%%%%%%%%
			for m = 1:length(this.tempo_lag_range)
				lag_candidate = this.tempo_lag_range(m);
				acf_comb = autocorrelation_comb(length(curr_feature_frame), lag_candidate);
				%%% CHECK THIS
				this.tempo_strength_data(k, m) = this.autocorrelation_data(k, :)*acf_comb;
			end

			% intermediate storage of tempo estimates
			% (storing the transpose of the two lists is equivalent to storing each
			% tempo/confidence tuple as a row)
			this.tempo_estimates{k} = [curr_tempo_estimates', curr_tempo_confidences'];
		end
	end

	%%%%%%%%%%%%%%%%%%%%
	% Compute phase estimates for each tempo hypothesis
	% following procedure in 'Context-Dependent Beat tracking'
	%%%%%%%%%%%%%%%%%%%%%
	function compute_phase_estimates(this)
		for k = 1:this.num_feature_frames
			curr_tempo_estimates = this.tempo_estimates{k}(:, 1);
			curr_tempo_confidences = this.tempo_estimates{k}(:, 2);

			% if there were no significant peaks in the autocorrelation, just skip
			% calculating phase alignments
			if isempty(curr_tempo_estimates)
				this.tempo_phase_estimates{k} = [];
				continue;
			end;

			curr_feature_frame = this.get_feature_frame(k);

			% build up rows of estimates for each frame inside the following loop
			% this avoids having 'null estimates' which happens when you
			% preallocate the array with zeros() but then the findpeaks
			% doesn't find as many peaks as the maximum
			% This could also be done by removing all the zero rows after the
			% following loop, but I wasn't bothered

			% frame_estimates = zeros(length(tempo_estimates)*this.MAX_PHASE_PEAKS, 4)
			frame_estimates = zeros(1, 4);
			estimate_number = 0;

			for tempo_index = 1:length(curr_tempo_estimates)
				% in samples
				curr_tempo_estimate = curr_tempo_estimates(tempo_index);
				curr_tempo_confidence = curr_tempo_confidences(tempo_index);

				% make an impulse train for each tempo hypothesis,
				% then slide it along the detection function for one tempo period,
				% to find where it lines up the best
				% feature needs to be a column vector
				
				alignment_function = ...
					beat_alignment_function(curr_feature_frame, curr_tempo_estimate);

				% max distance between peaks depends on the tempo
				% assume that 'off by half a semiquaver' is the same beat
				% location
				[phase_confidences, phase_estimates]  = ...
					findpeaks(alignment_function, ...
						0:curr_tempo_estimate-1, ...
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
					estimate_number = estimate_number + 1;
					tempo_phase_estimate_tuple = [
						curr_tempo_estimate, ...
						curr_tempo_confidence, ...
						phase_estimates(phase_index), ...
						phase_confidences(phase_index)
					];

					frame_estimates(estimate_number, :) = tempo_phase_estimate_tuple;

				end
			end

			this.tempo_phase_estimates{k} = frame_estimates;
		end
	end

	% for each sample frame, plot the detection function,
	% ACF with peaks indicated,
	% and then a plot of superimposed tempo and phase estimates

	% also do a plot of tempo vs phase; plot different tempi on different
	% subplots
	function plot_sample_intermediate_data(this, sample_frames)
		time_axis = this.feature_time_axis;
		for sample_frame = sample_frames
			if sample_frame > this.num_feature_frames
				warning('Cant plot sample frame %d; there are only %d feature frames\n', ...
					sample_frame, this.num_feature_frames);
				continue;
			end

			curr_frame = this.get_feature_frame(sample_frame);

% 			subplot(3, 1, 1);

% 			plot(time_axis, curr_frame);
% 
% 			title(sprintf('Detection function: t=%3.2f s', ...
% 				(sample_frame-1)*this.feature_hop_size/this.feature_sample_rate));
% 			xlabel('Time (seconds)');

			acf_data = this.autocorrelation_data(sample_frame, :);

			figure;

			plot(0:length(acf_data)-1, acf_data);
			hold on;
			title(sprintf('ACF: t=%3.2f s', ...
				(sample_frame-1)*this.feature_hop_size/this.feature_sample_rate));
			xlabel('Lag (samples)');
			ylabel('Similarity');
			stem(this.min_lag_samples, 1);
			stem(this.max_lag_samples, 1);

			curr_tp_estimates = this.tempo_phase_estimates{sample_frame};

			% add indications of peaks picked
			for estimate_idx = 1:size(curr_tp_estimates, 1)
				tempo_estimate = curr_tp_estimates(estimate_idx, 1);
				tempo_confidence = curr_tp_estimates(estimate_idx, 2);
				if tempo_estimate ~= 0
					stem(tempo_estimate, tempo_confidence);
				end
			end

			hold on;

			%plot autocorrelation comb function
			plot(this.tempo_lag_range, this.tempo_strength_data(sample_frame, :));
			%title(sprintf('Tempo strength'));
			%xlabel('Tempo (lag)');

 			figure;

			% plot esimated beat locations as impulses with height equal to
			% the alignment confidence, which can be overlaid on top of the feature
			% vector

			% First we pick only the 'most confident' estimates
			estimates_to_plot = 4;

			% add a column representing some sort of combined confidence of tempo
			% and phase, and then plot the most confident of these

			phase_confidences = curr_tp_estimates(:, 4);

			combined_confidences = phase_confidences;

			% add this column to the estimate array and sort estimate tuples
			% by combined confidence, in descending order
			sorted_estimates = sortrows([curr_tp_estimates, combined_confidences], -5);

			% Here's where we actually calculate and store the beat locations

			estimated_beat_locs = zeros(length(time_axis), estimates_to_plot);

			for estimate_idx = 1:min(estimates_to_plot, size(sorted_estimates, 1))
				tempo_estimate = sorted_estimates(estimate_idx, 1);
				phase_estimate = sorted_estimates(estimate_idx, 3);
				combined_confidence = sorted_estimates(estimate_idx, 5);

				% calculate beat locations by adding multiples of the tempo
				% hypothesis to the beat alignment estimate
				estimated_beat_locs(end+phase_estimate:-tempo_estimate:1, ...
					estimate_idx) = combined_confidence;
			end

			plot(time_axis, estimated_beat_locs);
			hold on; 
			plot(time_axis, curr_frame);
			title('Most confident beat locations');
			xlabel('Time (s)');
			ylabel('Confidence');
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
