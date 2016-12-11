% Class that performs tempo and phase estimation from a single dimensional
% feature, by using autocorrelation to detect likely periodicities (i.e.
% possible tempos), and then for each of those, using a sequence of impulses
% spaced at that tempo period to locate where the periodic peaks of the
% feature actually are (i.e. the beat phase)

% uses the second type of autocorrelation function which slides
% half the feature frame back over the other half

% Author: Max Fisher

classdef tpe_autocorrelation < tempo_phase_estimator

properties
	% intermediate storage of tempo estimates
	% this is a cell array in the same format as tempo_phase_estimates
	% but each element of this cell array is a list of tuples of the form
	% (tempo, tempo_confidence)
	tempo_estimates;
end


methods
	function initialise(this, params, feature_matrix, estimator_name)
		initialise@tempo_phase_estimator(this, params, feature_matrix, estimator_name);
		this.tempo_estimates = cell(this.num_feature_frames, this.num_features);
	end

	function compute_tempo_phase_estimates(this)
		this.compute_tempo_estimates
		this.compute_phase_estimates
	end

	% Performs an autocorrelation to detect the most likely periodicities in each frame,
	% and stores them in the intermediate cell array, tempo_estimates
	function compute_tempo_estimates(this)
		for k = 1:this.num_feature_frames
			tempo_list = this.params.tempo_lag_range';
			% this is a matrix, containing the frames for each feature
			curr_feature_frame = this.get_feature_frame(k);
			% used for the autocorrelation
			prev_nonoverlapping_samples = this.get_prev_nonoverlapping_frame(k);

			for n = 1:this.num_features
				feature_frame_n = curr_feature_frame(:, n);

				% if we have previous data, we can use it to improve the
				% autocorrelation estimates. This isn't a strict autocorrelation
				% any more, but we assume that the previous samples together
				% with the current frame have a constant tempo.
				if ~isempty(prev_nonoverlapping_samples)
					previous_data = prev_nonoverlapping_samples(:, n);
					acf = autocorrelation(feature_frame_n, previous_data);
				else
					% in the beginning when we don't have a previous frame, just
					% use an unbiased autocorrelation
					acf = autocorrelation(feature_frame_n);
				end

				% now pick peaks!
				% make sure to correct for 1-indexing in arrays

				% note that there may be peaks at smaller lag than max_bpm that
				% reflect the tatum (semiquaver/quaver) pulse

				% MinPeakProminence: necessary vertical descent on both sides
				% in order to count as a peak
				% since all values are less than 1, we make this 0.025

				% MinPeakDistance: horizonal separation between peaks
				% -> make it equivalent to 2x the maximum allowed tempo

				% don't sort in descending order: find faster tempos first
				% limit range of allowed tempos to have a minimum of 40BPM

				% heuristic: use peaks at half the given tempo to support the
				% given tempo. Do this by averaging the acf with a compressed
				% version
				doubled_acf = acf;

				for i = 1:length(acf)/2;
					doubled_acf(i) = acf(i)/2 + (acf(2*i - 1) + acf(2*i))/4;
				end

				% could be useful for compound time music?
				%tripled_acf = acf;

				%for i = 1:length(acf)/3;
				%	tripled_acf(i) = acf(i)/2 + ...
				%		(acf(3*i-2) + acf(3*i-1) + acf(3*i))/6;
				%end

				[tempo_confidences, tempo_peaks]  = ...
					findpeaks(doubled_acf(1+tempo_list), ...
						tempo_list, ...
						'MinPeakDistance', this.params.min_lag_samples/2, ...
						'MinPeakProminence', 0.025, ...
						'NPeaks', this.params.max_tempo_peaks, ...
						'SortStr', 'descend');

				% intermediate storage of tempo estimates
				this.tempo_estimates{k, n} = [tempo_peaks, tempo_confidences];
			end
		end
	end

	% Narrows down the search of possible tempos to only a small
	% subset identified by the tempo/phase estimator as likely.
	% Different features will have picked different peaks for their autocorrelation
	% functions, but if the different estimates are close enough, we treat them as
	% both estimating the same tempo.  The idea is that if a lot of different
	% features agree on the tempo, then it's more likely to be that tempo.

	% PARAMETERS:
	% tempo_estimates = cell(1, this.num_features)
	%	  is a cell matrix containing a set of tempo estimates (no confidences)
	%	  for the current frame. There is one set for each feature, which means the cell
	%	  matrix is 1 row (since it's only one frame's worth of data) by n columns, where
	%	  n is the number of features. The tempo estimates should be in SAMPLES.
	function clustered_tempos = cluster_tempo_estimates(this, tempo_estimates)

		% stores the set of tempos to search over in this iteration of the
		% Viterbi algorithm
		% each tempo is tagged with the number of the feature that it comes
		% from.
		candidate_tempos = zeros(this.num_features*this.params.max_tempo_peaks, 2);
		num_candidate_tempos = 0;

		for feature_idx = 1:this.num_features
			tempo_estimates_list = tempo_estimates{feature_idx};
			for estimate_idx = 1:length(tempo_estimates_list)
				num_candidate_tempos = num_candidate_tempos + 1;

				tempo_estimate = tempo_estimates_list(estimate_idx);
				candidate_tempos(num_candidate_tempos, 1) = tempo_estimate;
				candidate_tempos(num_candidate_tempos, 2) = feature_idx;
			end
		end
		% trim list to its actual length
		candidate_tempos = candidate_tempos(1:num_candidate_tempos, :);

		% now group together similar tempo estimates;
		% use mean shift clustering *BY BPM*.
		% 2BPM seems like a reasonable
		% cluster width

		sample_to_bpm_factor = this.params.feature_sample_rate*60;

		bpm_distance = @(x, y) sqrt(sum(sample_to_bpm_factor.*(1./x - 1./y)).^2);
		cluster_width = 4; % BPM
		clustered_tempos = mean_shift_cluster(candidate_tempos, ...
			bpm_distance, cluster_width);
		% what to do if there are no tempos?
	end


	%%%%%%%%%%%%%%%%%%%%
	% Compute phase estimates for each tempo hypothesis
	% following procedure in 'Context-Dependent Beat tracking'
	%%%%%%%%%%%%%%%%%%%%%
	function compute_phase_estimates(this)
		for k = 1:this.num_feature_frames
			curr_feature_frame_matrix = this.get_feature_frame(k);
			for n = 1:this.num_features

				feature_n_estimates = this.tempo_estimates{k, n};

				curr_tempo_estimates = feature_n_estimates(:, 1);
				curr_tempo_confidences = feature_n_estimates(:, 2);

				% if there were no significant peaks in the autocorrelation, just skip
				% calculating phase alignments
				if isempty(curr_tempo_estimates)
					this.tp_estimates{k, n} = [];
					continue
				end

				curr_feature_frame = curr_feature_frame_matrix(:, n);

				max_tempo_peaks = this.params.max_tempo_peaks;
				max_alignment_peaks = this.params.max_phase_peaks;

				frame_n_tp_estimates = zeros(max_tempo_peaks*max_alignment_peaks, 4);

				estimate_number = 0;

				for tempo_index = 1:length(curr_tempo_estimates)
					% in samples
					curr_tempo_estimate = curr_tempo_estimates(tempo_index);
					curr_tempo_confidence = curr_tempo_confidences(tempo_index);

					if curr_tempo_estimate < this.params.min_lag_samples/4
						% that's just too fast, skip this
						continue;
					end

					% make an impulse train for each tempo hypothesis,
					% then slide it along the detection function for one tempo period,
					% to find where it lines up the best
					% feature needs to be a column vector

					alignment_function = ...
						beat_alignment_function(curr_feature_frame, curr_tempo_estimate);

					% max distance between peaks depends on the tempo
					% assume that 'off by half a semiquaver' is the same beat
					% location
					[alignment_confidences, alignment_peaks]  = ...
						findpeaks(alignment_function, ...
							0:curr_tempo_estimate-1, ...
							'MinPeakDistance', curr_tempo_estimate/8, ...
							'NPeaks', max_alignment_peaks, ...
							'SortStr', 'descend');

					% make sign negative to reflect that the offsets are calculated
					% backwards from the end of the frame

					alignment_peaks = -1*alignment_peaks;

					% append tempo_alignment_estimate tuples to the kth array of the
					% tempo_alignment_estimates cell array
					for alignment_index = 1:length(alignment_peaks)
						estimate_number = estimate_number + 1;
						frame_n_tp_estimates(estimate_number, :) = [
							curr_tempo_estimate, ...
							curr_tempo_confidence, ...
							alignment_peaks(alignment_index), ...
							alignment_confidences(alignment_index)
						];
					end
				end

				% trim to size and put it all into the big cell matrix
				frame_n_tp_estimates = frame_n_tp_estimates(1:estimate_number, :);
				this.tp_estimates{k, n} = frame_n_tp_estimates;
			end
		end
	end

	% for each sample frame, plot the detection function,
	% ACF with peaks indicated,
	% and then a plot of superimposed tempo and phase estimates

	% also do a plot of tempo vs phase; plot different tempi on different
	% subplots
	function plot_sample_intermediate_data(this, sample_frames, feature_idx)

		for sample_frame = sample_frames
			time_axis = this.params.feature_frame_start_time(sample_frame) + ...
				this.feature_time_axis;
			if sample_frame > this.num_feature_frames
				warning('Cant plot sample frame %d; there are only %d feature frames\n', ...
					sample_frame, this.num_feature_frames);
				continue;
			end

			curr_frame_matrix = this.get_feature_frame(sample_frame);
			curr_frame = curr_frame_matrix(:, feature_idx);

% 			subplot(3, 1, 1);

% 			plot(time_axis, curr_frame);

% 			title(sprintf('Detection function: t=%3.2f s', ...
% 				(sample_frame-1)*this.feature_hop_size/this.feature_sample_rate));
% 			xlabel('Time (seconds)');

			acf_data = autocorrelation(curr_frame);

			figure;

			plot((0:length(acf_data)-1)/this.params.feature_sample_rate, acf_data);
			hold on;
			title(sprintf('Feature %d ACF: t=%3.2f s', feature_idx, ...
				this.params.frame_end_time(sample_frame)));
			xlabel('Lag (s)');
			ylabel('Similarity');
			stem(this.params.min_lag, 1, 'black', 'x');
			stem(this.params.max_lag, 1, 'black', 'x');

			curr_tp_estimates = this.tp_estimates{sample_frame, frame_idx};

			% add indications of peaks picked
			for estimate_idx = 1:size(curr_tp_estimates, 1)
				tempo_estimate = curr_tp_estimates(estimate_idx, 1)/...
					this.params.feature_sample_rate;

				tempo_confidence = curr_tp_estimates_s(estimate_idx, 2);

				if tempo_estimate ~= 0
					stem(tempo_estimate, tempo_confidence);
				end
			end

 			figure;

			% plot esimated beat locations as impulses with height equal to
			% the alignment confidence, which can be overlaid on top of the feature
			% vector

			% First we pick only the 'most confident' estimates
			estimates_to_plot = 4;

			% sort estimate tuples by phase confidence, in descending order
			sorted_estimates = sortrows(curr_tp_estimates, -4);

			% Here's where we actually calculate and store the beat locations
			% Do it and plot each hypothesised set of beat locations
			% separately
			for estimate_idx = 1:min(estimates_to_plot, size(sorted_estimates, 1))
				subplot(estimates_to_plot, 1, estimate_idx);
				plot(time_axis, curr_frame);
				hold on;

				tempo_estimate = sorted_estimates(estimate_idx, 1);
				phase_estimate = sorted_estimates(estimate_idx, 3);
				confidence = sorted_estimates(estimate_idx, 4);

				% calculate beat locations by adding multiples of the tempo
				% hypothesis to the beat alignment estimate
				estimated_beat_locs = zeros(length(time_axis), estimates_to_plot);
				estimated_beat_locs(end+phase_estimate:-tempo_estimate:1, ...
					estimate_idx) = confidence;

				plot(time_axis, estimated_beat_locs);

				estimated_period = tempo_estimate/this.params.feature_sample_rate;
				estimated_bpm = 60/estimated_period;
				title(sprintf('Beats: %.1f BPM (%.2f s), confidence = %.2f', ...
					estimated_bpm, estimated_period, confidence));
				xlabel('Time (s)');
				ylabel('Confidence');
			end
		end
	end

end
end
