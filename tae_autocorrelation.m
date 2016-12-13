% Class that performs tempo and beat alignment estimation from a single dimensional
% feature, by using autocorrelation to detect likely periodicities (i.e.
% possible tempos), and then for each of those, using a sequence of impulses
% spaced at that tempo period to locate where the periodic peaks of the
% feature actually are (i.e. the beat alignment)

% uses the second type of autocorrelation function which slides
% half the feature frame back over the other half

% Author: Max Fisher

classdef tae_autocorrelation < tempo_alignment_estimator

methods

	function tae = tae_autocorrelation(params, feature_matrix, name)
		% superclass constructor
		tae = tae@tempo_alignment_estimator(params, feature_matrix, name);
	end


	% comes up with a set of tempo/beat alignment estimates and confidences for each
	% input feature, in the given feature frame.
	% Rows of the feature frame represent feature samples calculated over time, while
	% columns represent different features.
	function frame_estimates = pick_tempo_and_alignment_estimates(this, frame_number)
		if ~(frame_number >= 1)
			error('frame number must be positive');
		end

		% frames of feature_win_length samples.
		% The last sample of past_frame immediately preceeds the first sample of
		% feature_frame. This is used to improve the quality of the autocorrelation.
		feature_frame = this.get_feature_frame(frame_number);
		%past_frame = this.get_prev_nonoverlapping_frame(frame_number);
		%
		%if ~isequal(size(feature_frame), size(past_frame))
		%	error('feature_frame and past_frame must have the same size');
		%end

		% An equivalent operation is just to take a longer feature frame size
		% and split it in half.


		num_features = size(feature_frame, 2);

		% holds lists of tempo/beat alignment estimates for each feature, in this frame
		% the maximum length of each list is max_estimates
		max_estimates_per_feature = ...
			this.params.max_tempo_peaks*this.params.max_alignment_peaks;

		frame_estimates = cell(1, num_features);

		% for each feature, find the max_tempo_peaks most likely tempos (peaks of
		% autocorrelation), and then for each of those tempos, find the
		% max_alignment_peaks most likely beat alignments (peaks of beat alignment
		% function).

		for n = 1:num_features
			% holds all tempo/beat alignment estimate tuples for this feature and frame
			feature_n_estimate_tuples = zeros(max_estimates_per_feature, 4);
			estimate_number = 0;

			[tempo_estimates, tempo_confidences] = ...
				this.pick_likely_tempos(feature_frame(:, n));

			for tempo_idx = 1:length(tempo_estimates)
				[align_estimates, align_confidences] = this.pick_likely_beat_alignments(...
					feature_frame(:, n), tempo_estimates(tempo_idx));

				for alignment_idx = 1:length(align_estimates)
					estimate_number = estimate_number + 1;
					feature_n_estimate_tuples(estimate_number, :) = [
						tempo_estimates(tempo_idx), ...
						tempo_confidences(tempo_idx), ...
						align_estimates(alignment_idx), ...
						align_confidences(alignment_idx)];
				end
			end

			% trim to size
			feature_n_estimate_tuples = feature_n_estimate_tuples(1:estimate_number, :);
			frame_estimates{n} = feature_n_estimate_tuples;

			% save internally for plotting, data output, etc.
			this.tempo_alignment_estimates{frame_number, n} = feature_n_estimate_tuples;

		end

	end

	% Given a feature frame for a single feature, determine which are the most likely
	% tempos in that frame, using an autocorrelation of the current feature frame
	% autocorrelation.m for more details.

	% Returns a list of (tempo, confidence) tuples, where the tempo is measured in
	% samples of ACF lag, and the confidence is the height of the autocorrelation
	% function at that lag.

	function [tempos, confidences] = pick_likely_tempos(this, feature_frame)
		if size(feature_frame, 2) ~= 1
			warning('extra columns of feature_frame are ignored');
		end

		% if we have previous data, we can use it to improve the
		% autocorrelation estimates. This isn't a strict autocorrelation
		% any more, but we assume that the previous samples together
		% with the current frame have a constant tempo.
		% see autocorrelation.m for more information

		% Equivalently, we can use a longer feature frame size and just split it
		% into halves.
		acf = autocorrelation(feature_frame);

		% now pick peaks! (making sure to correct for 1-indexing in arrays)

		% note that there may be peaks at smaller lag than max_bpm that
		% reflect the tatum (semiquaver/quaver) pulse, and peaks at larger lag than
		% min_bpm that reflect bar-length periodicities. It would be nice to
		% integrate these somehow into the algorithm.

		% HEURISTIC: use peaks at half the given tempo to support the
		% given tempo. Do this by averaging the acf with a compressed
		% version
		%doubled_acf = acf;

		%for i = 1:length(acf)/2;
		%	doubled_acf(i) = acf(i)/2 + (acf(2*i - 1) + acf(2*i))/4;
		%end

		% could be useful for compound time music?
		%tripled_acf = acf;

		%for i = 1:length(acf)/3;
		%	tripled_acf(i) = acf(i)/2 + ...
		%		(acf(3*i-2) + acf(3*i-1) + acf(3*i))/6;
		%end


		% MinPeakProminence: necessary vertical descent on both sides
		% in order to count as a peak
		% since all values are less than 1, we make this 0.025

		% MinPeakDistance: horizonal separation between peaks
		% -> make it equivalent to 2x the maximum allowed tempo

		allowed_tempos = this.params.tempo_lag_range';
		[confidences, tempos] = findpeaks(...
			acf(1+allowed_tempos), ...
			allowed_tempos, ...
			'MinPeakDistance', this.params.min_lag_samples/2, ...
			'MinPeakProminence', 0.025, ...
			'MinPeakHeight', 0.2, ...
			'NPeaks', this.params.max_tempo_peaks, ...
			'SortStr', 'descend');
	end


	% For the given tempo estimate, determine likely beat alignments, by finding the
	% shift for which a set of impulses spaced at that tempo estimate best correlates
	% with the given feature frame.
	% This vaguely follows the procedure from 'Context-Dependent Beat tracking'
	function [alignments, confidences] = pick_likely_beat_alignments(this, ...
			feature_frame, tempo_estimate)

		if size(feature_frame, 2) ~= 1
			warning('extra columns of feature_frame are ignored');
		end

		% make an impulse train for each tempo hypothesis,
		% then slide it along the detection function for one tempo period,
		% to find where it lines up the best
		% feature needs to be a column vector

		alignment_function = beat_alignment_function(feature_frame, tempo_estimate);

		% max distance between peaks depends on the tempo: we're making the
		% assumption that we don't need to distnguish between beat locations closer
		% than half a semiquaver at the given tempo estimate


		[confidences, alignments]  = findpeaks(...
			alignment_function, ...
			0:(tempo_estimate-1), ...
			'MinPeakDistance', tempo_estimate/8, ...
			'NPeaks', this.params.max_alignment_peaks, ...
			'SortStr', 'descend');

		% we make sign of the beat alignment negative, to reflect that the offsets
	 	% are calculated backwards from the end of the current feature frame
		alignments = -1*alignments;

	end


	% Narrows down the search of possible tempos to only a small
	% subset identified by the tempo/alignment estimator as likely.
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

	% for each sample frame, plot the detection function,
	% ACF with peaks indicated,
	% and then a plot of superimposed tempo and beat alignment estimates

	% also do a plot of tempo vs beat alignment; plot different tempi on different
	% subplots

	function plot_frame_data(this, frame_number, feature_number)

		time_axis = this.params.feature_frame_start_time(frame_number) + ...
			this.feature_time_axis;
		if frame_number > this.num_feature_frames
			warning('Cant plot sample frame %d; there are only %d feature frames\n', ...
				frame_number, this.num_feature_frames);
			return;
		end

		curr_frame_matrix = this.get_feature_frame(frame_number);
		curr_frame = curr_frame_matrix(:, feature_number);
		past_frame_matrix = this.get_prev_nonoverlapping_frame(frame_number);
		past_frame = past_frame_matrix(:, feature_number);

% 		subplot(3, 1, 1);

% 		plot(time_axis, curr_frame);

% 		title(sprintf('Detection function: t=%3.2f s', ...
% 			(frame_number-1)*this.feature_hop_size/this.feature_sample_rate));
% 		xlabel('Time (seconds)');

		acf_data = autocorrelation(curr_frame, past_frame);

		figure;

		plot((0:length(acf_data)-1)/this.params.feature_sample_rate, acf_data);
		hold on;
		title(sprintf('Feature %d ACF: t=%3.2f s', feature_number, ...
			this.params.frame_end_time(frame_number)));
		xlabel('Lag (s)');
		ylabel('Similarity');
		stem(this.params.min_lag, 1, 'black', 'x');
		stem(this.params.max_lag, 1, 'black', 'x');

		estimates = this.tempo_alignment_estimates{frame_number, feature_number};

		% add indications of peaks picked
		for estimate_idx = 1:size(estimates, 1)
			tempo_estimate = estimates(estimate_idx, 1)/this.params.feature_sample_rate;

			tempo_confidence = estimates(estimate_idx, 2);

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

		% sort estimate tuples by alignment confidence, in descending order
		sorted_estimates = sortrows(estimates, -4);

		% Here's where we actually calculate and store the beat locations
		% Do it and plot each hypothesised set of beat locations
		% separately
		for estimate_idx = 1:min(estimates_to_plot, size(sorted_estimates, 1))
			subplot(estimates_to_plot, 1, estimate_idx);
			plot(time_axis, curr_frame);
			hold on;

			tempo_estimate = sorted_estimates(estimate_idx, 1);
			alignment_estimate = sorted_estimates(estimate_idx, 3);
			confidence = sorted_estimates(estimate_idx, 4);

			% calculate beat locations by adding multiples of the tempo
			% hypothesis to the beat alignment estimate
			estimated_beat_locs = zeros(length(time_axis), estimates_to_plot);
			estimated_beat_locs(end+alignment_estimate:-tempo_estimate:1, ...
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

end % methods

end % classdef
