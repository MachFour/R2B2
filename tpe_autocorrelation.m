
% Class that performs tempo and phase estimation from a single dimensional
% feature, by using autocorrelation to detect likely periodicities (i.e.
% possible tempos), and then for each of those, using a sequence of impulses
% spaced at that tempo period to locate where the periodic peaks of the
% feature actually are (i.e. the beat phase)

% uses the second type of autocorrelation function which slides
% half the feature frame back over the other half

% Author: Max Fisher

classdef tpe_autocorrelation < tempo_phase_estimator

properties (Constant)
	% how many peaks to pick from the autocorrelation function, per frame
	MAX_TEMPO_PEAKS = 8;
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
end


methods
	function initialise(this, feature_data, feature_sample_rate, estimator_name)
		initialise@tempo_phase_estimator(this, feature_data, feature_sample_rate, estimator_name);
		this.autocorrelation_data = zeros(this.num_feature_frames, this.feature_win_length/2);
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
			% this.feature_win_length

			% we use the new method of autocorrelation

			% it's a column vector
			curr_feature_frame = this.get_feature_frame(k);
			first_half = curr_feature_frame(1:this.feature_win_length/2);
			second_half = curr_feature_frame(this.feature_win_length/2+1:end);

			acf2 = autocorrelation(second_half, first_half);
			this.autocorrelation_data(k, :) = acf2;

			% now pick peaks!
			% make sure to correct for 1-indexing in arrays

			% note that there may be peaks at smaller lag than MAX_BPM that reflect the tatum
			% (semiquaver/quaver) pulse

			% MinPeakProminence: necessary vertical descent on both sides 
			% in order to count as a peak
			% since all values are less than 1, we make this 0.025

			% MinPeakDistance: horizonal separation between peaks
			% -> make it equivalent to 2x the maximum allowed tempo

			% don't sort in descending order: find faster tempos first
			[curr_tempo_confidences, curr_tempo_estimates]  = ...
				findpeaks(this.autocorrelation_data(k, :),  ...
					0:this.feature_win_length/2-1, ...
					'MinPeakDistance', this.min_lag_samples/2, ...
					'MinPeakProminence', 0.025, ...
					'NPeaks', this.MAX_TEMPO_PEAKS, ...
					'SortStr', 'descend');

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

			frame_estimates = zeros(1, 4);
			% same but in seconds
			frame_estimates_s = zeros(1, 4);
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
					tp_estimate_tuple = [
						curr_tempo_estimate, ...
						curr_tempo_confidence, ...
						phase_estimates(phase_index), ...
						phase_confidences(phase_index)
					];

					tp_estimate_tuple_s = [
						curr_tempo_estimate/this.feature_sample_rate, ...
						curr_tempo_confidence, ...
						phase_estimates(phase_index)/this.feature_sample_rate, ...
						phase_confidences(phase_index)
					];

					frame_estimates(estimate_number, :) = tp_estimate_tuple;
					frame_estimates_s(estimate_number, :) = tp_estimate_tuple_s;

				end
			end

			this.tempo_phase_estimates{k} = frame_estimates;
			this.tempo_phase_estimates_s{k} = frame_estimates_s;
		end
	end

	% for each sample frame, plot the detection function,
	% ACF with peaks indicated,
	% and then a plot of superimposed tempo and phase estimates

	% also do a plot of tempo vs phase; plot different tempi on different
	% subplots
	function plot_sample_intermediate_data(this, sample_frames)

		for sample_frame = sample_frames
			time_axis = this.get_frame_end_time(sample_frame) + this.feature_time_axis;
			if sample_frame > this.num_feature_frames
				warning('Cant plot sample frame %d; there are only %d feature frames\n', ...
					sample_frame, this.num_feature_frames);
				continue;
			end

			curr_frame = this.get_feature_frame(sample_frame);

% 			subplot(3, 1, 1);

% 			plot(time_axis, curr_frame);

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

				estimated_period = tempo_estimate/this.feature_sample_rate;
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
