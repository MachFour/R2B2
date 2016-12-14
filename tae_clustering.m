% Class that performs tempo and beat alignment estimation from a single dimensional
% feature, by using autocorrelation to detect likely periodicities (i.e.
% possible tempos), and then for each of those, using a sequence of impulses
% spaced at that tempo period to locate where the periodic peaks of the
% feature actually are (i.e. the beat alignment)

% uses the second type of autocorrelation function which slides
% half the feature frame back over the other half

% estimates for each feature are created via clustering

% Author: Max Fisher

classdef tae_clustering < tempo_alignment_estimator

methods

	function tae = tae_clustering(params, feature_matrix, name)
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

		feature_frame = this.get_feature_frame(frame_number);
		num_features = size(feature_frame, 2);

		feature_tempo_estimates = cell(1, num_features);
		feature_beat_alignment_estimates = cell(1, num_features);

		for n = 1:num_features
			[feature_tempo_estimates{n}, ~] = ...
				this.pick_likely_tempos(feature_frame(:, n));
			
			[feature_beat_alignment_estimates{n}, ~] = ...
				this.pick_likely_beat_alignments(feature_frame(:, n));
		end
		
		% distance function
		alignment_dist = @(x, y) sqrt(sum(abs(x - y).^2));
		alignment_min_sep = 3; % samples
		
		sample_to_bpm_factor = this.params.feature_sample_rate*60;
		bpm_dist = @(x, y) sqrt(sum(sample_to_bpm_factor.*(1./x - 1./y)).^2);
		tempo_min_sep = 4; % BPM
		
		tempo_cluster_dict = tae_clustering.cluster_estimate_lists(...
			feature_tempo_estimates, bpm_dist, tempo_min_sep);
		beat_alignment_cluster_dict = tae_clustering.cluster_estimate_lists(...
			feature_beat_alignment_estimates, alignment_dist, alignment_min_sep);
		
		% make estimate tuples as the cross product of all clustered tempos and
		% all clustered_beat_alignments
		
		% but what if there are no estimates of one or the other type?
		
		frame_tempos = tempo_cluster_dict.keys();
		frame_alignments = beat_alignment_cluster_dict.keys();
		
		num_frame_tempos = length(frame_tempos);
		num_frame_alignments = length(frame_alignments);
		
		frame_estimate_list = zeros(num_frame_tempos*num_frame_alignments, 4);
		
		combined_estimate_number = 0;
		
		for tempo_idx = 1:num_frame_tempos
			tempo_estimate = frame_tempos{tempo_idx};
			tempo_popularity = tempo_cluster_dict(tempo_estimate);
			
			for align_idx = 1:num_frame_alignments
				alignment_estimate = frame_alignments{align_idx};
				alignment_popularity = beat_alignment_cluster_dict(alignment_estimate);
				
				% only allow tempos and beat alignments that more than one
				% feature votes for, and also the alignment estimate has to be
				% within one beat period of the end of the frame
				if tempo_popularity > 1 && alignment_popularity > 1 ...
						&& tempo_estimate > this.params.min_lag_samples/4 ...
						&& tempo_estimate > abs(alignment_estimate)
					combined_estimate_number = combined_estimate_number + 1;
					frame_estimate_list(combined_estimate_number, :) = [
						round(tempo_estimate), tempo_popularity, ...
						round(alignment_estimate), alignment_popularity];
				end
			end
		end
		
		frame_estimate_list = frame_estimate_list(1:combined_estimate_number, :);

		%keep return type the same as before
		frame_estimates = cell(1, 1);
		
		frame_estimates{1} = frame_estimate_list;
		
		% save internally for plotting, data output, etc.
		this.tempo_alignment_estimates{frame_number, 1} = frame_estimate_list;

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
	% peaks of the feature frame. The feature frame must therefore be an onset
	% detection function, so that it peaks at significant events in the audio.
	function [alignments, confidences] = pick_likely_beat_alignments(this, feature_frame)
		if size(feature_frame, 2) ~= 1
			warning('extra columns of feature_frame are ignored');
		end

		[confidences, alignments]  = findpeaks(...
			feature_frame, ...
			0:this.params.feature_win_length -1, ...
			'NPeaks', this.params.max_alignment_peaks, ...
			'SortStr', 'descend');

		% we make sign of the beat alignment negative, to reflect that the offsets
	 	% are calculated backwards from the end of the current feature frame
		alignments = -1*alignments;
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

methods (Static)
	
	% takes in a single dimensional cell array of estimates for each feature, 
	% be they tempo or phase, and turns the whole array into a dictionary whose 
	% keys are the clustered estimates (using the mean shift dict cluster function),
	% and whose values for each key are the indices/lists of the input cell array 
	% that contained an estimate inside that cluster.
	
	% min_sep and dist are as for mean_shift_dict_cluster.m

	function clustered_estimates = cluster_estimate_lists(estimate_lists, dist, min_sep)
		% value is a list of integers, so has to be specified as 'any' type
		unclustered_estimates = containers.Map('KeyType', 'double', 'ValueType', 'any');
		
		for list_idx = 1:length(estimate_lists)
			estimate_list = estimate_lists{list_idx};
			
			for estimate_idx = 1:length(estimate_list)
				estimate = estimate_list(estimate_idx);
				if unclustered_estimates.isKey(estimate)
					% a previous feature (or features) also had this estimate
					other_voters = unclustered_estimates(estimate);
					% append current feature to voter list
					unclustered_estimates(estimate) = [other_voters; list_idx];
				else
					% start a new list of features voting for this estimate
					unclustered_estimates(estimate) = [list_idx];
				end
			end
		end
		% cluster tempo votes
		clustered_estimates = mean_shift_dict_cluster(unclustered_estimates, ...
			dist, min_sep);

		% for each estimate (key) in the clustered estimates dict,
		% replace its value (list of voting features) by just the number of features.
		% Alternatives: sum of confidences? Or actually use which features vote for
		% it?

		cluster_centres = clustered_estimates.keys();
		for cluster_idx = 1:clustered_estimates.length
			cluster_centre = cluster_centres{cluster_idx};
			% replace its value (which is a list) by the length of that list
			num_voters = length(clustered_estimates(cluster_centre));
			clustered_estimates(cluster_centre) = num_voters;
		end
	end
end % methods (Static)

end % classdef
