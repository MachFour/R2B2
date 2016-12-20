%cluster_estimates.m
% takes an array of estimate lists from each feature, for a particular frame of
% features. Estimates of tempos that are close to each other (for different
% features) will be clustered together, then all beat alignments for tempos in each
% cluster are clustered separately. Confidences are discarded.

%Author: Max Fisher


% PARAMETERS:
%	tempo_alignment_estimates = cell(1, num_features)
%		is a cell array as produced by a tempo_alignment estimator, whose elements
%		are a list of estimate tuples for a particular feature frame, in the form
%		(tempo_estimate, tempo_confidence, alignment_estimate, alignment_confidence)
%	tempo_sep, alignment_sep
% 		are positive real numbers specifying the minimum separation between
% 		distinct tempo and beat alignment clusters.

% RETURNS:
% 	clusters = zeros(n, 2)
%		is a single list of clustered tempo and alignment estimates, produced as
%		described above.

function clustered_estimates = cluster_estimates(tempo_alignment_estimates, ...
		max_tempo_peaks, max_alignment_peaks, tempo_sep, alignment_sep)
	num_features = size(tempo_alignment_estimates, 2);

	if num_features == 0
		clustered_estimates = [];
		return
	elseif num_features == 1
		% no clustering SHOULD be needed (assuming peak picking was sensible)
		clustered_estimates = tempo_alignment_estimates{1, 1};
		return
	end

	% add all tempo/alignment estimate pairs to a single list,
	% discard confidence and which feature had which estimate

	% each tempo estimate has its alignment attached
	estimate_pairs = zeros(num_features*max_tempo_peaks*max_alignment_peaks, 2);
	% just tempos only, discarding alignment data
	estimate_tempos = zeros(num_features*max_tempo_peaks, 2);
	num_estimate_pairs = 0;
	num_estimate_tempos = 0;

	for n = 1:num_features
		feature_estimates = tempo_alignment_estimates{n};
		last_estimated_tempo = 0;

		for estimate_idx = 1:size(feature_estimates, 1);
			tempo = feature_estimates(estimate_idx, 1);
			beat_alignment = feature_estimates(estimate_idx, 3);

			num_estimate_pairs = num_estimate_pairs + 1;
			estimate_pairs(num_estimate_pairs, :) = [tempo, beat_alignment];

			if tempo ~= last_estimated_tempo
				num_estimate_tempos = num_estimate_tempos + 1;
				% could use the tempo confidence as the metadata here?
				estimate_tempos(num_estimate_tempos, :) = [tempo, n];
				last_estimated_tempo = tempo;
		end
	end

	estimate_pairs = estimate_pairs(1:num_estimate_pairs, :);
	estimate_tempos = estimate_tempos(1:num_estimate_tempos, :);


	% tempo clusters will be pretty much the same, the only difference will be a
	% slight change in the means due to having multiple identical tempos with
	% different beat alignment, but it should not cause a huge problem, since the
	% number of beat alignment peaks is limited.
	[clustered_tempos, voting_features] = ...
		mean_shift_cluster(estimate_tempos, tempo_sep);
	[~, alignments_for_tempo] = ...
		mean_shift_cluster(estimate_pairs, tempo_sep);

	% alignments_for_tempo is a cell array containing all alignments estimated for
	% the given tempo. Now cluster these, and put into clustered_ta_estimates.

	clustered_estimates = zeros(num_features*max_tempo_peaks*max_alignment_peaks, 4);
	num_clustered_estimates = 0;


	num_tempo_clusters = length(clustered_tempos);
	for tempo_idx = 1:num_tempo_clusters
		alignments = alignments_for_tempo{tempo_idx};
		num_voting_features = length(voting_features{tempo_idx});

		tempo_cluster_mean = clustered_tempos(tempo_idx);
		tempo_cluster_confidence = num_voting_features;


		% give the alignments themselves as the metadata, so we can keep track of
		% the number of alignments per cluster.
		[clustered_alignments, alignments_in_cluster] = ...
			mean_shift_cluster([alignments, alignments], alignment_sep);

		num_alignment_clusters = length(clustered_alignments);
		for align_idx = 1:num_alignment_clusters
			alignment_cluster_mean = clustered_alignments(align_idx);
			alignment_cluster_confidence = length(alignments_in_cluster{align_idx});

			% only include estimate if several features voted for it.
			% round tempo and alignment estimates to nearest sample values.
			if tempo_cluster_confidence > 1 && alignment_cluster_confidence > 1
				num_clustered_estimates = num_clustered_estimates + 1;
				clustered_estimates(num_clustered_estimates, :) = [
					round(tempo_cluster_mean), ...
					tempo_cluster_confidence, ...
					round(alignment_cluster_mean), ...
					alignment_cluster_confidence
					];
			end
		end
	end

	% trim to length
	clustered_estimates = clustered_estimates(1:num_clustered_estimates, :);

end


