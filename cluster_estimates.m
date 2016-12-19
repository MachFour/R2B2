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

	max_estimates = num_features*max_tempo_peaks*max_alignment_peaks;
	estimate_pairs = zeros(max_estimates, 2);
	num_estimates = 0;

	for n = 1:num_features
		feature_estimates = tempo_alignment_estimates{n};
		for estimate_idx = 1:size(feature_estimates, 1);
			tempo_estimate = feature_estimates(estimate_idx, 1);
			beat_alignment_estimate = feature_estimates(estimate_idx, 3);

			num_estimates = num_estimates + 1;
			estimate_pairs(num_estimates, :) = ...
				[tempo_estimate, beat_alignment_estimate];

		end
	end

	estimate_pairs = estimate_pairs(1:num_estimates, :);

	[clustered_tempos, alignments_for_tempo] = ...
		mean_shift_cluster(estimate_pairs, tempo_sep);

	% alignments_for_tempo is a cell array containing all alignments estimated for
	% the given tempo. Now cluster these, and put into clustered_ta_estimates.

	clustered_estimates = zeros(max_estimates, 4);
	num_clustered_estimates = 0;


	num_tempo_clusters = length(clustered_tempos);
	for tempo_idx = 1:num_tempo_clusters
		alignments = alignments_for_tempo{tempo_idx};
		% number of features that chose this tempo is approximately
		% round(length(alignments)/params.max_alignment_peaks).

		tempo_cluster_mean = clustered_tempos(tempo_idx);
		tempo_cluster_confidence = round(length(alignments)/max_alignment_peaks);


		% give the alignments themselves as the metadata, so we can keep track of
		% the number of alignments per cluster.
		[clustered_alignments, alignments_in_cluster] = ...
			mean_shift_cluster([alignments, alignments], alignment_sep);

		num_alignment_clusters = length(clustered_alignments);
		for align_idx = 1:num_alignment_clusters
			alignment_cluster_mean = clustered_alignments(align_idx);
			alignment_cluster_confidence = length(alignments_in_cluster{align_idx});

			num_clustered_estimates = num_clustered_estimates + 1;
			clustered_estimates(num_clustered_estimates, :) = [
				tempo_cluster_mean, ...
				tempo_cluster_confidence, ...
				alignment_cluster_mean, ...
				alignment_cluster_confidence
				];
		end
	end
	
	% trim and do rounding to nearest sample

	clustered_estimates = round(clustered_estimates(1:num_clustered_estimates, :));

end


