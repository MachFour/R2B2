
% 1 dimensional mean shift clustering, which groups points that are close together
% into clusters, represented by their arithmetic mean.
% Clusters are returned in the form of a 2x1 cell matrix, where the
% first column holds the mean of the cluster, and the second point holds a
% list of the metadata for each of the original datapoints contained in that cluster.
% Distance between clusters is measured using the standard euclidean distance metric

% from https://en.wikipedia.org/wiki/Mean_shift

% PARAMETERS:
% data = matrix(n, 2)
% 	is the 1 dimensional array of data points to cluster, together with an extra
% 	column of 'metadata' about the data point.
% min_separation
% 	is a nonnegative real number representing how close two datapoints distinct
% 	clusters can be to each other


function [cluster_means, cluster_metadata] = mean_shift_cluster(data, min_separation)
	if ~isequal(size(data, 2), 2)
		error('data is of incorrect size');
	end

	n = size(data, 1);

	if n == 0
		cluster_means = [];
		cluster_metadata = {};
		return
	end

	if n == 1
		cluster_means = data(1, 1);
		cluster_metadata = {data(1, 2)};
		return
	end

	% initialise
	original_sorted_data = sortrows(data, 1);
	% copy data
	cluster_data = original_sorted_data;

	% iterate until convergence (do while loop)
	eps = 1e-9;
	while 1
		% set x <- m(x) (see wiki page)
		old_cluster_data = cluster_data;
		for i = 1:n
			x_i = old_cluster_data(i, 1);

			% find the set of data points within min_separation
			% since this is in 1 dimension we can just search forward and backward to
			% find the bounds
			high_idx = i;
			low_idx = i;
			% search for highest and lowest index bounds such that the
			% distance between the data point and the value at that index
			% is less than min_separation
			while high_idx < n && distance(old_cluster_data(high_idx + 1, 1), x_i) <= min_separation - eps
				high_idx = high_idx + 1;
			end
			while low_idx > 1 && distance(old_cluster_data(low_idx - 1, 1), x_i) <= min_separation - eps
				low_idx = low_idx - 1;
			end

			nearby_data = low_idx:high_idx;

			average_value = sum(old_cluster_data(nearby_data, 1))/length(nearby_data);
			cluster_data(i, 1) = average_value;

			% keep the metadata
			cluster_data(i, 2) = old_cluster_data(i, 2);
		end
		% test for convergence
		if distance(cluster_data(:, 1), old_cluster_data(:, 1))/n < eps
			break
		end
	end

	% now create an array of the cluster means, and a cell array to hold the list of
	% metadata values for each of the datapoints in the respective clusters.

	cluster_means = zeros(n, 1);
	cluster_metadata = cell(n, 1);
	num_clusters = 0;

	i = 1;
	while i <= n
		current_cluster_mean = cluster_data(i, 1);
		current_cluster_metadata = zeros(n, 1);
		current_cluster_size = 0;

		% for as long as we can increment i and the key cluster stays the same, we're
		% still in the same cluster

		% this loop will always run once, since i = j initially.

		for j = i:n
			cluster_j_mean = cluster_data(j);
			if distance(cluster_j_mean, current_cluster_mean) <= min_separation - eps
				% same cluster, add jth key's value to current cluster values
				current_cluster_size = current_cluster_size + 1;
				current_cluster_metadata(current_cluster_size) = ...
					original_sorted_data(j, 2);
				% update start index of next outer loop
				i = i + 1;
			else
				% we're out of the cluster (since cluster_data is sorted)
				break
			end
		end
		
		% Trim metadata list to size
		current_cluster_metadata = current_cluster_metadata(1:current_cluster_size);

		% make sure to save the cluster!
		num_clusters = num_clusters + 1;
		cluster_means(num_clusters) = current_cluster_mean;
		cluster_metadata{num_clusters} = current_cluster_metadata;
	end

	% trim to correct size
	cluster_means = cluster_means(1:num_clusters);
	cluster_metadata = cluster_metadata(1:num_clusters);

end

function d = distance(x, y)
	d = sqrt(sum((x - y).^2));
end
