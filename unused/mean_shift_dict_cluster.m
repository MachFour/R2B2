
% 1 dimensional mean shift clustering, which groups points that are close together
% into clusters, represented by their arithmetic mean.
% This function accepts inputs in the form of key/value map, where the keys are the
% data to be clustered, and the value is some sort of identifier. When keys are
% clustered together, they are replaced in the map by a key equal to the cluster
% centre, and whose value is an array formed by the concatenation of the values of
% the original keys. This way, the original data can be identified in clusters.

% from https://en.wikipedia.org/wiki/Mean_shift

% PARAMETERS:
% data = containers.Map('KeyType', 'double', 'ValueType', 'any')
% 	contains the data to be clustered in the form of keys, whose values can be used
% 	to identify the original data when it is clustered. In the indended application,
% 	'ValueType' is a vector of integers. That way the limitation of having unique keys
% 	in the map can be avoided by just having that key's value being an array of
% 	identifiers.
% dist = @(x, y)
%	is a function handle describing how to measure the distance between keys x and y
%	in the data dictionary. It should return a nonnegative real number.
%	e.g. dist = @(x, y) sqrt(sum(abs(x - y).^2));
% min_separation
% 	is a nonnegative real number representing how close points in distinct clusters
% 	can be.


function clusters = mean_shift_dict_cluster(data, distance, min_separation)
	if ~isa(data, 'containers.Map')
		error('Data to cluster must be a map');
	end

	keys = data.keys;
	n = length(keys);

	if n <= 1
		% nothing to cluster; copy input to output.
		clusters = containers.Map(data.keys, data.values);
		return;
	end

	% get keys into a column vector
	unsorted_keys = zeros(n, 1);
	for key_idx = 1:n
		unsorted_keys(key_idx) = keys{key_idx};
	end

	% initialise
	original_sorted_keys = sort(unsorted_keys);
	cluster_data = original_sorted_keys;

	% iterate until convergence (do while loop)
	eps = 1e-9;
	while 1
		% set x <- m(x) (see wiki page)
		old_cluster_data = cluster_data;
		for i = 1:n
			x_i = old_cluster_data(i);

			% find the set of data points within min_separation
			% since data is 1 dimensional and sorted, we can just search forward and backward to
			% find the bounds
			high_idx = i;
			low_idx = i;
			% search for highest and lowest index bounds such that the
			% distance between the data point and the value at that index
			% is less than min_separation
			while high_idx < n && distance(old_cluster_data(high_idx + 1), x_i) <= min_separation - eps
				high_idx = high_idx + 1;
			end
			while low_idx > 1 && distance(old_cluster_data(low_idx - 1), x_i) <= min_separation - eps
				low_idx = low_idx - 1;
			end

			nearby_data = low_idx:high_idx;

			average_value = sum(old_cluster_data(nearby_data, 1))/length(nearby_data);
			cluster_data(i) = average_value;
		end
		% test for convergence
		if distance(cluster_data(:), old_cluster_data(:))/n < eps
			break;
		end
	end

	% now create a map where the keys are the mean of each cluster, and the values
	% are the contatenation of the values of each original key making up the cluster

	clusters = containers.Map('KeyType', 'double', 'ValueType', 'any');
	% the following variables need to keep their values across loop
	% iterations:
	% 	current_cluster_idx - index into rows of clusters cell array
	% 	current_cluster_size - index into clusters{current_cluster_idx, 2} array
	% 	current_cluster_mean - what the value of the current cluster is

	i = 1;
	while i <= n
		cluster_i_mean = cluster_data(i);
		values_in_cluster = [];

		% for as long as we can increment i and the key cluster stays the same, we're
		% still in the same cluster

		% this loop will always run once, since i = j initially.
		for j = i:n
			jth_original_key = original_sorted_keys(j);
			cluster_j_mean = cluster_data(j);
			if distance(cluster_j_mean, cluster_i_mean) <= min_separation - eps
				% same cluster, add jth key's value to cluster values
				values_in_cluster = [values_in_cluster; data(jth_original_key)];
			else
				% we're out of the cluster (since cluster_data is sorted)
				i = j - 1;
				break;
			end
		end
		% here, j is equal to the index of the last successfully processed key 
		% now move to the next cluster
		i = i + 1;

		% make sure to save the cluster!
		clusters(cluster_i_mean) = values_in_cluster;

	end

end

