
% 1 dimensional mean shift clustering, which groups points that are close together
% into clusters, represented by their arithmetic mean.
% Clusters are returned in the form of a 2x1 cell matrix, where the
% first column holds the mean of the cluster, and the second point holds a list
% of the original data rows that are contained in that cluster.

% from https://en.wikipedia.org/wiki/Mean_shift

% PARAMETERS:
% data = matrix(n, k)
% 	is the 1 dimensional cell array of data points to cluster. The extra
% 	k-1 columns can contain other information that is kept along with the data.
% dist = @(x, y)
%	is a function handle describing how to measure the distance
%	between the values x = data(a, 1) and y = data(b, 1) of elements a and b
%	of the data array. It should return a nonnegative real number.
%	e.g. dist = @(x, y) sqrt(sum(abs(x - y).^2));
% max_dist
% 	is a nonnegative real number representing how far two datapoints in the same
%	cluster can be.


function clusters = mean_shift_cluster(data, dist, max_dist)
	n = size(data, 1);
	k = size(data, 2);

	if n == 1
		clusters = {data{1}, data(1, 1:k)};
		return;
	end

	% initialise
	original_sorted_data = sortrows(data, 1);
	cluster_data = original_sorted_data;

	% iterate until convergence (do while loop)
	eps = 1e-9;
	while 1
		% set x <- m(x) (see wiki page)
		old_cluster_data = cluster_data;
		for i = 1:n
			x_i = old_cluster_data(i, 1);

			% find the set of data points within max_dist
			% since this is in 1 dimension we can just search forward and backward to
			% find the bounds
			high_idx = i;
			low_idx = i;
			% search for highest and lowest index bounds such that the
			% distance between the data point and the value at that index
			% is less than max_dist
			while high_idx < n && dist(old_cluster_data(high_idx + 1, 1), x_i) < max_dist
				high_idx = high_idx + 1;
			end
			while low_idx > 1 && dist(old_cluster_data(low_idx - 1, 1), x_i) < max_dist
				low_idx = low_idx - 1;
			end

			nearby_data = low_idx:high_idx;

			average_value = sum(old_cluster_data(nearby_data, 1))/length(nearby_data);
			cluster_data(i, 1) = average_value;

			% keep the extra information columns
			cluster_data(i, 2:k) = old_cluster_data(i, 2:k);
		end
		% test for convergence
		if dist(cluster_data(1:n, 1), old_cluster_data(1:n, 1))/n < eps
			break;
		end
	end

	% now create a cell array where the first column holds the mean of each
	% cluster, and the second column holds a list of the original data rows
	% in that cluster.
	% Problem: not sure how to preallocate array size for rows in the second column.

	clusters = cell(n, 2);
	% the following variables need to keep their values across loop
	% iterations:
	% 	current_cluster_idx - index into rows of clusters cell array
	% 	current_cluster_size - index into clusters{current_cluster_idx, 2} array
	% 	current_cluster_mean - what the value of the current cluster is

	for data_index = 1:n
		% which cluster are we adding to (represented by the mean)
		current_data_row = original_sorted_data(data_index, 1:k);
		% which cluster is this row of data in? (represented by the mean)
		data_cluster = cluster_data(data_index, 1);

		if data_index == 1 || dist(data_cluster, current_cluster_mean) > eps
			% make a new cluster, with a size of 1
			if data_index == 1
				current_cluster_idx = 1;
			else
				% increment cluster
				current_cluster_idx = current_cluster_idx + 1;
			end
			current_cluster_size = 1;
			current_cluster_mean = cluster_data(data_index, 1);

			clusters{current_cluster_idx, 1} = current_cluster_mean;
			clusters{current_cluster_idx, 2} = current_data_row;
		else
			% it's a part of the existing cluster, add the data row
			current_cluster_size = current_cluster_size + 1;
			% get the existing rows and add one
			current_cluster_data_rows = clusters{current_cluster_idx, 2};
			current_cluster_data_rows(current_cluster_size, 1:k) = current_data_row;
			% put them back
			clusters{current_cluster_idx, 2} = current_cluster_data_rows;
		end
	end
	% trim zero rows
	clusters = clusters(1:current_cluster_idx, 1:2);
end

