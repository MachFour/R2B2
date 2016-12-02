
% 1 dimensional mean shift clustering, which takes a vector of data, and groups
% points that are close together into clusters, represented by their mean.
% Clusters are returned in the form of a vector of pairs (a 2x1 matrix), where the
% first column is the mean of the cluster, and the second point is how many of the
% original datapoints are represented by that cluster.

% PARAMETERS:
% data = zeros(n, 1)
% 	is the 1 dimensional vector of data points to cluster
% max_dist
% 	is a number representing how far two datapoints in the same cluster can be.
% 	Distance between two datapoints x and y is measured by abs(x - y).

function clusters = mean_shift_cluster(data, max_dist)
	% from https://en.wikipedia.org/wiki/Mean_shift

	n = length(data);

	if n ~= size(data, 1)
		warning('First dimension of data is not the longest');
	end
	clusters = zeros(n, 2);

	% sort the data to make clustering easier
	sorted_data = sort(data, 1);

	means = zeros(n, 1);

	convergence_threshold = 0.01;

	while 1
		for data_idx = 1:n
			datapoint = sorted_data(data_idx);

			% find the set of data points within max_dist
			% since this is in 1 dimension we can just search forward and backward to
			% find the bounds
			high_bound = data_idx;
			low_bound = data_idx;
			while high_bound < n && abs(sorted_data(high_bound) - datapoint) < max_dist
				high_bound = high_bound + 1;
			end
			while low_bound > 1 && abs(datapoint - sorted_data(low_bound)) < max_dist
				low_bound = low_bound - 1;
			end

			neighbourhood = low_bound:high_bound;
			neighbourhood_size = length(neighbourhood) + 1;

			means(data_idx) = sum(sorted_data(neighbourhood))/neighbourhood_size;
		end
		% test for convergence
		if sum(abs(means - sorted_data)^2)/n > convergence_threshold
			break;
		else
			% set x <- m(x)
			sorted_data = means;
			continue;
		end
	end
	clusters = means;
end

