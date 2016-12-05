data = [
1.0000,    1.0000
1.1000,    2.0000
1.0500,    3.0000
3.0000,    4.0000
3.0000,    5.0000
1.5000,    6.0000
10.0000,    7.0000
10.1000,    8.0000
10.2000,    9.0000
10.3000,   10.0000
20.0000,   11.0000
30.0000,   12.0000
];

desired_output = {
1.05, [1, 1; 1.05, 3; 1.1, 2]
1.5,  [1.5, 6]
3,    [3,4;3,5]
10.1  [10, 7; 10.1, 8; 10.2, 9]
10.3, [10.3, 10]
20,   [20,11]
30,   [30,12]
};

dist = @(x, y) sqrt(sum(abs(x - y).^2));
max_dist = 0.1;
clusters = mean_shift_cluster(data, dist, max_dist);

if isequal(clusters, desired_output)
	disp('Test passed');
else
	disp('Test failed');
end
