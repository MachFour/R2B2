data = containers.Map('keyType', 'double', 'ValueType', 'any');
data(1)    =   1;
data(1.1)  =   2;
data(1.05) =   3;
data(3)    =   4;
data(3)    =   [data(3); 5];
data(1.5)  =   6;
data(10)   =   7;
data(10.1) =   8;
data(10.2) =   9;
data(10.3) =  10;
data(20)   =  11;
data(30)   =  12;

% min separation = 0.1
desired_output01 = containers.Map([1.05, 1.5, 3, 10, 10.1, 10.2, 10.3, 20, 30], ...
	{[1; 3; 2], [6], [4; 5], [7], [8], [9], [10], [11], [12]});
% min separation = 0.11
desired_output011 = containers.Map([1.05, 1.5, 3, 10.149999999999999, 20, 30], ...
	{[1; 3; 2], [6], [4; 5], [7; 8; 9; 10], [11], [12]});

dist = @(x, y) sqrt(sum(abs(x - y).^2));
clusters01 = mean_shift_dict_cluster(data, dist, 0.1);

clusters011 = mean_shift_dict_cluster(data, dist, 0.11);

if isequal(clusters01.keys, desired_output01.keys) ...
	&& isequal(clusters01.values, desired_output01.values) ...
	&& isequal(clusters011.keys, desired_output011.keys) ...
	&& isequal(clusters011.values, desired_output011.values)
	disp('Test passed');
else
	disp('Test failed');
end
