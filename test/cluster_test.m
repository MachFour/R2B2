data = [
1,		1;
1.1,	2;
1.05,	3;
3,		4;
3,		5;
1.5,	6;
10,		7;
10.1,	8;
10.2,	9;
10.3,	10;
20,		11;
30,		12
];

% min separation = 0.1
desired_clusters01 = [
1.05
1.5
3
10
10.1
10.2
10.3
20
30
];

desired_metadata01 = {
[1; 3; 2]
6
[4; 5]
7
8
9
10
11
12
};

% min separation = 0.11
desired_clusters011 = [
1.05
1.5
3
10.149999999999999
20
30
];

desired_metadata011 = {
[1; 3; 2]
6
[4; 5]
[7; 8; 9; 10]
11
12
};


[clusters01, metadata01] = mean_shift_cluster(data, 0.1);

[clusters011, metadata011] = mean_shift_cluster(data, 0.11);

if isequal(clusters01, desired_clusters01) ...
	&& isequal(metadata01, desired_metadata01) ...
	&& isequal(clusters011, desired_clusters011) ...
	&& isequal(metadata011, desired_metadata011)
	disp('Test passed');
else
	disp('Test failed');
end
