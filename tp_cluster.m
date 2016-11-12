% tp_cluster.m
% Author:	   Daniel parker
% Description:  tp_cluster = 'tempo phase cluster' is a data structure that
%			   holds (tempo, phase, confidence) points for use in clustering
%			   methods and cluster storage. it is accompanied by some
%			   useful methods to keep usage simple. not sure yet if all
%			   the getters are needed but the confidence ones are
%			   certainly used...

classdef tp_cluster < handle

properties (Constant)
	T_I = 1;	% tempo index
	P_I = 3;	% phase index
	C_I = 4;	% confidence index

end % properties (Constant)

properties
	tp_ests;	% a cell array indexed by feature, containing a matrix
				% of points with rows of the form [tempo, phase, confidence]

	n_f;		% the number of features present in this cluster

	n_pts;	  % an array containing the number of points for each
				% feature. saves typing when looping

	non_empty_features; % the set of non-empty features. useful for indexing
end % properties

properties (Dependent)
	min_conf;   % minimum confidence
	max_conf;   % maximum confidence
	tot_conf;   % sum of all the confidences accross all features

	min_phase;
	max_phase;

	min_tempo;  % these aren't used anywhere yet so we might not need them
	max_tempo;

end % properties (Dependent)

methods
	% ==== constructor ====
	function initialise(this, num_features)
		this.n_f     = num_features;
		this.tp_ests = cell(this.n_f, 1);
		this.n_pts   = zeros(1, this.n_f);
	end

	% ==== getters ====
	function c = get.min_conf(this)
		confs = zeros(this.n_f, 1);
		for i = this.non_empty_features
			confs(i) = min(this.tp_ests{i}(:, this.C_I));
		end
		c = min(confs);
	end

	function c = get.max_conf(this)
		confs = zeros(this.n_f, 1);
		for i = this.non_empty_features
			confs(i) = max(this.tp_ests{i}(:, this.C_I));
		end
		c = max(confs);
	end

	function c = get.tot_conf(this)
		confs = zeros(this.n_f, 1);
		for i = this.non_empty_features
			confs(i) = sum(this.tp_ests{i}(:, this.C_I));
		end
		c = sum(confs);
	end

	function p = get.min_phase(this)
		phases = zeros(this.n_f, 1);
		for i = this.non_empty_features
			phases(i) = min(this.tp_ests{i}(:, this.P_I));
		end
		p = min(phases);
	end

	function p = get.max_phase(this)
		phases = zeros(this.n_f, 1);
		for i = this.non_empty_features
			phases(i) = max(this.tp_ests{i}(:, this.P_I));
		end
		p = max(phases);
	end

	function t = get.min_tempo(this)
		phases = zeros(this.n_f, 1);
		for i = this.non_empty_features
			phases(i) = min(this.tp_ests{i}(:, this.T_I));
		end
		t = min(phases);
	end

	function t = get.max_tempo(this)
		phases = zeros(this.n_f, 1);
		for i = this.non_empty_features
			phases(i) = max(this.tp_ests{i}(:, this.T_I));
		end
		t = max(phases);
	end

	% ==== general methods ====

	% adds a point to tp_ests assuming point is formatted as
	% [tempo, phase, confidence]
	function add_point(this, feature, point)
		this.n_pts(feature) = this.n_pts(feature) + 1;
		this.tp_ests{feature} = [this.tp_ests{feature}; point];
		if ~any(this.non_empty_features == feature)
			this.non_empty_features = [this.non_empty_features, feature];
		end
	end

	% for examinining individual points
	function z = get_point(this, feature, index)
		z = this.tp_ests{feature}(index, :);
	end

end % methods

end % classdef
