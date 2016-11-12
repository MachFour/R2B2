% Aggregator.m
% author:	   Daniel Parker
% Description:  This class defines how the tempo and phase estimates from
%			   the previous stage are 'aggregated' together to produce a
%			   single contextually influenced output of the true phase and
%			   tempo. clustering methods are used to find the most likely
%			   tempo hypothesis, and then phase data is analysed to
%			   measure evidence for different tempo hypotheses. once a
%			   verdict is reached, a (tempo, phase) point is outputted and
%			   then fed back into the aggregator to influence the next
%			   aggregation window.

classdef aggregator < handle

properties (Constant)
	T_I = 1;	% tempo index
	P_I = 3;	% phase index
	C_I = 4;	% confidence index

	eps = 20;

	num_features = 4;
end % properties (Constant

properties
	% this cell array records the tempo and phase estimates accross all
	% the features and is indexed by frame
	estimate_windows;
	w_i; % the current window index

	% allowed harmonics for 'tempo' and 'phase seperation'
	% need lots of harmonics now since the tempo can get real slow...
	harmonics = [1,2,3,4,5,6,8,12];

	% this stores the hypothesis cluster 'tp_matrix' as generated by the
	% hypothesis cluster class, for each window. indexed by window
	hypothesis_data;

	% the output data points from the aggregator as rows:
	% [tempo, phase, probability] indexed by window. used in the feedback
	% mechanism
	tp_outputs;

	% how many frames there are in a window. used to compute the feedback
	% point
	window_length;

end % properties

properties (Dependent)
	% the most recent tempo_phase_estimate
	curr_tp_estimate;

end % properties (Dependent)

methods
	% ==== constructor ====
	function initialise(this, window_length)
		% no data in anything to start with
		this.w_i				= 0;
		this.estimate_windows   = cell(1);
		this.hypothesis_data	= cell(1);
		this.tp_outputs		 = [];
		this.window_length	  = window_length;
	end

	% ==== getters ====
	function point = get.curr_tp_estimate(this)
		point = this.tp_outputs(this.w_i,:);
	end

	% ==== general methods ====
	% this is the function that is called each time a new set of estimates
	% comes in. upon calling, the aggregator class will aggregate the data
	% in window together and then compute output estimate
	function add_window(this, window, window_index)
		this.w_i = this.w_i + 1;

		if this.w_i == 1
			n_f = this.num_features;
		else % we add an extra feature for our feedback point
			n_f = this.num_features + 1;
		end

		% set the dimension to the number of features in the window
		this.estimate_windows{this.w_i} = cell(n_f, 1);

		% transform the data into what is used by this class
		for i = 1:n_f
			if i <= this.num_features
				this.estimate_windows{this.w_i}{i} = ...
					window{i}.tempo_phase_estimates{window_index};
			else
				% here we compute the estimate for the feedback observation
				past_tempo = this.tp_outputs(this.w_i - 1,1);
				past_phase = this.tp_outputs(this.w_i - 1,2);
				past_confidence = this.tp_outputs(this.w_i - 1,3);
				curr_phase = past_phase + past_tempo;
				while curr_phase < this.window_length/8 && past_tempo ~= 0
					curr_phase = curr_phase + past_tempo;
				end
				% now we are past the end of the next window. shave off a
				% past_tempo and subtract the window length to get a
				% negative phase
				curr_phase = curr_phase - past_tempo - this.window_length/8;
				% ^ distance from the end of the next window
				this.estimate_windows{this.w_i}{i} = ...
					[past_tempo, 0, curr_phase, past_confidence];
			end
		end

		% instantiate the hypothesis cluster
		this.hypothesis_data{this.w_i} = hypothesis_cluster;
		this.hypothesis_data{this.w_i}.initialise(this.estimate_windows{this.w_i}, ...
			this.harmonics, n_f);

		[tempo_out, harmonic_index] = this.compute_winning_tempo();
		phase_out = this.compute_winning_phase(harmonic_index);
		conf_out  = this.compute_winning_confidence();

		this.tp_outputs = [this.tp_outputs; tempo_out, phase_out, conf_out];
	end

	% computes the winning tempo using collect evidence
	function [t, i_h] = compute_winning_tempo(this)
		e = this.collect_evidence();

		% find the harmonic index coresponding to the max
		i_h = find(e == max(e));
		if length(i_h) > 1 || isempty(i_h)
			t = 0;
			i_h = []; % to stop glitches at the end...
		else
			t = this.hypothesis_data{this.w_i}.tempo_c_o_m/this.harmonics(i_h);
		end
	end

	% this function collects evidence for each tempo hypothesis based on
	% the phase structure and probabilities calculated in the
	% hypothesis_cluster class. the winning tempo is taken to be the
	% harmonic with the highest 'evidence' measure
	%
	% this is kind of arbitrary at the moment and could do with some
	% thinking about...
	function evidence_list = collect_evidence(this)
		evidence_list = zeros(length(this.harmonics), 1);

		for i = 1:length(evidence_list)
			evidence_list(i) = ...
				sum(this.harmonics.*...
				this.hypothesis_data{this.w_i}.P_p_given_t_at_h(i,:))*...
				this.hypothesis_data{this.w_i}.t_clusters{i}.tot_conf;
		end

		% weigh tempos by the log-normal distribution
		harmonics_list = (1./this.harmonics)*...
			this.hypothesis_data{this.w_i}.tempo_c_o_m;
		% this is arbitrary at the moment. i just drew a distribution that
		% looked good
		w = lognpdf(harmonics_list/172.2, 0.15, 0.7);
		evidence_list = evidence_list.*w';
	end

	% this function finds the most confident 'epsilon' grouping of phases at
	% the winning tempo harmonic, and then computes a weighted mean of them
	% as the output phase. i'm not sure that this will work based on brief
	% testing but its a start!
	function phi = compute_winning_phase(this, index)
		% find the cluster containing all the phase points at the winning
		% tempo harmonic
		if isempty(index) % to stop glitches at the end...
			phi = 0;
		else
			cluster = this.hypothesis_data{this.w_i}.t_clusters{index};

			top_score = 0;

			% find the most confident epsilon ball
			for i = cluster.non_empty_features
				for j = 1:cluster.n_pts(i)
					seed = cluster.get_point(i, j);
					test_cluster = this.points_in_ball(seed, index);
					if test_cluster.tot_conf > top_score
						top_score = test_cluster.tot_conf;
						most_confident_cluster = test_cluster;
					end
				end
			end

			% take the weighted mean of that ball
			weighted_phases = zeros(sum(most_confident_cluster.n_pts), 1);
			k = 1;
			for i = most_confident_cluster.non_empty_features
				for j = 1:most_confident_cluster.n_pts(i)
					point = most_confident_cluster.get_point(i, j);
					weighted_phases(k) = point(this.P_I)*point(this.C_I);
					k = k + 1;
				end
			end

			phi = sum(weighted_phases)/most_confident_cluster.tot_conf;
		end
	end

	function cluster = points_in_ball(this, seed, index)
		% we check all the points in the tempo cluster and if they're
		% within eps of the seed, we add them to to the cluster
		tempo_cluster = this.hypothesis_data{this.w_i}.t_clusters{index};

		cluster = tp_cluster;
		cluster.initialise(tempo_cluster.n_f);

		% find the most confident epsilon ball
		for i = tempo_cluster.non_empty_features
			for j = 1:tempo_cluster.n_pts(i)
				point = tempo_cluster.get_point(i, j);
				if abs(point(this.P_I) - seed(this.P_I)) < this.eps/2
					cluster.add_point(i, point);
				end
			end
		end
	end

	function c = compute_winning_confidence(this)
		e = this.collect_evidence();
		K = 0.77; % the feedback constant
		c = K*(max(e)/sum(e))*this.hypothesis_data{this.w_i}.P_t_clustering;
	end

end % methods

end % classdef