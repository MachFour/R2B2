% hypothesis_cluster.m
% Author:	   Daniel Parker
% Description: This class takes in a window of tempo and phase estimates,
%			   computes the general tempo hypothesis by clustering and
%			   then clusters the most likely tempo cluster by phase and
%			   stores this information in a matrix of clusters, indexed by
%			   harmonics. to use this class, simply instantiate to cluster
%			   the current window of estimates.

classdef hypothesis_cluster < handle

properties (Constant)
	T_I = 1;	% tempo index
	P_I = 3;	% phase index
	C_I = 4;	% confidence index

	eps_tempo = 20; % the size of the epsilon ball used for tempo in frames
	eps_phase = 20; % ditto for phase
end % properties (Constant)

properties
	n_f;		% the number of features

	harmonics;  % the array of allowed harmonics

	t_clusters; % single row of tempo clusters indexed by harmonic index

	tp_matrix;  % the matrix of clusters indexed {i, j}. columns corespond to
				% tempo harmonics and rows corespond to phase harmonics
				% at a particular tempo harmonic. if harmonics = [1,2,3]
				% then tp_matrix will look like
				% (1 1   1) < allowed divisions of phase at a tempo harmonic
				% (2 3/2 0)
				% (3 0   0)
				%  ^ fundamental tempo (j)
				% this matrix will always be triangular in shape so be
				% careful how you index it

	tested_tempos = []; % used in clustering to improve efficiency
	tested_phases = []; % ditto for phases

	P_t_clustering;	 % the probability of the winning tempo cluster
						% i.e. the sum of confidences in that cluster

	P_p_given_t_at_h;   % this is the probability of a phase clustering
						% given a tempo harmonic. i.e. it is the sum of
						% confidences of the winning phase cluster, divided
						% by the total confidence of that harmonic. it is a
						% matrix with the structure:
						% (1 0   0) < allowed divisions of phase at a tempo harmonic
						% (2 1   0)
						% (3 3/2 1)
						%  ^ fundamental tempo (j)
						% this is for easy use
end % properties

properties (Dependent)
	n;			  % the number of harmonics ~ dimensions of the cluster
					% matrix

	tempo_c_o_m;	% the tempo centre of mass at the harmonic. to get the
					% centre of a particular harmonic, divide by that
					% harmonic

	non_empty_t_clusters;   % a list of all non empty harmonics in the
							% tempo cluster. useful for indexing

end % properties (Dependent)

methods
	% ==== constructor ====
	function initialise(this, data_window, harmonics, num_features)
		this.n_f		= num_features;
		this.harmonics  = harmonics;

		this.t_clusters = cell(1, this.n);  % we probably don't need to do this
		this.tp_matrix  = cell(this.n);

		for j = 1:this.n
			for i = 1:this.n
				if i + j <= this.n + 1  % keep matrix triangular
					this.tp_matrix{i,j} = tp_cluster;
					this.tp_matrix{i,j}.initialise(num_features);
				end
			end
			% initialise the tempo clusters
			this.t_clusters{j} = tp_cluster;
			this.t_clusters{j}.initialise(num_features);
		end

		this.P_p_given_t_at_h = zeros(this.n);

		this.cluster_by_tempo(data_window);
		this.cluster_by_phase();
	end

	% ==== getters ====
	function n = get.n(this)
		n = length(this.harmonics);
	end

	% computes the tempo (period) centre of mass. divide this by the
	% harmonics to get the tempo desired
	function t = get.tempo_c_o_m(this)
		weighted_sum = 0;
		sum_of_confs = 0;

		for i = 1:this.n
			h = this.harmonics(i);
			for j = this.t_clusters{i}.non_empty_features
				feature_data = this.t_clusters{i}.tp_ests{j};
				weighted_sum = weighted_sum + ...
					sum(h*feature_data(:,this.T_I).*feature_data(:,this.C_I));
			end
			sum_of_confs = sum_of_confs + this.t_clusters{i}.tot_conf;
		end
		t = weighted_sum/sum_of_confs;
	end

	function s = get.non_empty_t_clusters(this)
		s = [];
		for i = 1:this.n
			if ~isempty(this.t_clusters{i}.non_empty_features)
				s = [s, i];
			end
		end
	end

	% ==== general methods ====

	% also used to make indexing easier less bug prone. this function may
	% not actually end up getting used but the getter above is for sure...
	function s = non_empty_p_clusters(this, t_harmonic)
		s = [];
		for i = 1:this.n
			if i + t_harmonic < this.n + 1
				if ~isempty(this.tp_matrix{i, t_harmonic}.non_empty_features)
					s = [s; i];
				end
			end
		end
	end

	% this method chooses a winning harmonic cluster by trying to generate
	% clusters using each point as a seed, rating them by confidence and then
	% keeping the most confident cluster as the winner
	function cluster_by_tempo(this, data_window)
		top_score = 0; % for choosing a winning cluster

		for i = 1:this.n_f % iterate over each feature
			curr_feature_data = data_window{i};
			sz = size(curr_feature_data);
			for j = 1:sz(1) % loop over each row in the window table
				curr_point = curr_feature_data(j,:);

				% this stops us from re-testing a bunch of points that all
				% have the same tempo. makes the clustering much faster
				if ~any(this.tested_tempos == curr_point(this.T_I))
					% generate a cluster centred around curr_point
					clusters = this.tempo_query(curr_point, data_window);

					% measure its performance by computing the sum of
					% confidences over the whole cluster set
					curr_score = 0;
					for k=1:this.n % loop over the harmonics
						curr_score = curr_score + clusters{k}.tot_conf;
					end
					if curr_score > top_score
						top_score = curr_score;
						this.t_clusters = clusters;
					end
				end

				% update the list of tempos we've tested
				if isempty(this.tested_tempos)
					this.tested_tempos = curr_point(this.T_I);
				else
					this.tested_tempos = [this.tested_tempos, curr_point(this.T_I)];
				end
			end
		end
		this.P_t_clustering = top_score;
	end

	% produces a 'cluster_set' organised by harmonics for a given seed
	% point
	function cluster_set = tempo_query(this, cluster_seed, tempo_window)
		% initialise an empty cluster set indexed by harmonics
		cluster_set = cell(1, this.n);
		for i = 1:this.n
			cluster_set{i} = tp_cluster;
			cluster_set{i}.initialise(this.n_f);
		end

		for i = 1:this.n_f % iterate over all features
			curr_feature_data = tempo_window{i};
			sz = size(curr_feature_data);
			for j = 1:sz(1) % iterate over all estimates
				curr_point = curr_feature_data(j,:);

				% we may want to change this later to include harmonics on
				% the other side but for now this is 'precise' enough
				if  curr_point(this.T_I) <= cluster_seed(this.T_I)
					[d, i_h] = this.d_T(curr_point(this.T_I), cluster_seed(this.T_I));
					if d < this.eps_tempo/2 % is it in the epsilon ball?
						% append the current point to the cluster
						cluster_set{i_h}.add_point(i, curr_point);
					end
				end
			end
		end
	end

	% this function is the harmonic metric we use to measure the distance
	% between two tempos. index coresponds to the index of the harmonic
	% used to make the fold
	function [dist, index] = d_T(this, x, x_ref)
		[dist, index] = min(abs(x*this.harmonics - x_ref));
	end

	% for phase clustering. excuse the nested loops :P
	% this function uses t_clusters once they have been computed and then
	% splits each tempo harmonic by phase into clusters that may give
	% 'evidence' to other harmonics. for example if the fundamental is a
	% period of 1s with phase groupings seperated by 0.5s, then this would
	% give evidence for the tempo period of 0.5s. this function goes an
	% finds all the maximal groups at each harmonic that show the temporal
	% seperation as other tempo harmonics, and stores all of these clusters
	% in tp_matrix. it also computes the probability of a phase clustering
	% given a tempo harmonic, which is stored in the matrix
	% P_p_given_t_at_h
	function cluster_by_phase(this)
		for i = this.non_empty_t_clusters % iterate over non-empty tempo harmonics
			% now we check phase harmonics at tempo harmonic i
			for j=this.non_empty_t_clusters
				if i + j <= this.n + 1  % keeping it triangular
					top_score = 0; % for updating the winner

					for k=this.t_clusters{i}.non_empty_features
						for m=1:this.t_clusters{i}.n_pts(k)
							cluster_seed = this.t_clusters{i}.get_point(k, m);

							% phase offset is the spacing allowed between
							% phase points in the same epsilon ball in the
							% tempo cluster at harmonic i
							%
							% this isn't quite right and needs fixing...
							phase_offset = this.tempo_c_o_m/this.harmonics(j);
							cluster = this.phase_query(cluster_seed, ...
								i, phase_offset);

							% we rank each cluster by their global confidence
							% and keep only the largest
							curr_score = cluster.tot_conf;
							if curr_score > top_score
								top_score = curr_score;
								% phase is rows, columns are tempo harmonics
								this.tp_matrix{j,i} = cluster;

								% update probability matrix
								this.P_p_given_t_at_h(i + j - 1,i) = ...
									top_score/this.t_clusters{i}.tot_conf;
							end
						end
					end
				end
			end
		end
	end

	% given a phase seed, tempo harmonic and phase offset, this function
	% finds all the points in t_clusters(tempo harmonic) that are an
	% integer multiple of the offset away from the phase seed
	function cluster = phase_query(this, cluster_seed, i_tempo_harmonic, ...
			phase_offset)
		% initialise an empty cluster set
		cluster = tp_cluster;
		cluster.initialise(this.n_f);

		% this is the tempo cluster that the cluster_seed came from
		phase_data = this.t_clusters{i_tempo_harmonic};

		for i = phase_data.non_empty_features % iterate over each feature
			for j = 1:phase_data.n_pts(i) % iterate over all estimates
				curr_point = phase_data.get_point(i, j);

				% here we only want to look at points that are further down
				% the tempo phase plot, than the cluster seed
				if curr_point(this.P_I) <= cluster_seed(this.P_I)
					% compute the actual number of levels that we need to
					% check
					num_pts = ceil((cluster_seed(this.P_I) - ...
						phase_data.min_phase)/phase_offset);
					d = this.d_phi(curr_point(this.P_I), ...
						cluster_seed(this.P_I), num_pts, phase_offset);
					if d < this.eps_phase/2 % is it in the epsilon ball?
						cluster.add_point(i, curr_point);
					end
				end
			end
		end
	end

	% this is the metric that we use to compute the distance between two
	% phase points modulo the offset
	function dist = d_phi(this, y, y_ref, num_pts, offset)
		h = 0:(num_pts - 1);
		dist = min(abs(y + h*offset - y_ref));
	end

end % methods

end % classdef
