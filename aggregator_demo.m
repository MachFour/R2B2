% Class that performs tempo and phase estimation from a single dimensional
% feature, by using autocorrelation to detect likely periodicities (i.e.
% possible tempos), and then for each of those, using a sequence of impulses
% spaced at that tempo period to locate where the periodic peaks of the
% feature actually are (i.e. the beat phase)

classdef aggregator_demo < aggregator

properties (Constant)
	% for plotting: colours are used to distinguish clusters only
	COLOURS = {'r', 'g', 'b', 'c', 'm', 'y', 'k', 'w'};
	% for plotting: point styles are used to distinguish features only
	POINT_STYLE = {'+', 'o', '*', 'x', 's', 'd', '^', '>', '<'};

	% assumed variables: these should really be got from the feature_extractor
	% which is being used!
	audio_sample_rate = 44100;
	audio_hop_size = 512; % samples
	feature_upsampling = 2;
	feature_sample_rate = audio_sample_rate/audio_hop_size*feature_upsampling;
end

methods
	function plot_cluster_set(this, cluster_set, phases)
		hold on
		for i = 1:length(cluster_set)
			if ~isempty(cluster_set{i})
				for j=cluster_set{i}.non_empty_features
					feature_data = cluster_set{i}.tp_ests{j};
					if phases
						tempos = 0.02*i + feature_data(:,cluster_set{i}.T_I)/this.feature_sample_rate;
					else
						tempos = feature_data(:,cluster_set{i}.T_I)/this.feature_sample_rate;
					end
					phases = feature_data(:,cluster_set{i}.P_I)/this.feature_sample_rate;
					scatter(tempos, phases, this.POINT_STYLE{j}, this.COLOURS{i});
				end
			end
		end
		hold off
	end

	% not fixed...
	function plot_cluster(this, cluster, colour)
		hold on
		for i = 1:this.test_aggregator.num_features
			feature_data = cluster{i}.tempo_phase_estimates{this.test_frame};
			if ~isempty(feature_data)
				tempos = feature_data(:,this.test_aggregator.TEMPO_INDEX)/this.feature_sample_rate;
				phases = feature_data(:,this.test_aggregator.PHASE_INDEX)/this.feature_sample_rate;
				scatter(tempos, phases, this.POINT_STYLE{i}, colour);
			end
		end
		hold off
	end

	function plot_tempo_clustering_demo(this, hypothesis_cluster, window)
		figure

		% zero out the phases
		h = hypothesis_cluster;
		for i=h.non_empty_t_clusters
			for j=h.t_clusters{i}.non_empty_features
				for k=1:h.t_clusters{i}.n_pts(j)
					h.t_clusters{i}.tp_ests{j}(k,h.P_I) = 0;
				end
			end
		end

		% plot all the data points on a number line
		hold on
		for i=1:h.n_f
			feature_data = window{i};
			if ~isempty(feature_data)
				tempos = feature_data(:,h.T_I)/this.feature_sample_rate;
				sz = size(feature_data);
				phases = zeros(sz(1), 1);
				scatter(tempos, phases, this.POINT_STYLE{i}, 'k');
			end
		end
		hold off

		this.plot_cluster_set(h.t_clusters, false);

		hold on
		x = h.tempo_c_o_m/this.feature_sample_rate;
		eps = h.eps_tempo;
		for i=h.non_empty_t_clusters
			k = h.harmonics(i);
			line([x/k, x/k],[-0.1, 0.1], 'Color', 'Green');
			line([x/k + eps/(2*k*this.feature_sample_rate), ...
				x/k + eps/(2*k*this.feature_sample_rate)],[-0.05, 0.05], 'Color', 'Red');
			line([x/k - eps/(2*k*this.feature_sample_rate), ...
				x/k - eps/(2*k*this.feature_sample_rate)],[-0.05, 0.05], 'Color', 'Red');
		end
		hold off
	end

	function plot_phase_cluster_sets(this, hypothesis_cluster)
		figure
		hold on
		for i=hypothesis_cluster.non_empty_t_clusters
			this.plot_cluster_set(hypothesis_cluster.tp_matrix(:,i), true);
		end
		hold off
	end
end % methods

end % classdef
