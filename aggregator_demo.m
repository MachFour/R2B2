% Class that performs tempo and phase estimation from a single dimensional
% feature, by using autocorrelation to detect likely periodicities (i.e.
% possible tempos), and then for each of those, using a sequence of impulses
% spaced at that tempo period to locate where the periodic peaks of the
% feature actually are (i.e. the beat phase)

classdef aggregator_demo < aggregator

properties (Constant)
    % for plotting: colours are used to distinguish clusters only
    COLOURS     = {'r', 'g', 'b', 'c', 'm', 'y', 'k', 'w'};
    % for plotting: point styles are used to distinguish features only
    POINT_STYLE = {'+', 'o', '*', 'x', 's', 'd', '^', '>', '<'};
end

properties
    % for testing purposes
    test_frame; 
    
    % for plotting (comes from max...)
    sample_rate; % in Hz
    
    test_aggregator;
end

methods
	function initialise(this, tp_data, test_frame)
        this.test_aggregator = aggregator; 
		this.test_aggregator.initialise(tp_data, test_frame);
        
        this.sample_rate = tp_data{1}.feature_sample_rate;
        this.test_frame = test_frame;
    end

    function plot_tempo_phase_data(this)
        hold on
        for i = 1:this.test_aggregator.num_features
            feature_data = this.test_aggregator.tp_data{i}.tempo_phase_estimates{this.test_frame};
            if ~isempty(feature_data)
                tempos = feature_data(:,this.test_aggregator.TEMPO_INDEX)/this.sample_rate;
                phases = feature_data(:,this.test_aggregator.PHASE_INDEX)/this.sample_rate;
                scatter(tempos, phases, this.POINT_STYLE{i}, 'k');
            end
        end
        hold off
    end
    
    function plot_cluster_set(this, cluster_set, phases)
        hold on
        for i = 1:length(cluster_set)
            for j = 1:this.test_aggregator.num_features
                feature_data = cluster_set{i}{j}.tempo_phase_estimates{this.test_frame};
                if ~isempty(feature_data)
                    if phases 
                        sz = size(feature_data);
                        tempos = 0.02*i + feature_data(:,this.test_aggregator.TEMPO_INDEX)/this.sample_rate;
                    else
                        tempos = feature_data(:,this.test_aggregator.TEMPO_INDEX)/this.sample_rate;
                    end
                    phases = feature_data(:,this.test_aggregator.PHASE_INDEX)/this.sample_rate;
                    scatter(tempos, phases, this.POINT_STYLE{j}, this.COLOURS{i});
                end
            end
        end        
        hold off
    end
    
    function plot_cluster(this, cluster, colour)
        hold on
        for i = 1:this.test_aggregator.num_features
            feature_data = cluster{i}.tempo_phase_estimates{this.test_frame};
            if ~isempty(feature_data)
                tempos = feature_data(:,this.test_aggregator.TEMPO_INDEX)/this.sample_rate;
                phases = feature_data(:,this.test_aggregator.PHASE_INDEX)/this.sample_rate;
                scatter(tempos, phases, this.POINT_STYLE{i}, colour);
            end
        end
        hold off
    end
    
    function plot_tempo_clustering_demo(this)
        t = this.test_aggregator.tempo_cluster_set;
        for i=1:length(this.test_aggregator.tempo_harmonics)
            t{i} = this.test_aggregator.zero_out_data(this.test_aggregator.tempo_cluster_set{i}, ...
                this.test_aggregator.PHASE_INDEX);
        end
        
        hold on
        for i = 1:this.test_aggregator.num_features
            feature_data = this.test_aggregator.tp_data{i}.tempo_phase_estimates{this.test_frame};
            if ~isempty(feature_data)
                tempos = feature_data(:,this.test_aggregator.TEMPO_INDEX)/this.sample_rate;
                sz = size(feature_data);
                phases = zeros(sz(1), 1);
                scatter(tempos, phases, this.POINT_STYLE{i}, 'k');
            end
        end
        hold off
        
        this.plot_cluster_set(t, false);
        
        hold on
            x = this.test_aggregator.centre_of_mass(this.test_aggregator.tempo_cluster_set,...
                this.test_aggregator.TEMPO_INDEX)/this.sample_rate;
            eps = this.test_aggregator.eps_tempo;
            for h=this.test_aggregator.tempo_harmonics
                line([x/h, x/h],[-0.1, 0.1], 'Color', 'Green');
                line([x/h + eps/(2*h*this.sample_rate), ...
                    x/h + eps/(2*h*this.sample_rate)],[-0.05, 0.05], 'Color', 'Red');
                line([x/h - eps/(2*h*this.sample_rate), ...
                    x/h - eps/(2*h*this.sample_rate)],[-0.05, 0.05], 'Color', 'Red');
            end
        hold off
    end
    
    function plot_phase_cluster_sets(this)
        for i=1:length(this.test_aggregator.tempo_harmonics)
            hold on
            this.plot_cluster_set({this.test_aggregator.phase_cluster_matrix{i,:}}, true);
            hold off

            %display(i);
            %this.test_aggregator.print_points_in_cluster_set({this.test_aggregator.phase_cluster_matrix{i,:}});
            
            %x = this.test_aggregator.centre_of_mass({this.test_aggregator.phase_cluster_matrix{i,:}}, this.test_aggregator.PHASE_INDEX);
            %line([0.3,1.3], [x/this.sample_rate,x/this.sample_rate], 'Color', 'green');
        end
    end
end % methods

end % classdef
