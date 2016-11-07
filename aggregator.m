%                           aggregator.m
% author:       daniel parker
% description:  this class defines how the tempo and phase estimates from 
%               the previous stage are 'aggregated' together to produce a 
%               single contextually influenced output of the true phase and 
%               tempo. clustering methods are used to find the most likely
%               tempo hypothesis, and then phase data is analysed to
%               measure evidence for different tempo hypotheses. once a
%               verdict is reached, a (tempo, phase) point is outputted and
%               then fed back into the aggregator to influence the next
%               aggregation window.

classdef aggregator < handle

properties (Constant)
    T_I = 1;    % tempo index
    P_I = 3;    % phase index
    C_I = 4;    % confidence index
    
    eps = 10;
end % properties (Constant

properties
    % this cell array records the tempo and phase estimates accross all
    % the features and is indexed by frame
    estimate_windows;
    w_i; % the current window index
    
    % allowed harmonics for 'tempo' and 'phase seperation'
    harmonics = [1,2,3,4]; 
    
    % this stores the hypothesis cluster 'tp_matrix' as generated by the
    % hypothesis cluster class, for each window. indexed by window
    hypothesis_data; 
   
    % the output data points from the aggregator as rows:
    % [tempo, phase, probability] indexed by window. used in the feedback 
    % mechanism
    tp_outputs;
    
end % properties

properties (Dependent)    
    % as measured from the size of the incoming data set
    num_features;

    % the most recent tempo_phase_estimate
    curr_tp_estimate;
    
end % properties (Dependent)

methods
    % ==== constructor ====
	function initialise(this)
        % no data in anything to start with
        this.w_i                = 0;
        this.estimate_windows   = cell(1);
        this.hypothesis_data    = cell(1);
        this.tp_outputs         = [];
    end

    % ==== getters ====
    function n = get.num_features(this)
        n = length(this.estimate_windows(this.w_i));
    end
    
    function point = get.curr_tp_estimate(this)
        point = this.tp_outputs(this.w_i,:);
    end
    
    % ==== general methods ====
    % this is the function that is called each time a new set of estimates
    % comes in. upon calling, the aggregator class will aggregate the data
    % in window together and then compute output estimate 
    function add_window(this, window, window_index)
        this.w_i = this.w_i + 1;
        % set the dimension to the number of features in the window
        this.estimate_windows{this.w_i} = cell(length(window));
        
        % transform the data into what is used by this class
        for i=1:this.num_features
            this.estimate_windows{this.w_i}{i} = ...
                window{i}.tempo_phase_estimates{window_index};
        end
        
        % instantiate the hypothesis cluster
        this.hypothesis_data{this.w_i} = hypothesis_cluster;
        this.hypothesis_data{this.w_i}.initialise(this.estimate_windows{this.w_i}, ...
            this.harmonics, this.num_features);
        
        [tempo_out, harmonic_index] = this.compute_winning_tempo();
        phase_out = this.compute_winning_phase(tempo_out, harmonic_index);
        conf_out  = this.compute_winning_confidence();
        
        this.tp_outputs = [this.tp_outputs; tempo_out, phase_out, conf_out];
    end
    
    % computes the winning tempo using collect evidence
    function [t, i_h] = compute_winning_tempo(this)
        e = this.collect_evidence();
        
        % find the harmonic index coresponding to the max
        i_h = find(e == max(e)); 
        t = this.hypothesis_data{this.w_i}.tempo_c_o_m/this.harmonics(i_h);
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
        for i=1:length(evidence_list)
            evidence_list(i) = ...
                sum(this.harmonics.*...
                this.hypothesis_data{this.w_i}.P_p_given_t_at_h(i,:))*...
                this.hypothesis_data{this.w_i}.t_clusters{i}.tot_conf;
        end
    end
    
    % this function finds the most confident 'epsilon' grouping of phases at
    % the winning tempo harmonic, and then computes a weighted mean of them
    % as the output phase. i'm not sure that this will work based on brief
    % testing but its a start!
    function phi = compute_winning_phase(this, index)
        % find the cluster containing all the phase points at the winning
        % tempo harmonic
        cluster = this.hypothesis_data{this.w_i}.t_clusters{index};
        
        test_cluster = tp_cluster;
        test_cluster.initialise(this.num_features)
        top_score = 0;
        
        % find the most confident epsilon ball
        for i=cluster.non_empty_features
            for j=1:cluster.n_pts(i)
                seed = cluster.get_point(i, j);
                most_confident_cluster = this.points_in_ball(seed, index);
                if test_cluster.tot_conf > top_score
                    top_score = test_cluster.tot_conf;
                    most_confident_cluster = test_cluster;
                end
            end
        end
        
        % take the weighted mean of that ball
        weighted_phases = zeros(sum(most_confident_cluster.n_pts), 1);
        k = 1;
        for i=most_confident_cluster.non_empty_features
            for j=1:most_confident_cluster.n_pts(i)
                point = most_confident_cluster.get_point(i, j);
                weighted_phases(k) = point(this.P_I)*point(this.C_I);
                k = k + 1;
            end
        end
        
        phi = sum(weighted_phases)/most_confident_cluster.tot_conf;
        
    end
    
    function cluster = points_in_ball(seed, index)
        % we check all the points in the tempo cluster and if they're
        % within eps of the seed, we add them to to the cluster
        tempo_cluster = this.hypothesis_data{this.w_i}.t_clusters{index};
        
        cluster = tp_cluster;
        cluster.initialise(this.num_features)
        
        % find the most confident epsilon ball
        for i=cluster.non_empty_features
            for j=1:cluster.n_pts(i)
                point = tempo_cluster.get_point(i, j);
                if abs(point(this.P_I) - seed(this.P_I)) < this.eps/2
                    cluster.add_point(point);
                end
            end
        end
    end
    
    % not sure what to do here just yet...
    function c = compute_winning_confidence(this)
        c = 1;
    end
    
end % methods
    
end % classdef