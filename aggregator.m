% this class defines how the tempo and phase estimates from the previous
% stage are 'aggregated' together to produce a single contextually 
% influenced output of the true phase and tempo

classdef aggregator < handle

properties (Constant)
    % matrix columns coresponding to different types of data
    TEMPO_INDEX = 1;
    PHASE_INDEX = 3;
    CONFIDENCE_INDEX = 4;
end % properties (Constant

properties
    % this cell array contains the tempo and phase estimates accross all
    % the features
    tp_data;
    
    % the frame index we're up to -> remove for final implementation
    frame;
    
    % BPM ranges to allow when detecting tempo
	min_bpm = 40;
	max_bpm = 240;
    
    % epsilon ball radii used when clustering tempo and phase
    eps_tempo = 20;
    eps_phase = 20;
    
    % allowed harmonics for tempo and phase
    tempo_harmonics = [1,2,3,4];
    
    tempo_cluster_set;
    phase_cluster_matrix;
    
    evidence_list = [0,0,0,0];
end % properties

properties (Dependent)    
    % BPM ranges translated to autocorrelation time lag (seconds)
	% max lag needs to also be smaller than the frame time divided by 4, for
	% the tempo strength function to work properly
	max_lag;
	min_lag;
    
    num_features;
    
    tempo_centre_of_mass;
end % properties (Dependent)

methods
	function initialise(this, tp_data, frame)
        this.tp_data = tp_data;
        this.frame = frame;
	end

	% === GETTERS (for dependent properties) === 
	function l = get.max_lag(this)
		l = 60/this.min_bpm;
    end
    function l = get.min_lag(this)
		l = 60/this.max_bpm;
    end
    function n = get.num_features(this)
        n = length(this.tp_data);
    end
    function c = get.tempo_centre_of_mass(this)
        c = this.centre_of_mass(this.tempo_cluster_set, this.TEMPO_INDEX);
    end
    
    % === HELPER FUNCTIONS ===
    % returns the highest value found at index loc in the cluster given to
    % it (e.g. used to find the highest tempo estimate in a cluster).
    function largest = get_max(this, cluster, loc)
        feature_max = zeros(this.num_features, 1);
        for i = 1:this.num_features
            feature_data = cluster{i}.tempo_phase_estimates{this.frame};
            if ~isempty(feature_data)
                feature_max(i) = max(feature_data(:,loc));
            end
        end
        largest = max(feature_max);
    end
    
    % returns the smallest value found at index loc in the cluster given to
    % it.
    function smallest = get_min(this, cluster, loc)
        feature_min = zeros(this.num_features, 1);
        for i = 1:this.num_features
            feature_data = cluster{i}.tempo_phase_estimates{this.frame};
            if ~isempty(feature_data)
                feature_min(i) = min(feature_data(:,loc));
            end
        end
        smallest = min(feature_min);
    end
    
    % sums all values at index loc accross all features in a given cluster
    function s = get_sum(this, cluster, loc)
        feature_sum = zeros(this.num_features, 1);
        for i = 1:this.num_features
            feature_data = cluster{i}.tempo_phase_estimates{this.frame};
            if ~isempty(feature_data)
                feature_sum(i) = sum(feature_data(:,loc));
            end
        end
        s = sum(feature_sum);
    end
    
    % test whether a cluster is empty. useful for finding the highest
    % harmonic present. this is needed because isempty({{} {} {}}) returns
    % false but we consider the result empty in this case
    function bool = is_empty(this, cluster)
        bool = true;
        for i = 1:this.num_features
            if ~isempty(cluster{i}.tempo_phase_estimates{this.frame}) 
                bool = false;
                break;
            end
        end
    end
    
    % for debugging
    function print_points_in_cluster_set(this, cluster_set)
        for i=1:(length(cluster_set))
            for j=1:this.num_features
                feature_data = cluster_set{i}{j}.tempo_phase_estimates{this.frame};
                if ~isempty(feature_data)
                    fprintf('cluster: %d \t feature: %d\n', i, j);
                end
                disp(feature_data);
            end
        end
    end
    
    % used to compute either the tempo centre of mass or phase centre of
    % mass. centre of mass is defined to be a weighted mean of all the
    % points in a cluster once folded down to the fundamental tempo.
    % weights are simply the confidences. the input cluster_set is a list
    % of clusters where each index of this list coresponds to the same
    % index of the harmonic. for example if cluster_set = {c1, c2, c3} and
    % harmonics = [2, 4, 5] then c1 is paired up with 2, c2 with 4, etc.
    function x = centre_of_mass(this, cluster_set, loc)
        harmonics = this.tempo_harmonics;
    
        weighted_sum = 0;
        sum_of_confs = 0;
        
        for i=1:length(harmonics)
            for j=1:this.num_features
                feature_data = cluster_set{i}{j}.tempo_phase_estimates{this.frame};
                if ~isempty(feature_data)
                    if loc == this.TEMPO_INDEX
                        weighted_sum = weighted_sum + ...
                            sum(harmonics(i)*feature_data(:,loc).*feature_data(:,this.CONFIDENCE_INDEX));
                    elseif loc == this.PHASE_INDEX
                        % needs to be thought about...
                        weighted_sum = weighted_sum;
                    end
                end
            end
            sum_of_confs = sum_of_confs + this.get_sum(cluster_set{i}, this.CONFIDENCE_INDEX);
        end
        x = weighted_sum/sum_of_confs;
    end
    
    % this function will zero out one of the columns of the estimates
    % matrix. this is useful for when we want to plot in just one dimension
    % only
    function cluster_out = zero_out_data(this, cluster_in, loc)
        cluster_out = cluster_in; % copy over the structure
        for i = 1:this.num_features
            feature_data = cluster_in{i}.tempo_phase_estimates{this.frame};
            sz = size(feature_data);
            feature_data(:,loc) = zeros(sz(1), 1); % overwrite with zeros
            cluster_out{i}.tempo_phase_estimates{this.frame} = feature_data;
        end
    end    
        
    % === AGGREGATION METHODS ===
    % this method chooses a winning cluster by trying to generate clusters
    % using each point as a seed, rating them by confidence and then
    % keeping the most confident cluster as the winner
    function cluster_by_tempo(this)
        top_score = 0;
        for i=1:this.num_features
            curr_feature_data = this.tp_data{i}.tempo_phase_estimates{this.frame};
            sz = size(curr_feature_data);
            for j=1:sz(1) % loop over each row in curr_feature_data
                curr_point = curr_feature_data(j,:);
                
                % generate a cluster centred around curr_point
                test_cluster_set = this.tempo_query(curr_point);
                
                % measure its performance by computing the sum of
                % confidences over the whole cluster set
                curr_score = 0;
                for k=1:length(test_cluster_set)
                    curr_score = curr_score + ...
                        this.get_sum(test_cluster_set{k}, this.CONFIDENCE_INDEX);
                end
                if curr_score > top_score
                    top_score = curr_score;
                    this.tempo_cluster_set = test_cluster_set;
                end
            end
        end
    end
    
    % produces a 'cluster_set' organised by harmonics for a given seed
    % point
    function cluster_set = tempo_query(this, cluster_seed)
        % initialise an empty cluster set
        cluster_set = cell(length(this.tempo_harmonics), 1);
        for i=1:length(cluster_set)
            cluster_set{i} = this.new_cluster;
        end
    
        for i=1:this.num_features % iterate over all features
            curr_feature_data = this.tp_data{i}.tempo_phase_estimates{this.frame};
            sz = size(curr_feature_data);
            for j=1:sz(1) % iterate over all estimates
                curr_point = curr_feature_data(j,:);
                if  curr_point(this.TEMPO_INDEX) <= cluster_seed(this.TEMPO_INDEX)
                    [fold, harmonic] = this.fold_tempo(curr_point(this.TEMPO_INDEX), ...
                        cluster_seed(this.TEMPO_INDEX));
                    if fold > 0 % that is, if a folding actually occured...
                        index = find(this.tempo_harmonics==harmonic);
                        % append the current point to the cluster
                        % associated with the harmonic used
                        cluster_set{index}{i}.tempo_phase_estimates{this.frame} =  ...
                            [cluster_set{index}{i}.tempo_phase_estimates{this.frame}; curr_point];
                    end
                end
            end
        end   
    end

    % this function is just used to instantiate a new cluster so that the
    % data structure looks the same as the tempo_phase_estimates structure
    function cluster = new_cluster(this)
        cluster = cell(this.num_features, 1);
        for i=1:this.num_features
            cluster{i} = cluster_from_feature;
            % this is unecessarily big. figure out how to preserve the data
            % structure efficiently later...
            cluster{i}.tempo_phase_estimates = cell(this.frame, 1);
        end
    end
    
    % this function 'folds' the point x around the harmonics of x_ref until
    % it finds an appropriate epsilon ball to sit inside
    function [fold, harmonic] = fold_tempo(this, x, x_ref)
        fold = -1; 
        % if fold is negative, then x could not be folded into the harmonic 
        % ball set centred at x_ref        
        for harmonic = this.tempo_harmonics
            if (x*harmonic < x_ref + this.eps_tempo/2) && ...
                    (x*harmonic > x_ref - this.eps_tempo/2)
                fold = x*harmonic;
                break;
            end
        end
    end
    
    % for phase clustering
    function cluster_by_phase(this)
        this.phase_cluster_matrix = cell(length(this.tempo_harmonics));
        for i=1:length(this.tempo_harmonics)
            for j=1:length(this.tempo_harmonics)
                this.phase_cluster_matrix{i,j} = this.new_cluster;
            end
        end
        
        % we don't need to check all the tempo harmonics...
        for i=1:length(this.tempo_harmonics)
            for j=1:length(this.tempo_harmonics)
                top_score = 0;
                for k=1:this.num_features
                    feature_data = this.tempo_cluster_set{this.tempo_harmonics(i)}{k}.tempo_phase_estimates{this.frame};
                    sz = size(feature_data);
                    for m=1:sz(1)
                        cluster_seed = feature_data(m,:);
                        cluster = this.phase_query(cluster_seed, this.tempo_harmonics(j), i);
                        
                        % we rank each cluster by their global confidence
                        % and keep only the largest
                        curr_score = this.get_sum(cluster, this.CONFIDENCE_INDEX);
                        if curr_score > top_score
                            top_score = curr_score;
                            this.phase_cluster_matrix{i, j} = cluster;
                        end
                    end
                end
            end
        end
    end
    
    function total_size = get_cluster_size(this, cluster)
        total_size = 0;
        for i=1:this.num_features
            sz = size(cluster{i}.tempo_phase_estimates{this.frame});
            total_size = total_size + sz(1);
        end
    end
    
    % for a given tempo harmonic and phase harmonic this function finds the
    % largest subcluster such that phase centered epsilon intervals are
    % equally space by factors of the current tempo harmonic
    function cluster = phase_query(this, cluster_seed, tempo_harmonic, harmonic_index)
        % initialise an empty cluster set
        cluster = this.new_cluster;
    
        phase_data = this.tempo_cluster_set{harmonic_index};
        
        for i=1:this.num_features % iterate over all features
            curr_feature_data = phase_data{i}.tempo_phase_estimates{this.frame};
            sz = size(curr_feature_data);
            for j=1:sz(1) % iterate over all estimates
                curr_point = curr_feature_data(j,:);
                if  curr_point(this.PHASE_INDEX) <= cluster_seed(this.PHASE_INDEX)
                    offset =  this.tempo_centre_of_mass/tempo_harmonic;
                    num_pts = ceil((this.get_max(phase_data, this.PHASE_INDEX)...
                        - this.get_min(phase_data, this.PHASE_INDEX))/offset);
                    folded = this.fold_phase(curr_point(this.PHASE_INDEX), ...
                        cluster_seed(this.PHASE_INDEX), num_pts, offset);
                    if folded % that is, if a folding actually occured...
                        cluster{i}.tempo_phase_estimates{this.frame} =  ...
                            [cluster{i}.tempo_phase_estimates{this.frame}; curr_point];
                    end
                end
            end
        end 
    end
    
     % this function 'folds' the point x around the harmonics of x_ref until
    % it finds an appropriate epsilon ball to sit inside
    function folded = fold_phase(this, x, x_ref, num_pts, offset)
        folded = false; 
        % if fold is negative, then x could not be folded into the harmonic 
        % ball set centred at x_ref and will subsequently be excluded from 
        % the cluster    
        for h = 0:(num_pts - 1) 
            % negative signs are needed as phase is always negative
            if (x + h*offset < x_ref + this.eps_phase/2) && ...
                    (x + h*offset > x_ref - this.eps_phase/2)
                folded = true;
                break;
            end
        end
    end
    
    % this is the evidence metric for tempo. it's quite complicated and
    % will likely be a point of change in this approach
    function collect_evidence(this)
        for i=1:length(this.evidence_list)
            weighted_cluster_confidence = 0;
            for j=1:length(this.tempo_harmonics)
                if j < i 
                    % for tempo harmonic i we want the confidence of all 
                    % higher tempo harmonics that have a harmonic phase 
                    % structure that matches harmonic i
                    % we weight it by the harmonic number since we're more
                    % likely to see more phase harmonics at the fundamental
                    % than we are at higher harmonics <- THIS MAY CHANGE
                    % we take the ratio of a tempo clusters total
                    % confidence and the confidence of the phase harmonic
                    % sub-cluster
                    weighted_cluster_confidence = weighted_cluster_confidence...
                        + this.tempo_harmonics(j)*this.get_sum(...
                        this.phase_cluster_matrix{j,i},...
                        this.CONFIDENCE_INDEX)/this.get_sum(...
                        this.tempo_cluster_set{j},...
                        this.CONFIDENCE_INDEX);
                end
            end

            % we need the 1 here otherwise the fundamental gets no score
            % ever...
            this.evidence_list(i) = (1 + weighted_cluster_confidence)*this.get_sum(this.tempo_cluster_set{i}, this.CONFIDENCE_INDEX);
        end
    end
    
end % methods

% these methods come from any super-classes 
methods (Abstract)
    
end % methods (Abstract)
    
end % classdef