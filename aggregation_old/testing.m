% this is just the guts of r2b2 plus some extras for easy use
audio_filename  = 'siarus.wav';
audio_directory = 'aggregation_old/test_songs'; 

audio_file_path = strcat(audio_directory, '/', audio_filename);
[audio_data, audio_sample_rate] = audioread(audio_file_path);

% give audio to feature class
feature1 = odf_klapuri;
feature1.initialise(audio_data, audio_sample_rate, ...
    'klapuris_feature');

feature1.compute_feature;

% make a tempo_phase estimator for each channel of feature1, since
% each channel is analysed independently
tp_estimator = cell(feature1.num_feature_channels, 1);

for n = 1:feature1.num_feature_channels
    tp_estimator{n} = tpe_autocorrelation;
    tp_estimator{n}.initialise(feature1.feature_matrix(:, n), ...
        feature1.feature_sample_rate, sprintf('acf-feature1-ch%d', n));
end

% pick some sample frames to plot
sample_frames = round(linspace(1, tp_estimator{1}.num_feature_frames, 4));
sample_frames = [sample_frames(2:end-1)];

for n = 1:feature1.num_feature_channels
    tp_estimator{n}.compute_tempo_phase_estimates;
end

%% AGGREGATION METHODS
%% CLUSTERING
clf
a = aggregator_demo;
a.initialise(tp_estimator, 10);
a.plot_tempo_phase_data();

% test tempo clustering
a.test_aggregator.cluster_by_tempo();
figure
a.plot_tempo_clustering_demo();

% test phase clustering
a.test_aggregator.cluster_by_phase();
figure
a.plot_phase_cluster_sets();

%% CHOOSING A WINNER
a.test_aggregator.collect_evidence();
%display(a.test_aggregator.evidence_list);

%% PLOTTING TEMPO VS FRAME NO AND THE WINNING PHASE
tempos = [];
phases = [];
for i=1:length(tp_estimator{1}.tempo_phase_estimates)
    a = aggregator_demo;
    a.initialise(tp_estimator, i);
    
    a.test_aggregator.cluster_by_tempo();
    a.test_aggregator.cluster_by_phase();
    
    a.test_aggregator.collect_evidence();
    index = find(a.test_aggregator.evidence_list==max(a.test_aggregator.evidence_list));
    if length(index) > 1
        tempo = 0;
    else    
        tempo = a.test_aggregator.tempo_centre_of_mass/a.test_aggregator.tempo_harmonics(index);
        tempo = 60/(tempo/a.sample_rate);
    end
    tempos = [tempos, tempo];
    
    % now we get the phase!
    top_score = 0;
    for j=1:a.test_aggregator.num_features
        cluster_sum = 0;
        if length(index) == 1
            cluster_sum = a.test_aggregator.get_sum(a.test_aggregator.phase_cluster_matrix{index,j}, 4);
        end
        if cluster_sum > top_score
            top_score = cluster_sum;
            phase_index = j;
        end
    end
    
    % now we take a weighted mean of the phases in that sub-sub-cluster
    %weighted_phases = [];
    %for k = 1:a.test_aggregator.num_features
    %    if length(index) == 1
    %        sz = size(a.test_aggregator.phase_cluster_matrix{index,phase_index}{k}.tempo_phase_estimates{i});
    %        for j = 1:sz(1)
    %            feature_data = a.test_aggregator.phase_cluster_matrix{index,phase_index}{k}.tempo_phase_estimates{i};
    %            if ~isempty(feature_data)
    %                weighted_phases = [weighted_phases, feature_data(j, 3)*feature_data(j, 4)];
    %            end
    %        end
    %    else
    %        weighted_phases = [weighted_phases, 0];
    %    end
    %end
    %if top_score == 0
    %    top_score = 1;
    %    weighted_phases = 0*weighted_phases;
    %end
    
    % in this approach we just take the most confident phase of the group
    top_score = 0;
    winning_phase = 0;
    for k=1:a.test_aggregator.num_features
        if length(index) == 1
            sz = size(a.test_aggregator.phase_cluster_matrix{index,phase_index}{k}.tempo_phase_estimates{i});
            for j = 1:sz(1)
                feature_data = a.test_aggregator.phase_cluster_matrix{index,phase_index}{k}.tempo_phase_estimates{i};
                if ~isempty(feature_data)
                    if feature_data(j, 4) > top_score
                        top_score = feature_data(j, 4);
                        winning_phase = feature_data(j, 3);
                    end
                end
            end
        end
    end
    
    %winning_phase = sum(weighted_phases/top_score); % we take the weigthed mean
    phases = [phases, winning_phase/(a.sample_rate)];
end
scatter(1:length(tempos), tempos, 'b');

%% WRITE TO ANNOTATION FILE...
filename = 'siarus.txt';
outfile = fopen(filename, 'w+');
for i=1:length(tempos) % or in otherwords, the number of 'windows'
    curr_tempo_period = 60/tempos(i); % tempos is in bpm
    num_beats = floor(5.94*0.125/curr_tempo_period);
    for j=1:num_beats
        curr_beat_time = 5.94 + (i - 1)*5.94*0.125 + j*curr_tempo_period + phases(i);
        fprintf(outfile, '%.3f\n', curr_beat_time);
    end
end