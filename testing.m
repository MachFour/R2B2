% this is just the guts of r2b2 plus some extras for easy use
audio_filename  = 'mollasses.wav';
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
display(a.test_aggregator.evidence_list);

%% PLOTTING TEMPO VS FRAME NO
tempos = [];
for i=1:length(tp_estimator{1}.tempo_phase_estimates)
    a = aggregator_demo;
    a.initialise(tp_estimator, i);
    
    a.test_aggregator.cluster_by_tempo();
    a.test_aggregator.cluster_by_phase();
    
    a.test_aggregator.collect_evidence();
    index = find(a.test_aggregator.evidence_list==max(a.test_aggregator.evidence_list));
    tempo = a.test_aggregator.tempo_centre_of_mass/a.test_aggregator.tempo_harmonics(index);
    tempo = 60/(tempo/a.sample_rate);
    tempos = [tempos, tempo];
end
scatter(1:length(tempos), tempos, 'b');