% r2b2.m
% Author: Max Fisher

function r2b2(audio_filename, audio_directory, data_output_directory)
	if nargin < 2
		audio_directory = '.';
	end
	if nargin < 3
		output_data = 0;
	else
		output_data = 1;
	end
	
	audio_file_path = strcat(audio_directory, '/', audio_filename);
	[audio_data, audio_sample_rate] = audioread(audio_file_path);
	
	% give audio to feature class
	feature1 = odf_klapuri;
	feature1.initialise(audio_data, audio_sample_rate, ...
		'klapuris_feature');
	
	feature1.compute_feature;
	
	% make a tempo_phase estimator for each channel of feature1, since
	% each channel is analysed independently
	tp_estimator = cell(feature1.num_feature_channels);
	

	
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
		tp_estimator{n}.plot_sample_intermediate_data(sample_frames);
		
		if output_data
			tp_estimator{n}.output_tp_estimates(data_output_directory);
		end
	end
		
	

	
	
