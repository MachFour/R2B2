% r2b2.m
% This script shows how to use the different parts of the R2B2 system together to produce
% useful beat-tracking outputs

% Author: Max Fisher

function r2b2(audio_filename, audio_directory, data_output_directory)

	plot_data = 0;
	output_audio = 1;


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
	tp_estimator = cell(feature1.num_feature_channels, 1);

	for n = 1:feature1.num_feature_channels
		tp_estimator{n} = tpe_autocorrelation;
		tp_estimator{n}.initialise(feature1.feature_matrix(:, n), ...
			feature1.feature_sample_rate, sprintf('acf-feature1-ch%d', n));
	end

	% plot all the frames!
	num_sample_frames = tp_estimator{1}.num_feature_frames/10;
	sample_frames = round(linspace(1, ...
		tp_estimator{1}.num_feature_frames, num_sample_frames+2));
	%plot in reverse order so they'll appear on top of each other in
	%correct order
	sample_frames = [sample_frames(end-1:-1:1)];

	% don't plot anything
	sample_frames = [];
	if ~isempty(sample_frames)
		disp('Sample frames:');
		disp(sample_frames);
	end

	for n = 1:feature1.num_feature_channels
		tp_estimator{n}.compute_tempo_phase_estimates;
		if output_data
			tp_estimator{n}.output_tempo_phase_data(data_output_directory);
		end
	end

	% Viterbi algorithm

	num_features = feature1.num_feature_channels;
	num_feature_frames = tp_estimator{1}.num_feature_frames;
	params = processing_params;

	viterbi = bp_viterbi;
	viterbi.initialise(params, 'Viterwheats', num_features);

	% iterate over each frame, giving the viterbi algorithm its observations
	for k = 1:num_feature_frames
		kth_frame_tp_estimates = cell(1, num_features);
		kth_frame_features = zeros(params.feature_win_length, num_features);
		for n = 1:num_features
			% just take all the estimates, no clustering.
			feature_n_tp_estimates = tp_estimator{n}.tempo_phase_estimates{k};
			% put in a 2 column matrix of tempo and phase
			if ~isempty(feature_n_tp_estimates)
				kth_frame_tp_estimates{1, n} = ...
					[feature_n_tp_estimates(:, 1), feature_n_tp_estimates(:, 3)];
			end
			kth_frame_features(:, n) = tp_estimator{n}.get_feature_frame(k);
		end

		viterbi.step_frame(kth_frame_features, kth_frame_tp_estimates);
	end

	beat_times = viterbi.compute_beat_times;



	if output_data == 1
		% write beat times out to an annotation file
		filename = strrep(audio_filename, '.wav', '.txt');
		outfile = fopen(strcat(data_output_directory, '/', filename), 'w+');
		for timestamp = beat_times
			fprintf(outfile, '%.3f\n', timestamp);
		end
	else
		disp('Beat times');
		disp(beat_times);
	end

	if plot_data == 1
		% do scatter plot of all estimates in sample frames
		for k = sample_frames
			figure; hold on;
			for n = 1:feature1.num_feature_channels
				estimates_n = tp_estimator{n}.tempo_phase_estimates{k};
				%plot in seconds
				if isempty(estimates_n)
					warning('No estimates for feature %d in frame %d', n, k);
				else
					feature_fs = tp_estimator{n}.feature_sample_rate;
					tempos = estimates_n(:, 1);
					phases = estimates_n(:, 3);
					confidences = estimates_n(:, 4);
					stem3(tempos/feature_fs, phases/feature_fs, confidences, 'filled');
				end

			end
			title(sprintf('Scatterplot of tempo/phase estimates at %.2f s', ...
				tp_estimator{1}.get_frame_end_time(k)));
			xlabel('Tempo period (seconds)');
			ylabel('Offset from frame end (seconds)');
            zlabel('Confidence');
		end
	end

	if output_audio == 1
		% output audio annotations
		blip_track = mkblips(beat_times, audio_sample_rate, length(audio_data));

		%mix audio channels together if stereo
		if size(audio_data, 2) == 2
			audio_data = ((audio_data(:,1) + audio_data(:,2)) / 2);
		end

		% mix blips with audio, make blips slightly softer
		annotated_audio = (3*audio_data + 2*blip_track)/5;
		annotated_audio_filename = strcat('annotated_', audio_filename);
		annotated_audio_path = strcat(audio_directory, '/', annotated_audio_filename);
		audiowrite(annotated_audio_path, annotated_audio, audio_sample_rate);
	end

end
