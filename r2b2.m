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

	%
	% ALGORITHM PARAMETERS
	%

	% (the following three values are commonly used in speech processing)
	audio_win_time = 20/1000; %seconds,
	audio_win_overlap_proportion = 0.5;
	audio_win_type = @hann;

	min_bpm = 20;
	max_bpm = 240;

	max_tempo_peaks = 8;
	max_alignment_peaks = 4;

	feature_win_time = 6; %seconds
	feature_win_overlap_proportion = 0.75;
	% or maybe use a window that weights recent samples more than older (by a
	% few seconds) samples
	feature_win_type = 'rect';

	feature_upsample_factor = 2;

	params = processing_params(audio_sample_rate, audio_win_time, ...
		audio_win_overlap_proportion, audio_win_type, feature_win_time, ...
		feature_win_overlap_proportion, feature_win_type, feature_upsample_factor, ...
		min_bpm, max_bpm, max_tempo_peaks, max_alignment_peaks);

	disp('algorithm parameters');
	disp(params);

	%
	% FEATURE CALCULATION
	%

	feature1 = odf_klapuri;
	feature1.initialise(audio_data, audio_sample_rate, ...
		'klapuris_feature');
	feature1.compute_feature;

	num_features = feature1.num_feature_channels;

	%
	% VITERBI ALGORITHM
	%

	ta_estimator = tae_autocorrelation(params, feature1.feature_matrix, 'tae-acf');
	viterbi = bp_viterbi(params, 'Viterwheats', num_features);

	viterbi.initialise

	% iterate over each frame, giving the viterbi algorithm its observations
	num_feature_frames = ta_estimator.num_feature_frames;

	for k = 1:num_feature_frames
		% get all features' estimates at once
		kth_feature_frame_matrix = ta_estimator.get_feature_frame(k);
		kth_frame_ta_estimates = ta_estimator.pick_tempo_and_alignment_estimates(k);
		viterbi.step_frame(kth_feature_frame_matrix, kth_frame_ta_estimates);
	end

	beat_times = viterbi.compute_beat_times;

	%
	% DATA OUTPUT
	%

	if output_data
		% write beat times out to an annotation file
		filename = strrep(audio_filename, '.wav', '.txt');
		outfile = fopen(strcat(data_output_directory, '/', filename), 'w+');
		for timestamp = beat_times
			fprintf(outfile, '%.3f\n', timestamp);
		end

		ta_estimator.output_tempo_alignment_data(data_output_directory);

	else
		disp('Beat times');
		disp(beat_times);
		disp('Winning states');
		for i = 1:length(viterbi.winning_states)
			if isa(viterbi.winning_states{i}, 'model_state')
				fprintf('frame %d: bpm = %.1f, observed time = %.2f, prob = %4g\n', ...
					viterbi.winning_states{i}.frame_number, ...
					viterbi.winning_states{i}.tempo_bpm, ...
					viterbi.winning_states{i}.beat_location, ...
					viterbi.winning_probabilities(i));
			end
		end
	end

	%
	% DATA PLOTTING
	%

	if plot_data == 1

		num_sample_frames = num_feature_frames/10;
		sample_frames = round(linspace(1, num_feature_frames, num_sample_frames+2));
		%plot in reverse order so they'll appear on top of each other in
		%correct order
		sample_frames = sample_frames(end-1:-1:1);

		% don't plot anything
		if ~isempty(sample_frames)
			disp('Sample frames:');
			disp(sample_frames);
		end


		% do scatter plot of all estimates in sample frames
		for k = sample_frames
			figure; hold on;
			for n = 1:num_features
				estimates_n = ta_estimator.tempo_alignment_estimates{k, n};
				%plot in seconds
				if isempty(estimates_n)
					warning('No estimates for feature %d in frame %d', n, k);
				else
					feature_fs = params.feature_sample_rate;
					tempos = estimates_n(:, 1);
					alignments = estimates_n(:, 3);
					confidences = estimates_n(:, 4);
					stem3(tempos/feature_fs, alignments/feature_fs, confidences, 'filled');
				end

			end
			title(sprintf('Scatterplot of tempo/beat alignment estimates at %.2f s', ...
				params.frame_end_time(k)));
			xlabel('Tempo period (seconds)');
			ylabel('Offset from frame end (seconds)');
            zlabel('Confidence');
		end
	end

	%
	% ANNOTETED AUDIO OUTPUT
	%

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
