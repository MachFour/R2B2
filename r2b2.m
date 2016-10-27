% r2b2.m
% Author: Max Fisher

function r2b2(wav_name, output_directory)
%%%%%%%%%%%%%%%
%% (independently) tweakable parameters
%%%%%%%%%%%%%%%

% fix up to trim off extension
WAV_NAME = wav_name;
WAV_EXT = '.wav';

% assume mono audio wav file in current dir.

[y, Fs] = audioread(strcat('./', WAV_NAME, WAV_EXT));

%%%%%%%%%%%%
%% Calculate onset detection function
%%%%%%%%%%%%%

[df, FEATURE_RATE] = odf_klapuri(y, fs);




feature_frame = @(g, n) g((n-1)*FEATURE_HOP_SIZE + 1: min(((n-1)*FEATURE_HOP_SIZE + FEATURE_WIN_LENGTH), end), :);

feature_time_axis = (1:FEATURE_WIN_LENGTH)/FEATURE_RATE;

sample_frames = round(linspace(1, NUM_FEATURE_FRAMES, 4));
sample_frames = sample_frames(2:end-1);

%sample_frames = [sample_frames ; sample_frames + 1];

% use to disable plotting;
%sample_frames = [];


%%%%%%%%%%%%%%%%%%%%
%% Processing of feature vector for periodicity
% In reality everything would be done in the previous loop,
% also tempo estimates would have to be made for past feature measurements
%%%%%%%%%%%%%%%%%%%%%%


tempo_estimates = cell(NUM_FEATURE_FRAMES, NUM_FEATURE_CHANNELS);
tempo_confidences = cell(NUM_FEATURE_FRAMES, NUM_FEATURE_CHANNELS);
beat_alignment_estimates = cell(NUM_FEATURE_FRAMES, NUM_FEATURE_CHANNELS);
beat_alignment_confidences = cell(NUM_FEATURE_FRAMES, NUM_FEATURE_CHANNELS);

NUM_TEMPO_ESTIMATES = 4;

tic;

time = toc;
fprintf('spent %f s calculating beat_alignments\n', time);

%%%%%%%%%%%%%%%%%%%%
%% Write beat data to files
%%%%%%%%%%%%%%%%%%%%%
tic;

% export data in the following format
% one file per feature channel
for c = 1:NUM_FEATURE_CHANNELS;
    filename_for_channel = strcat('./', WAV_NAME, BEAT_DATA_OUTPUT_SUFFIX, ...
        sprintf('-channel-%d.txt', c));
    outfile = fopen(filename_for_channel, 'w+');

    end;
end;

time = toc;
fprintf('spent %f s writing beat data\n', time);

% vim: tabstop=4 expandtab
