% r2b2.m
% Author: Max Fisher

%%%%%%%%%%%%%%%
%% (independently) tweakable parameters
%%%%%%%%%%%%%%%
WAV_NAME = 'walking-short';
WAV_EXT = '.wav';

% window parameters, in samples
% approximate time of audio to use for windowing.
% Window size in samples will be rounded to the nearest
% power of 2.

AUDIO_WIN_TIME_APPROX = 20/1000; %seconds,
AUDIO_WIN_OVERLAP_PERCENT = 50;
%these are common values in speech processing

AUDIO_WIN_TYPE = 'Hann';

% feature parameters

% minimum length of audio to use for tempo estimation at each time point (seconds)
% actual length of feature frames will be power of 2 closest to this.
FEATURE_TIME_APPROX = 5; %seconds

% Approximately how often to make estimates of new tempo by recalculating
% the autocorrelation. This is used to calculate the feature hop size
%TEMPO_UPDATE_FREQUENCY_APPROX = 2; % Hz

FEATURE_WIN_OVERLAP_PERCENT = 87.5;

% or maybe use a window that weights recent samples more than older (by a few seconds) samples
FEATURE_WIN_TYPE = 'Rect';

% how many bands to use for periodicity analysis
NUM_FEATURE_CHANNELS = 4;

% BPM ranges to allow when detecting tempo
MIN_BPM = 40;
MAX_BPM = 240;

% BPM ranges translated to autocorrelation lag
% max lag needs to also be smaller than the frame time divided by 4, for
% the tempo strength function to work properly
MAX_LAG = 60/MIN_BPM;
MIN_LAG = 60/MAX_BPM;


% output the feature values to a time-value CSV file with the following suffix
FEATURE_OUTPUT_SUFFIX = '-feature';

% output the tempo and beat alignment values to a text file with the following suffix
BEAT_DATA_OUTPUT_SUFFIX = '-beat-data';

%%%%%%%%%%%%%%%%%%%%%%%%
%% calculated parameters
% (these are completely derived from the above ones)
%%%%%%%%%%%%%%%%%%%%%%%%

% assume mono audio wav file in current dir.

[y, Fs] = audioread(strcat('./', WAV_NAME, WAV_EXT));

AUDIO_WIN_LENGTH = 2^round(log2(Fs*AUDIO_WIN_TIME_APPROX));
AUDIO_WIN_TIME = AUDIO_WIN_LENGTH/Fs;
% number of samples advanced by each window
AUDIO_HOP_SIZE = AUDIO_WIN_LENGTH*(1-AUDIO_WIN_OVERLAP_PERCENT/100);

NUM_AUDIO_FRAMES = ceil(length(y)/AUDIO_HOP_SIZE);

% sampling rate of the feature measurement
FEATURE_INTERP = 2;
FEATURE_RATE = FEATURE_INTERP*Fs/AUDIO_HOP_SIZE; % times 2 after upsampling

% how many past feature calculation samples to use when estimating current tempo
FEATURE_WIN_LENGTH = 2^round(log2(FEATURE_TIME_APPROX*FEATURE_RATE));
% actual length of audio used, as corresponding to the FEATURE_WIN_LENGTH
FEATURE_WIN_TIME = FEATURE_WIN_LENGTH/FEATURE_RATE;

%FEATURE_HOP_SIZE = 2^round(log2(FEATURE_RATE/TEMPO_UPDATE_FREQUENCY_APPROX));
FEATURE_HOP_SIZE = FEATURE_WIN_LENGTH*(1-FEATURE_WIN_OVERLAP_PERCENT/100);
TEMPO_UPDATE_FREQUENCY = FEATURE_RATE/FEATURE_HOP_SIZE; % Hz

NUM_FEATURE_FRAMES = ceil(NUM_AUDIO_FRAMES/FEATURE_HOP_SIZE);

fprintf('Audio sample rate: \t %3.2f kHz\n', Fs/1000);
fprintf('Audio window time: \t %f ms\n', AUDIO_WIN_TIME*1000);
fprintf('Audio window length: \t %d samples\n', AUDIO_WIN_LENGTH);
fprintf('Audio hop size: \t %d samples\n', AUDIO_HOP_SIZE);
fprintf('Number of audio frames: \t %d\n', NUM_AUDIO_FRAMES);
fprintf('\n');

fprintf('Feature sample rate: %f Hz', FEATURE_RATE);
if FEATURE_INTERP > 1
    fprintf(' (after interpolating by %d)\n', FEATURE_INTERP);
else
    fprintf('\n');
end

fprintf('Feature window time: \t %f s \n', FEATURE_WIN_TIME);
fprintf('Feature window length: \t %d samples\n', FEATURE_WIN_LENGTH);
fprintf('Feature hop size: %d samples\n', FEATURE_HOP_SIZE);
fprintf('Number of feature frames: \t %d\n', NUM_FEATURE_FRAMES);

fprintf('Tempo Update Frequency \t: %3.3f Hz\n', TEMPO_UPDATE_FREQUENCY);

if MAX_LAG > FEATURE_WIN_TIME/4 - 30/FEATURE_RATE
    disp('Warning:');
    disp('Reduced MAX_LAG to be 30 samples less than 1/4 of FEATURE_WIN_TIME');
    MAX_LAG = FEATURE_WIN_TIME/4 - 30/FEATURE_RATE;
    fprintf('New MAX_LAG: %f\n', MAX_LAG);
    MIN_TEMPO = 60/MAX_LAG;
    fprintf('New MIN_TEMPO: %f\n', MIN_TEMPO);
end;



%%%%%%%%%%%%%%
%% algorithm setup and variable initialisation
%%%%%%%%%%%%%%


%y[n] is the audio sample
%f[n] is the feature vector

%%%%%%%%%%%%%
%Mel filterbank
% used for MFCC's and also for Klapuri's features
%%%%%%%%%%%%
NUM_MEL_FILTERS = 36;
mel_filterbank = melfilter(Fs, AUDIO_WIN_LENGTH, NUM_MEL_FILTERS);
% returns a NUM_MEL_FILTERS*AUDIO_WIN_LENGTH/2 size matrix,
% perfect for left multiplying the DFT magnitude spectrum by

compressed_mel_band_energy = zeros(NUM_AUDIO_FRAMES, NUM_MEL_FILTERS);

% for now, store all the feature calculations in df[k, n], but only use the last
% FEATURE_WIN_LENGTH many to compute tempo at each time step

% differenced feature
df = zeros(NUM_AUDIO_FRAMES, NUM_FEATURE_CHANNELS);
% previous log power spectrq for differencing
prev_log_power_spectrum = zeros(AUDIO_WIN_LENGTH, 1);
% second difference
%prev2_log_power_spectrum = zeros(AUDIO_WIN_LENGTH, 1);

% choice of windowing function.
win = window(@hann, AUDIO_WIN_LENGTH);

% handle to get the correct set of samples for frame n; 1-indexed.
audio_frame = @(n) y(((n-1)*AUDIO_HOP_SIZE + 1): min(((n-1)*AUDIO_HOP_SIZE + AUDIO_WIN_LENGTH), end), 1);
% handle to get sample numbers for a given frame

% low pass filter for filtering the log magnitude spectrum
% [4] specifies a 10Hz cutoff frequency for the filter, but we're filtering the spectrum...
% what's the sampling rate?
% -> related to the frequency resolution
delta_f = Fs/AUDIO_WIN_LENGTH; %  ~= 44Hz -> this is the distance (in Hz frequency) between DFT samples
% -> thus the frequency 'sample rate' is 1/44 (Hz)^-1 (i.e. seconds!)
% -> 44Hz! ...  does this mean that frequency resolution is worse than an entire octave 
% at the lower end of hearing?? Yes -> we should really use the Constant-Q transform here
% -> Or filter with logarithmically spaced bandpass filters, like Klapuri


LP_cutoff_freq = 15; %Hz

LP_cutoff_freq_normalised = LP_cutoff_freq/(FEATURE_RATE/2);
LP_order = 6;
[LP_num, LP_denom] = butter(LP_order, LP_cutoff_freq_normalised);
% normalise by ratio of frequencies (again following [4])
LP_num = LP_num/LP_cutoff_freq_normalised;

feature_frame = @(g, n) g((n-1)*FEATURE_HOP_SIZE + 1: min(((n-1)*FEATURE_HOP_SIZE + FEATURE_WIN_LENGTH), end), :);

feature_time_axis = (1:FEATURE_WIN_LENGTH)/FEATURE_RATE;

sample_frames = round(linspace(1, NUM_FEATURE_FRAMES, 4));
sample_frames = sample_frames(2:end-1);

%sample_frames = [sample_frames ; sample_frames + 1];

% use to disable plotting;
%sample_frames = [];

%%%%%%%%%%%%%%%%%
%% Processing of raw audio frames
%%%%%%%%%%%%%%%%%
tic;

% we compute Klapuri's feature from [4]
for k = 1:NUM_AUDIO_FRAMES
    % last frames will not be complete; if so, truncate the window
    curr_frame = audio_frame(k).*win(1:length(audio_frame(k)));
    mean = sum(curr_frame)/length(curr_frame);
    variance = sum(curr_frame.^2)/(length(curr_frame) - 1);
    % normalise to variance 1
    curr_frame = (curr_frame - mean)/sqrt(variance);
    
    power_spectrum = abs(fft(curr_frame, AUDIO_WIN_LENGTH));
    % throw away (symmetrical) right half of spectrum
    power_spectrum = power_spectrum(1:AUDIO_WIN_LENGTH/2);

    mel_band_energy = mel_filterbank * power_spectrum;

    % calculate log magnitude energy spectrum with parameter mu = 100 (
    mu = 100;

    compressed_mel_band_energy(k, :) = log(1 + mu*mel_band_energy)/log(1 + mu);
end;
time = toc;
fprintf('spent %f s calculating band energies\n', time);

tic;
% here Klapuri performs a 2x interpolation of the compressed mel band energy over time before calculating the difference between frames
% have to find a way of doing this inside the loop
compressed_mel_band_energy = upsample(compressed_mel_band_energy, 2);
compressed_mel_band_energy = filter(LP_num, LP_denom, compressed_mel_band_energy, [], 1);
time = toc;

fprintf('spent %f s interpolating feature\n', time);

for k = 2:NUM_AUDIO_FRAMES
    % take the half-wave rectified difference between this and prev frame's
    % log power spectrum, to create a feature that peaks at the onset time

    % Bock: second difference so that we take the difference with a 50%
    % overlapping window, but have double the feature rate
    
    unrectified_difference = compressed_mel_band_energy(k, :) - compressed_mel_band_energy(k-1, :);
    HWR_difference = unrectified_difference;
    negative_entries = find(HWR_difference < 0);
    HWR_difference(negative_entries) = 0;

    % divide up filters evenly into NUM_FEATURE_CHANNELS channels
    filter_divisions = round(linspace(1,NUM_MEL_FILTERS, NUM_FEATURE_CHANNELS+1));
    
    for c = 1:NUM_FEATURE_CHANNELS
        mel_filters_in_channel = filter_divisions(c):filter_divisions(c+1);
        undifferenced_band_energy = sum(compressed_mel_band_energy(mel_filters_in_channel));
        differenced_band_energy = sum(HWR_difference(mel_filters_in_channel));
        % weighted sum the result with the undifferentiated log_power_spectrum
        % makes analysis better according to [4]
        lambda = 1;%0.8;
        df(k, c) = (1 - lambda)*undifferenced_band_energy + ...
            lambda*differenced_band_energy;
    end;
    
%     if find(k == sample_frames)
%         figure;
%         subplot(2, 1, 1);
%         plot(compressed_mel_band_energy(k, :));
%         title(sprintf('mel band energies, k = %d', k));
%         subplot(2, 1, 2);
%         plot(HWR_difference);
%         title(sprintf('HWR difference of mel band energies, k = %d', k));
%     end;

end;
time = toc;
fprintf('spent %f s calculating feature differential \n', time);

%%%%%%%%%%%%%%%%%%%%
%% Write feature data to files
%%%%%%%%%%%%%%%%%%%%%

% export f to time instants to import into Sonic Visualiser
% f[k] represents the difference between frame k and frame k-1.
% However the frames overlap - at what point should this measurement be
% attributed to?
% in terms of time, we say that this value corresponds to the end of frame k-1.
% Then we put f(k) at (k-1)/FEATURE_RATE seconds

% outfile = cell(FEATURE_CHANNELS, 1);
% for band = 1:FEATURE_CHANNELS
%     filename_for_band = strcat('./', WAV_NAME, FEATURE_OUTPUT_SUFFIX, sprintf('-%d.txt', band));
%     outfile{band} = fopen(filename_for_band, 'w+');
% end
% 
% for k = 1:length(f)
%     for band = 1:FEATURE_CHANNELS
%         fprintf(outfile{band}, '%f\t%f\n', (k-1)/FEATURE_RATE, f(k, band));
%     end
% end;

%%%%%%%%%%%%%%%%%%%%
%% Processing of feature vector for periodicity
% In reality everything would be done in the previous loop,
% also tempo estimates would have to be made for past feature measurements
%%%%%%%%%%%%%%%%%%%%%%

% beat hypothesis data
% the nth index of the array in each cell corresponds to a single
% (tempo, tempo confidence, beat location, beat location confidence) tuple
% which forms a complete hypothesis of a particular feature at a particular
% frame

tempo_estimates = cell(NUM_FEATURE_FRAMES, NUM_FEATURE_CHANNELS);
tempo_confidences = cell(NUM_FEATURE_FRAMES, NUM_FEATURE_CHANNELS);
beat_alignment_estimates = cell(NUM_FEATURE_FRAMES, NUM_FEATURE_CHANNELS);
beat_alignment_confidences = cell(NUM_FEATURE_FRAMES, NUM_FEATURE_CHANNELS);

NUM_TEMPO_ESTIMATES = 4;

tic;
for k = 1:NUM_FEATURE_FRAMES

    % go through and calculate autocorrelations for each slice of 
    % FEATURE_LENGTH detection function features
    % no windowing so far
    curr_df_frame = feature_frame(df, k); % no windowing so far
    frame_length = length(curr_df_frame);
    if frame_length < FEATURE_WIN_LENGTH/2
        % probs end of song, ignore
        continue;
    end;

    % actual lag time, trimmed to length of current ODF frame
    curr_time_axis = feature_time_axis(1:frame_length);

    % lag converted to tempo in BPM, flip to make increasing
    %curr_tempo_axis = fliplr(60./curr_time_axis);

    % separate autocorrelations for each channel

    acf = zeros(frame_length, NUM_FEATURE_CHANNELS);
    for c = 1:NUM_FEATURE_CHANNELS
        ac = autocorrelation(curr_df_frame(:, c));
        acf(:, c) = ac;
    end;

    % range of lags to consider (in samples)
    candidate_tempos = round(FEATURE_RATE*MIN_LAG) : ...
        min(round(FEATURE_RATE*MAX_LAG),frame_length);
    
    % now pick peaks!
    % note that there may be peaks at smaller lag than MAX_BPM that reflect the tatum
    % (semiquaver/quaver) pulse
    % find peaks until confidence is within a certain ratio of the highest  
    for c = 1:NUM_FEATURE_CHANNELS
        % MinPeakProminence: necessary vertical descent on both sides in order to count as
        % a peak
        % MinPeakDistance: horizonal separation between peaks.
      [tempo_confidences{k, c}, tempo_estimates{k, c}]  = ...
            findpeaks(acf(candidate_tempos, c), candidate_tempos, ...
                'MinPeakProminence', 0.05, ...
                'MinPeakDistance', MIN_LAG/2, ...
                'MinPeakHeight', 0, ...
                'NPeaks', NUM_TEMPO_ESTIMATES, ...
                'SortStr', 'descend');
        % remove peak at 0, and peaks with confidence below 0.
    end;
    
    % alternative: 'SHIFT-INVARIANT-COMB-FILTERBANK'
    % multiply with periodic comb for each possible tempo in range
    % range of comb separations to try (in samples)

    
    tempo_strength = zeros(length(candidate_tempos), NUM_FEATURE_CHANNELS);
    for c = 1:NUM_FEATURE_CHANNELS
        for m = 1:length(candidate_tempos)
            acf_comb = autocorrelation_comb(length(curr_time_axis), candidate_tempos(m));
          tempo_strength(m, :) = acf'*acf_comb;
        end;
    end;
%     
    % plot some example feature frames

    if find(k == sample_frames)
        figure;
        for c = 1:NUM_FEATURE_CHANNELS
            subplot(NUM_FEATURE_CHANNELS, 3, 3*(c-1)+1);
            plot(curr_time_axis, curr_df_frame(:, c));
            title(sprintf('Detection function: band=%d t=%3.2f s', c, (k-1)*FEATURE_HOP_SIZE/FEATURE_RATE))
            xlabel('Time (seconds)');
            ylabel('Energy');


            subplot(NUM_FEATURE_CHANNELS, 3, 3*(c-1)+2);
            plot(curr_time_axis, acf(:,c)); hold on;
            title(sprintf('ACF: band=%d t=%3.2f s', c, (k-1)*FEATURE_HOP_SIZE/FEATURE_RATE));
            xlabel('Lag (seconds)');
            ylabel('Similarity');
            stem(MIN_LAG, 1); stem(MAX_LAG, 1);

            estimates = tempo_estimates{k, c};
            confidences = tempo_confidences{k, c};

            if ~isempty(estimates)
                stem(curr_time_axis(estimates), confidences);
            end
            
            subplot(NUM_FEATURE_CHANNELS, 3, 3*(c-1)+3);
            plot((candidate_tempos/FEATURE_RATE), tempo_strength(:, c));
            title(sprintf('Tempo strength: band=%d', c));
            xlabel('Tempo (BPM)');
        end
    end;

end;
time=toc;
fprintf('spent %f s calculating autocorrelations\n', time);

%%%%%%%%%%%%%%%%%%%%
%% Compute phase estimates for each tempo hypothesis
% make sure to use the (differenced) DF
% following procedure in 'Context-Dependent Beat tracking'
%%%%%%%%%%%%%%%%%%%%%

tic;
for k=1:NUM_FEATURE_FRAMES
    % current detection funciton frame
    curr_df_frame = feature_frame(df, k); 
    % actual lag time, trimmed to length of current feature frame
    curr_time_axis = feature_time_axis(1:length(feature_frame(df, k)));
    
    ALIGNMENTS_PER_TEMPO = 3;
    % compute beat alignments for each tempo hypothesis in each channel,
    % pick ALIGNMENTS_PER_TEMPO highest peaks
       
    % plot phase alignments for different tempo hypothesis from this channel
    if find(k == sample_frames)
        figure;
        subplot(NUM_FEATURE_CHANNELS, NUM_TEMPO_ESTIMATES, 1);
    end;
    
    for c=1:NUM_FEATURE_CHANNELS
        % One iteration per frame per feature channel
        tempos = tempo_estimates{k, c};

        if isempty(tempos)
            continue;
        end;
        
        beat_alignment_estimates{k, c} = ...
            zeros(length(tempos), ALIGNMENTS_PER_TEMPO);
        beat_alignment_confidences{k, c} = ...
            zeros(length(tempos), ALIGNMENTS_PER_TEMPO);

        for tempo_index = 1:length(tempos)
            
            % in samples
            tempo_hypothesis = tempos(tempo_index);
            
            % make an impulse train for each tempo hypothesis,
            % then slide it along the detection function for one tempo period,
            % to find where it lines up the best
            % give feature as row vector, as it is flipped left-to-right
            alignment_function = ...
                beat_alignment_function(curr_df_frame(:, c)', tempo_hypothesis);

            % max distance between peaks depends on the tempo
            % assume that 'off by half a semiquaver' is the same beat
            % location
            [curr_alignment_confidences, curr_alignment_estimates]  = ...
                findpeaks(alignment_function, 1:tempo_hypothesis, ...
                    'MinPeakHeight', 0.01, ...
                    'MinPeakDistance', tempo_hypothesis/8, ...
                    'NPeaks', ALIGNMENTS_PER_TEMPO, ...
                    'SortStr', 'descend');


            % make sign negative to reflect that the offsets are calculated
            % backwards from the end of the frame

            curr_alignment_estimates = -1*curr_alignment_estimates;
            beat_locations = zeros(length(curr_df_frame(:, c)), ALIGNMENTS_PER_TEMPO);
            num_alignment_peaks = length(curr_alignment_confidences);

            % store calculated beat locations as impulses with height equal to
            % the alignment confidence, which can be overlaid on top of the feature
            % vector
            for alignment_index = 1:num_alignment_peaks
                curr_alignment_estimate = curr_alignment_estimates(alignment_index);
                curr_alignment_confidence = curr_alignment_confidences(alignment_index);
                % calculate beat locations by adding multiples of the tempo
                % hypothesis to the beat alignment estimate
                beat_locations(end+curr_alignment_estimate:-tempo_hypothesis:1, alignment_index) = ...
                    curr_alignment_confidence;
            end;

            if find(k == sample_frames)
                subplot(NUM_FEATURE_CHANNELS, NUM_TEMPO_ESTIMATES, ...
                    (c - 1)*NUM_TEMPO_ESTIMATES + tempo_index);
                
                plot(curr_time_axis, curr_df_frame(:, c)); hold on;
                ylim([0, max(curr_df_frame(:, c))]);
                nonzero_indices_1 = find(beat_locations(:, 1) ~=0);
                nonzero_indices_2 = find(beat_locations(:, 2) ~=0);
                nonzero_indices_3 = find(beat_locations(:, 3) ~=0); 
                stem(curr_time_axis(nonzero_indices_1)', beat_locations(nonzero_indices_1, 1), 'red');
                stem(curr_time_axis(nonzero_indices_2)', beat_locations(nonzero_indices_2, 2), 'blue');
                stem(curr_time_axis(nonzero_indices_3)', beat_locations(nonzero_indices_3, 3), 'green');
                %plot((1:tempo_hypothesis)/FEATURE_RATE, alignment_function);
                %hold on;
                title(sprintf('tempo = %2.3f, t = %3.2f s, band=%d', ...
                    tempo_hypothesis/FEATURE_RATE, ...
                    (k-1)*FEATURE_HOP_SIZE/FEATURE_RATE, ...
                    c)); 
                xlabel('Offset from end of frame (seconds)');
                ylabel('Strength');
            end;
            
            beat_alignment_confidences{k, c}(tempo_index, 1:num_alignment_peaks) = ...
                curr_alignment_confidences;
            beat_alignment_estimates{k, c}(tempo_index, 1:num_alignment_peaks) = ...
                curr_alignment_estimates;
        end;
    end;

end;

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

    for k = 1:NUM_FEATURE_FRAMES
        % set time for estimate estimates from frame k 
        % to be at the end of that frame
        % so first frame's estimates correspond to FEATURE_TIME
        frame_time = FEATURE_WIN_TIME + (k-1)/TEMPO_UPDATE_FREQUENCY;
        for tempo_index = 1:length(tempo_estimates{k, c})
            % tempo and alignment are expresseed in samples,
            % convert to time by dividing by feature rate
            curr_tempo = tempo_estimates{k, c}(tempo_index);
            curr_tempo_confidence = tempo_confidences{k, c}(tempo_index);
            curr_alignments = beat_alignment_estimates{k, c}(tempo_index, :);
            curr_alignment_confidences = beat_alignment_confidences{k, c}(tempo_index, :);

            % output one line per alignment estimate per tempo estimate
            % per frame per channel!
            for alignment_index = 1:length(curr_alignments)
                curr_alignment = curr_alignments(alignment_index);
                curr_alignment_confidence = curr_alignment_confidences(alignment_index);
                fprintf(outfile, strcat( ...
                    'frame=%d\t', 'time=%f\t', ...
                    'tempo=%f\t', 'tempo_confidence=%f\t', ...
                    'alignment=%f\t', 'alignment_confidence=%f', ...
                    '\n'), ...
                    k, frame_time, ...
                    curr_tempo/FEATURE_RATE, curr_tempo_confidence, ...
                    curr_alignment/FEATURE_RATE, curr_alignment_confidence);
            end;
        end
    end;
end;

time = toc;
fprintf('spent %f s writing beat data\n', time);

% vim: tabstop=4 expandtab
