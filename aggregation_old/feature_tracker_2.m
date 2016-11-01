% feature_tracker.m
% Author: Max Fisher

% Initial implementation of a simplified real-time beat tracker, which
% relies on a single audio feature, with the goal of integrating it into
% the Beaker ensemble learning framework [1].
% We use the Log Filtered Spectral Flux onset detection as proposed in [2],
% which takes ideas from [4] as well as other papers,
% and follow the general method outlined in [3] for tempo estimation.

% [1] Daniels, M. (2015).
% An Ensemble Framework for Real-Time Audio Beat Tracking.
% PhD dissertation, University of California, San Diego.

% [2] Böck, S. Krebs, F. Schedl, M. (2012).
% Evaluating the Online Capabilities of Onset Detection Methods.
% Published by the International Society of Music Information Retrieval
% (ISMIR), 2012.

% [3] Percival, G. Tzanetakis, G. (2014).
% Streamlined Tempo Estimation Based on Autocorrelation and Cross-Correlation
% With Pulses.
% IEEE/ACM Transactions on Audio, Speech and Language Processing,
% Vol. 22, No. 12, December 2014

% [4] Klapuri A., Eronen A., Astola, J. (2006).
% Analysis of the Metre of Acoustic Musical Signals.
% IEEE Transactions on Audio, Speech and Language Processing,
% Vol. 14, No. 1, January 2006

% Problems: ACF is like comb filtering with sine waves - fast but they're
% too 'blunt' to be useful -> need more spiky functions. (i.e. custom comb)
% DFT evaluates all frequencies equally spaced -> doesn't affect how
% hearing works. We need to have more information at lower frequencies
% -> constant-Q transform, or ERB bands. 

% For now; just sum up DFT Bins logarithmically.


%%%%%%%%%%%%%%%
%% (independently) tweakable parametersuntitled
%%%%%%%%%%%%%%%

% window parameters, in samples
% approximate time of audio to use for windowing.
% Window size in samples will be rounded to the nearest
% power of 2.
WAV_EXT = '.wav';
    

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
NUM_FEATURE_CHANNELS = 5;

% BPM ranges to allow when detecting tempo
MIN_BPM = 40;
MAX_BPM = 240;

% BPM ranges translated to autocorrelation lag
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
FEATURE_INTERP = 1;
FEATURE_RATE = FEATURE_INTERP*Fs/AUDIO_HOP_SIZE; % times 2 after upsampling

% how many past feature calculation samples to use when estimating current tempo
FEATURE_WIN_LENGTH = 2^round(log2(FEATURE_TIME_APPROX*Fs/AUDIO_HOP_SIZE));
% actual length of audio used, as corresponding to the FEATURE_WIN_LENGTH
FEATURE_WIN_TIME = FEATURE_WIN_LENGTH*AUDIO_HOP_SIZE/Fs;

%FEATURE_HOP_SIZE = 2^round(log2(FEATURE_RATE/TEMPO_UPDATE_FREQUENCY_APPROX));
FEATURE_HOP_SIZE = FEATURE_WIN_LENGTH*(1-FEATURE_WIN_OVERLAP_PERCENT/100);
TEMPO_UPDATE_FREQUENCY = FEATURE_RATE/FEATURE_HOP_SIZE; % Hz

NUM_FEATURE_FRAMES = ceil(NUM_AUDIO_FRAMES/FEATURE_HOP_SIZE);

disp('Audio sample rate: '); disp(Fs);
disp('Audio window time (seconds): '); disp(AUDIO_WIN_TIME);
disp('Audio window length (samples): '); disp(AUDIO_WIN_LENGTH);
disp('Audio hop size (samples): '); disp(AUDIO_HOP_SIZE);
disp('Number of audio frames: '); disp(NUM_AUDIO_FRAMES);

disp('Feature sample rate: '); disp(FEATURE_RATE);
disp('Feature window time (seconds): '); disp(FEATURE_WIN_TIME);
disp('Feature window length (samples): '); disp(FEATURE_WIN_LENGTH);
disp('Feature hop size (samples): '); disp(FEATURE_HOP_SIZE);
disp('Number of feature frames: '); disp(NUM_FEATURE_FRAMES);

disp('Tempo Update Frequency (Hz):'); disp(TEMPO_UPDATE_FREQUENCY);





%%%%%%%%%%%%%%
%% algorithm setup and variable initialisation
%%%%%%%%%%%%%%


%y[n] is the audio sample
%f[n] is the feature vector

% for now, store all the feature calculations in f[k, n], but only use the last
% FEATURE_WIN_LENGTH many to compute tempo at each time step
f = zeros(NUM_AUDIO_FRAMES, NUM_FEATURE_CHANNELS);
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
LP_cutoff_freq = 10; %Hz
% 1 corresponds to the maximum frequency 'sample rate' of pi == 1/(2*df)
% ......not sure about this part
LP_cutoff_freq_normalised = LP_cutoff_freq/(delta_f/2);
LP_order = 5; % was six, but halved in order to do reversed filtering
[LP_num, LP_denom] = butter(LP_order, LP_cutoff_freq_normalised);
% normalise by ratio of frequencies (again following [4])
LP_num = LP_num/LP_cutoff_freq_normalised;



feature_frame = @(g, n) g((n-1)*FEATURE_HOP_SIZE + 1: min(((n-1)*FEATURE_HOP_SIZE + FEATURE_WIN_LENGTH), end), :);

feature_time_axis = (1:FEATURE_WIN_LENGTH)/FEATURE_RATE;

sample_frames = round(linspace(1, NUM_FEATURE_FRAMES, 4));
sample_frames = sample_frames(2:end-1);

%sample_frames = [sample_frames ; sample_frames + 1];

% use to disable plotting;
sample_frames = [];

%%%%%%%%%%%%%%%%%
%% Processing of raw audio frames
%%%%%%%%%%%%%%%%%
band = cell(1, NUM_FEATURE_CHANNELS);
band_ends = [1, round(2.^linspace(round(log2(AUDIO_WIN_LENGTH)/4), log2(AUDIO_WIN_LENGTH/2), NUM_FEATURE_CHANNELS))];
disp('Audio frequency bands (DFT indexes):');
disp(band_ends);
for c = 1:NUM_FEATURE_CHANNELS
    % this excludes the lowest DFT band (i.e. DC up to Fs/NFFT Hz)
    band{1, c} = (band_ends(c)+1):band_ends(c+1);
end;
tic;
for k = 1:NUM_AUDIO_FRAMES
    % last frames will not be complete; if so, truncate the window
    curr_frame = audio_frame(k).*win(1:length(audio_frame(k)));
    mean = sum(curr_frame)/length(curr_frame);
    variance = sum(curr_frame.^2)/(length(curr_frame) - 1);
    curr_frame = (curr_frame  - mean)/sqrt(variance);
    
    % calculate log magnitude energy spectrum with parameter mu = 1000
    % (idea introduced in [4], value for mu from [3])
    mu = 1000;

    % following steps in [4] we interpolate the log power spectrum by a factor of 2,
    % (restoring the length of the log spectrum vector to AUDIO_WIN_LENGTH)
    % and low pass filter with a butterworth filter: 6th order; cutoff 10Hz.
    %filtered_power_spectrum = filter(LP_num, LP_denom, upsample(log_power_spectrum, 2));
    % jokes... LPF doesn't seem to to do much...
    
    
    % in the case of incomplete frames, pad out to AUDIO_WIN_LENGTH samples
    % use double length FFT and then throw away half the spectrum,
    % plus DC level?
    X = fft(curr_frame, AUDIO_WIN_LENGTH);
    log_power_spectrum = log(1 + mu*abs(X(1:AUDIO_WIN_LENGTH/2)))/log(1 + mu);

    % half-wave rectified difference between this and prev frame's
    % log power spectrum, to create a feature that peaks at the onset time
    if k > 1
        % second difference so that we take the difference with a 50%
        % overlapping window, but have double the feature rate
        unrectified_difference = log_power_spectrum - prev_log_power_spectrum;
        HWR_difference = unrectified_difference;
        negative_entries = find(HWR_difference < 0);
        HWR_difference(negative_entries) = 0;
    end
  
    %prev2_log_power_spectrum = prev_log_power_spectrum;
    prev_log_power_spectrum = log_power_spectrum;


    % now create 4 features based on summing the flux in adjacent bins, so
    % that we can do analysis in different frequency channels
    % use increasing number of bins as we increase in frequency, as this
    % vaguely reflects human hearing being logarithmic in pitch, but more
    % linear at lowest frequencies (ERB/critical Bands)
    % "Poor man's constant Q transform"
    % do we really need to normalise by length?
    for c = 1:NUM_FEATURE_CHANNELS
        f(k, c) = sum(log_power_spectrum(band{c}))/length(band{c});
        % weighted sum the result with the undifferentiated log_power_spectrum
        % makes analysis better according to [4]
        lambda = 0.8;
        if k > 2
            df(k, c) = sum((1 - lambda)*log_power_spectrum(band{c}) + ...
                lambda*HWR_difference(band{c}))/length(band{c});
        end;
    end
    
%     if find(k == sample_frames)
%         figure;
%         subplot(2, 1, 1);
%         plot(log_power_spectrum);
%         title(sprintf('log power spectrum, k = %d', k));
%         subplot(2, 1, 2);
%         plot(HWR_difference);
%         title(sprintf('HWR difference of log power spectrum, k = %d', k));
%     end;
    %f(k) = sum(log_power_spectrum);

end;
time = toc;
fprintf('spent %f s calculating feature\n', time);

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

% try upsampling feature?
% low pass
%f = filter(LP_num, LP_denom, upsample(f, 2));
% then filter in reverse to remove phase distortion
%f = filter(LP_num, LP_denom, flipud(f));
% if FEATURE_INTERP > 1
%     f_interp = interp(f, FEATURE_INTERP);
%     f = f_interp;
% end;

% beat hypothesis data
% the nth index of the array in each cell corresponds to a single
% (tempo, tempo confidence, beat location, beat location confidence) tuple
% which forms a complete hypothesis of a particular feature at a particular
% frame

tempo_estimates = cell(NUM_FEATURE_FRAMES, NUM_FEATURE_CHANNELS);
tempo_confidences = cell(NUM_FEATURE_FRAMES, NUM_FEATURE_CHANNELS);
beat_alignment_estimates = cell(NUM_FEATURE_FRAMES, NUM_FEATURE_CHANNELS);
beat_alignment_confidences = cell(NUM_FEATURE_FRAMES, NUM_FEATURE_CHANNELS);


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
        ac = autocorrelation_2(curr_df_frame(:, c));
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
                'NPeaks', 4, ...
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
            title(sprintf('Detection function: freq band=%d t=%3.1f s', c, (k-1)*FEATURE_HOP_SIZE/FEATURE_RATE))
            xlabel('Time (seconds)');
            ylabel('Energy');


            subplot(NUM_FEATURE_CHANNELS, 3, 3*(c-1)+2);
            plot(curr_time_axis, acf(:,c)); hold on;
            title(sprintf('Autocorr: freq band=%d t=%3.1f s', c, (k-1)*FEATURE_HOP_SIZE/FEATURE_RATE));
            xlabel('Lag (seconds)');
            ylabel('Similarity');
            stem(MIN_LAG, 1); stem(MAX_LAG, 1);

            estimates = tempo_estimates{k, c};
            confidences = tempo_confidences{k, c};

            if ~isempty(estimates)
                stem(curr_time_axis(estimates), confidences);
            end
            
            subplot(NUM_FEATURE_CHANNELS, 3, 3*(c-1)+3);
            plot(60./(candidate_tempos/FEATURE_RATE), tempo_strength(:, c));
            title(sprintf('Tempo strength: freq band=%d', c));
            xlabel('Tempo (BPM)');
            ylabel('Strength');
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
    for c=1:NUM_FEATURE_CHANNELS
        %% One iteration per frame per feature channel
        tempos = tempo_estimates{k, c};

        if isempty(tempos)
            continue;
        end;
        
        beat_alignment_estimates{k, c} = ...
            zeros(length(tempos), ALIGNMENTS_PER_TEMPO);
        beat_alignment_confidences{k, c} = ...
            zeros(length(tempos), ALIGNMENTS_PER_TEMPO);


        
        % plot phase alignments for different tempo hypothesis from this channel
        if find(k == sample_frames)
            figure;
            subplot(1+length(tempos), 1, 1);
            plot(curr_time_axis, curr_df_frame(:, c));
            title(sprintf('Log power: freq band=%d t=%3.1f s', c, (k-1)*FEATURE_HOP_SIZE/FEATURE_RATE))
            xlabel('Time (seconds)');
            ylabel('Energy in band');
        end;
        
        for tempo_index = 1:length(tempos)
            
            % in samples
            tempo_hypothesis = tempos(tempo_index);
            
            % make an impulse train for each tempo hypothesis,
            % then slide it along the detection function for one tempo period,
            % to find where it lines up the best
            % give feature as row vector, as it is flipped left-to-right
            alignment_function = ...
                beat_alignment_function(curr_df_frame(:, c)', tempo_hypothesis);

            [curr_alignment_confidences, curr_alignment_estimates]  = ...
                findpeaks(alignment_function, 1:tempo_hypothesis, ...
                    'MinPeakHeight', 0.01, ...
                    'MinPeakDistance', MIN_LAG/2, ...
                    'NPeaks', ALIGNMENTS_PER_TEMPO, ...
                    'SortStr', 'descend');
            beat_locations = zeros(length(curr_df_frame(:, c)), ALIGNMENTS_PER_TEMPO);
            num_alignment_peaks = length(curr_alignment_confidences);
            
            % store calculated beat locations as impulses with height equal to
            % the alignment confidence, which can be overlaid on top of the feature
            % vector
            for alignment_index = 1:num_alignment_peaks
                curr_alignment_estimate = curr_alignment_estimates(alignment_index);
                curr_alignment_confidence = curr_alignment_confidences(alignment_index);
                beat_locations(end-curr_alignment_estimate:-tempo_hypothesis:1, alignment_index) = ...
                    curr_alignment_confidence;
            end
                        
            if find(k == sample_frames)
                subplot(1+length(tempos), 1, 1+tempo_index);
                plot(curr_time_axis, curr_df_frame(:, c)); hold on;
                ylim([min(curr_df_frame(:, c)), max(curr_df_frame(:, c))]);
                stem(curr_time_axis, beat_locations(:, 1), 'red');
                stem(curr_time_axis, beat_locations(:, 2), 'blue');
                %plot((1:tempo_sample_spacing)/FEATURE_RATE, alignment_function);
                %hold on;
                title(sprintf('Beat alignment: tempo = %f', 60/(tempo_hypothesis/FEATURE_RATE))); 
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
feature_tempos = {0,0,0,0,0};
feature_phases = {0,0,0,0,0};

for c = 1:NUM_FEATURE_CHANNELS;
    %filename_for_channel = strcat('./', WAV_NAME, BEAT_DATA_OUTPUT_SUFFIX, ...
    %    sprintf('-channel-%d.txt', c));
    %outfile = fopen(filename_for_channel, 'w+');

    tempos = {};
    phases = {};
    t_confs = {};
    p_confs = {};
    
    for k = 1:NUM_FEATURE_FRAMES
        % set time for estimate estimates from frame k 
        % to be at the end of that frame
        % so first frame's estimates correspond to FEATURE_TIME
        frame_time = FEATURE_WIN_TIME + (k-1)/TEMPO_UPDATE_FREQUENCY;
        tempo = [];
        phase = [];
        tempo_confidence = [];
        phase_confidence = [];
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
                %fprintf(outfile, strcat( ...
                %    'frame=%d\t', 'time=%f\t', ...
                %    'tempo=%f\t', 'tempo_confidence=%f\t', ...
                %    'alignment=%f\t', 'alignment_confidence=%f', ...
                %    '\n'), ...
                %    k, frame_time, ...
                %    curr_tempo/FEATURE_RATE, curr_tempo_confidence, ...
                %    curr_alignment/FEATURE_RATE, curr_alignment_confidence);
                tempo = [tempo, curr_tempo/FEATURE_RATE];
                tempo_confidence = [tempo_confidence, curr_tempo_confidence];
                phase = [phase, curr_alignment/FEATURE_RATE];
                phase_confidence = [phase_confidence, curr_alignment_confidence];
            end;
        end
        tempos = [tempos, tempo];
        phases = [phases, phase];
        t_confs = [t_confs, tempo_confidence];
        p_confs = [p_confs, phase_confidence];
    end;
    feature_tempos(c) = {tempos};
    feature_phases(c) = {phases};
    feature_t_confs(c) = {t_confs};
    feature_p_confs(c) = {p_confs};
end;

time = toc;
fprintf('spent %f s writing beat data\n', time);

% vim: tabstop=4 expandtab
