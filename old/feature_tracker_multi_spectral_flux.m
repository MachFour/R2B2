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
%% (independently) tweakable parameters
%%%%%%%%%%%%%%%
WAV_NAME = 'a-new-sky-presets';
WAV_EXT = '.wav';

% window parameters, in samples
AUDIO_WIN_LENGTH = 1024;
% number of samples advanced by each window
AUDIO_HOP_SIZE = 512;
AUDIO_WIN_TYPE = 'Hann';

% feature parameters

% minimum length of audio to use for tempo estimation at each time point (seconds)
% actual length of feature frames will be power of 2 closest to this.
FEATURE_TIME_APPROX = 5; %seconds

% Approximately how often to make estimates of new tempo by recalculating
% the autocorrelation. This is used to calculate the feature hop size
TEMPO_UPDATE_FREQUENCY_APPROX = 2; % Hz

% or maybe use a window that weights recent samples more than older (by a few seconds) samples 
FEATURE_WIN_TYPE = 'Rect';

% how many bands to use for periodicity analysis
FEATURE_CHANNELS = 4;

% BPM ranges to allow when detecting tempo
MIN_TEMPO = 25;
MAX_TEMPO = 400;

% BPM ranges translated to autocorrelation lag
MAX_LAG = 60/MIN_TEMPO;
MIN_LAG = 60/MAX_TEMPO;


% output the feature values to a time-value CSV file with the following suffix
FEATURE_OUTPUT_SUFFIX = '-feature';

%%%%%%%%%%%%%%%%%%%%%%%%
%% setup and calculated parameters
%%%%%%%%%%%%%%%%%%%%%%%%

% assume mono audio wav file in current dir.


[y, Fs] = audioread(strcat('./', WAV_NAME, WAV_EXT));

NUM_AUDIO_FRAMES = ceil(length(y)/AUDIO_HOP_SIZE);

% sampling rate of the feature measurement
FEATURE_RATE = 2*Fs/AUDIO_HOP_SIZE; % times 2 after upsampling

% how many past feature calculation samples to use when estimating current tempo
FEATURE_WIN_LENGTH = 2^round(log2(FEATURE_TIME_APPROX*Fs/AUDIO_HOP_SIZE));
% actual length of audio used, as corresponding to the FEATURE_WIN_LENGTH
FEATURE_TIME = FEATURE_WIN_LENGTH*AUDIO_HOP_SIZE/Fs;

FEATURE_HOP_SIZE = 2^round(log2(FEATURE_RATE/TEMPO_UPDATE_FREQUENCY_APPROX));
TEMPO_UPDATE_FREQUENCY = FEATURE_RATE/FEATURE_HOP_SIZE; % Hz

NUM_FEATURE_FRAMES = ceil(NUM_AUDIO_FRAMES/FEATURE_HOP_SIZE);


%y[n] is the audio sample
%f[n] is the feature vector

% for now, store all the feature calculations, but only use the last
% FEATURE_LENGTH many to compute tempo at each time step
f = zeros(NUM_AUDIO_FRAMES, 1);

% multi channel feature
f_multi = zeros(NUM_AUDIO_FRAMES, 4);

% choice of windowing function.
win = window(@hann, AUDIO_WIN_LENGTH);

%log_power_spectrum = zeros(AUDIO_WIN_LENGTH);

% handle to get the correct set of samples for frame n; 1-indexed.
frame = @(n) y(((n-1)*AUDIO_HOP_SIZE + 1): min(((n-1)*AUDIO_HOP_SIZE + AUDIO_WIN_LENGTH), end), 1);
% handle to get sample numbers for a given frame

% low pass filter for filtering the log magnitude spectrum
% [4] specifies a 10Hz cutoff frequency for the filter, but we're filtering the spectrum...
% what's the sampling rate?
% -> related to the frequency resolution
df = Fs/AUDIO_WIN_LENGTH; %  ~= 44Hz -> this is the distance (in Hz frequency) between DFT samples
% -> thus the frequency 'sample rate' is 1/44 (Hz)^-1 (i.e. seconds!)
% -> 44Hz! ...  does this mean that frequency resolution is worse than an entire octave 
% at the lower end of hearing?? Yes -> we should really use the Constant-Q transform here
% -> Or filter with logarithmically spaced bandpass filters, like Klapuri
LP_cutoff_freq = 10; %Hz
% 1 corresponds to the maximum frequency 'sample rate' of pi == 1/(2*df)
% ......not sure about this part
LP_cutoff_freq_normalised = LP_cutoff_freq/(df/2);
LP_order = 6;
[LP_num, LP_denom] = butter(LP_order, LP_cutoff_freq_normalised);
% normalise by ratio of frequencies (again following [4])
LP_num = LP_num/LP_cutoff_freq_normalised;

%%%%%%%%%%%%%%%%%
%% Processing of raw audio frames
%%%%%%%%%%%%%%%%%

for k = 1:NUM_AUDIO_FRAMES
    % last frames will not be complete; if so, truncate the window
    curr_frame = frame(k).*win(1:length(frame(k)));

    % normalise and subtract mean? -> no, do it for the ACF estimation

    % calculate log magnitude energy spectrum with parameter mu = 1000
    % (idea introduced in [4], value for mu from [3])
    mu = 1000;

    % following steps in [4] we interpolate the log power spectrum by a factor of 2,
    % (restoring the length of the log spectrum vector to AUDIO_WIN_LENGTH)
    % and low pass filter with a butterworth filter: 6th order; cutoff 10Hz.
    %filtered_power_spectrum = filter(LP_num, LP_denom, upsample(log_power_spectrum, 2));

    % LPF doesn't seem to to do much...
    % in the case of incomplete frames, pad out to AUDIO_WIN_LENGTH samples
    % use double length FFT and then throw away half the spectrum,
    % plus DC level?
    X = fft(curr_frame, 2*AUDIO_WIN_LENGTH);
    log_power_spectrum = log(1 + mu*abs(X(1:AUDIO_WIN_LENGTH)))/log(1 + mu);

    % now create 4 features based on summing the flux in adjacent bins, so
    % that we can do analysis in different frequency channels
    % use increasing number of bins as we increase in frequency, as this
    % vaguely reflects human hearing being logarithmic in pitch, but more
    % linear at lowest frequencies (ERB Bands)
    % "Poor man's constant Q transform"
    % >> round(2.^linspace(0, 10, 5))
    % ans =  1          6           32          181        1024
    % >> round(2.^linspace(6, 10, 5))
    % ans =  64         128         256         512        1024
    band1 = 1:AUDIO_WIN_LENGTH/8;
    band2 = AUDIO_WIN_LENGTH/8+1:AUDIO_WIN_LENGTH/4;
    band3 = AUDIO_WIN_LENGTH/4+1:AUDIO_WIN_LENGTH/2;
    band4 = AUDIO_WIN_LENGTH/2+1:AUDIO_WIN_LENGTH;
    f_multi(k, 1) = sum(log_power_spectrum(band1))/length(band1);
    f_multi(k, 2) = sum(log_power_spectrum(band2))/length(band2);
    f_multi(k, 3) = sum(log_power_spectrum(band3))/length(band3);
    f_multi(k, 4) = sum(log_power_spectrum(band4))/length(band4);
    f(k) = sum(log_power_spectrum);

end;

% try upsampling feature?
% low pass
f = filter(LP_num, LP_denom, upsample(f, 2));
f_multi = filter(LP_num, LP_denom, upsample(f_multi, 2));

%%%%%%%%%%%%%%%%%%%%
%% Write feature data to file
%%%%%%%%%%%%%%%%%%%%%

% export f to time instants to import into Sonic Visualiser
% f[k] represents the difference between frame k and frame k-1.
% However the frames overlap - at what point should this measurement be
% attributed to?
% in terms of time, we say that this value corresponds to the end of frame k-1.
% Then we put f(k) at (k-1)/FEATURE_RATE seconds

outfile_0 = fopen(strcat('./', WAV_NAME, FEATURE_OUTPUT_SUFFIX, '-0.txt'), 'w+');
outfile_1 = fopen(strcat('./', WAV_NAME, FEATURE_OUTPUT_SUFFIX, '-1.txt'), 'w+');
outfile_2 = fopen(strcat('./', WAV_NAME, FEATURE_OUTPUT_SUFFIX, '-2.txt'), 'w+');
outfile_3 = fopen(strcat('./', WAV_NAME, FEATURE_OUTPUT_SUFFIX, '-3.txt'), 'w+');
outfile_4 = fopen(strcat('./', WAV_NAME, FEATURE_OUTPUT_SUFFIX, '-4.txt'), 'w+');
 

single_channel_outfile = fopen(strcat('./', WAV_NAME, FEATURE_OUTPUT_SUFFIX, '.txt'), 'w+');
for k = 1:length(f)
    fprintf(outfile_0, '%f\t%f\n', (k-1)/FEATURE_RATE, f(k));
    fprintf(outfile_1, '%f\t%f\n', (k-1)/FEATURE_RATE, f_multi(k, 1));
    fprintf(outfile_2, '%f\t%f\n', (k-1)/FEATURE_RATE, f_multi(k, 2));
    fprintf(outfile_3, '%f\t%f\n', (k-1)/FEATURE_RATE, f_multi(k, 3));
    fprintf(outfile_4, '%f\t%f\n', (k-1)/FEATURE_RATE, f_multi(k, 4));
end;

%%%%%%%%%%%%%%%%%%%%
%% Processing of feature vector for periodicity
% In reality everything would be done in the previous loop,
% also tempo estimates would have to be made for past feature measurements
%%%%%%%%%%%%%%%%%%%%%%

tempo_peak = zeros(NUM_FEATURE_FRAMES, 1);
tempo_confidence = zeros(NUM_FEATURE_FRAMES, 1);

feature_frame = @(n) f((n-1)*FEATURE_HOP_SIZE + 1: min(((n-1)*FEATURE_HOP_SIZE + FEATURE_WIN_LENGTH), end), 1);
feature_frame_multichannel = @(n) f_multi((n-1)*FEATURE_HOP_SIZE + 1: min(((n-1)*FEATURE_HOP_SIZE + FEATURE_WIN_LENGTH), end), :);

feature_time_axis = (1:FEATURE_WIN_LENGTH)/FEATURE_RATE;

%for k = 1:NUM_FEATURE_FRAMES
for k = round(linspace(1, NUM_FEATURE_FRAMES, 5))
    % go through and calculate autocorrelations for each slice of FEATURE_LENGTH features
    % no windowing so far
    curr_feature_frame_channels = feature_frame_multichannel(k); % no windowing so far
    ac_total = autocorrelation(feature_frame(k));
    ac_1 = autocorrelation(curr_feature_frame_channels(:, 1));
    ac_2 = autocorrelation(curr_feature_frame_channels(:, 2));
    ac_3 = autocorrelation(curr_feature_frame_channels(:, 3));
    ac_4 = autocorrelation(curr_feature_frame_channels(:, 4));


    % plot some example feature frames
    if find(k == round(linspace(1, NUM_FEATURE_FRAMES, 5)))
        figure;

        curr_time_axis = feature_time_axis(1:length(feature_frame(k)));
        curr_fk = feature_frame_multichannel(k);

        subplot(5, 2, 1); plot(curr_time_axis, feature_frame(k)); 
        title(sprintf('Total Summed LPS at feature frame %d', k))
        subplot(5, 2, 3); plot(curr_time_axis, curr_fk(:, 1));
        title(sprintf('Band 1 Summed LPS at feature frame %d', k))
        subplot(5, 2, 5); plot(curr_time_axis, curr_fk(:, 2));
        title(sprintf('Band 2 Summed LPS at feature frame %d', k))
        subplot(5, 2, 7); plot(curr_time_axis, curr_fk(:, 3)); 
        title(sprintf('Band 3 Summed LPS at feature frame %d', k))
        subplot(5, 2, 9); plot(curr_time_axis, curr_fk(:, 4));
        title(sprintf('Band 4 Summed LPS at feature frame %d', k))

        subplot(5, 2, 2); plot(curr_time_axis, ac_total);
        title(sprintf('Total autocorrelation for feature frame %d', k));
        hold on ; stem(MIN_LAG, 1); stem(MAX_LAG, 1);

        subplot(5, 2, 4); plot(curr_time_axis, ac_1);
        title(sprintf('Band 1 autocorrelation for feature frame %d', k));
        hold on ; stem(MIN_LAG, 1); stem(MAX_LAG, 1);

        subplot(5, 2, 6); plot(curr_time_axis, ac_2);
        title(sprintf('Band 2 autocorrelation for feature frame %d', k));
        hold on ; stem(MIN_LAG, 1); stem(MAX_LAG, 1);

        subplot(5, 2, 8); plot(curr_time_axis, ac_3);
        title(sprintf('Band 3 autocorrelation for feature frame %d', k));
        hold on ; stem(MIN_LAG, 1); stem(MAX_LAG, 1);

        subplot(5, 2, 10); plot(curr_time_axis, ac_4);
        title(sprintf('Band 4 autocorrelation for feature frame %d', k));
        hold on ; stem(MIN_LAG, 1); stem(MAX_LAG, 1);
    end

    % normalise to have area?? 1
    %tempo_peaks = find(ac == max(ac));
    %tempo_confidences = max(ac);

    % find peaks until confidence is within a certain ratio of the highest
    % ac(tempo_peak) = 0;
    % tempo_peak = find(ac == max(ac)) -> find second peak, make sure it's far enough away
end

% vim: tabstop=4 expandtab


