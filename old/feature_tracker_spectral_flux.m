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


%%%%%%%%%%%%%%%
%% (independently) tweakable parameters
%%%%%%%%%%%%%%%
WAV_NAME = 'walking-short';
WAV_EXT = '.wav';

% window parameters, in samples
AUDIO_WIN_LENGTH = 1024;
% number of samples advanced by each window
AUDIO_HOP_SIZE = 512;
AUDIO_WIN_TYPE = 'Hann';

% feature parameters
% how many past feature calculation samples to use when estimating current tempo
FEATURE_WIN_LENGTH = 1024; % half that of [3]
FEATURE_HOP_SIZE = 64; % half that of [3] (to match win length)
% or maybe use a window that weights recent samples more than older (by a few seconds) samples 
FEATURE_WIN_TYPE = 'Rect';

% BPM ranges to allow when detecting tempo
MIN_TEMPO = 40;
MAX_TEMPO = 280;

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
FEATURE_RATE = Fs/AUDIO_HOP_SIZE;

% what length of audio to use for tempo estimation at each time point (seconds)
FEATURE_TIME = FEATURE_WIN_LENGTH*AUDIO_HOP_SIZE/Fs;

% How often to make estimates of new tempo by recalculating the autocorrelation
TEMPO_UPDATE_FREQUENCY = FEATURE_RATE/FEATURE_HOP_SIZE; % Hz
NUM_FEATURE_FRAMES = ceil(NUM_AUDIO_FRAMES/FEATURE_HOP_SIZE);


%y[n] is the audio sample
%f[n] is the feature vector

% for now, store all the feature calculations, but only use the last
% FEATURE_LENGTH many to compute tempo at each time step
f = zeros(NUM_AUDIO_FRAMES, 1);

% choice of windowing function.
win = window(@hann, AUDIO_WIN_LENGTH);

prev_log_power_spectrum = zeros(AUDIO_WIN_LENGTH, 1);

% handle to get the correct set of samples for frame n; 1-indexed.
frame = @(n) y(((n-1)*AUDIO_HOP_SIZE + 1): min(((n-1)*AUDIO_HOP_SIZE + AUDIO_WIN_LENGTH), end), 1);

% low pass filter for filtering the log magnitude spectrum
% [4] specifies a 10Hz cutoff frequency for the filter, but we're filtering the spectrum...
% what's the sampling rate?
% -> related to the frequency resolution
df = Fs/AUDIO_WIN_LENGTH; %  ~= 44Hz -> this is the distance (in Hz frequency) between DFT samples
% -> thus the frequency 'sample rate' is 1/44 (Hz)^-1 (i.e. seconds!)
% -> 44Hz! ...  does this mean that frequency resolution is worse than an entire octave 
% at the lower end of hearing?? Yes -> we should really use the Constant-Q transform here
% -> Or filter with logarithmically spaced bandpass filters, like Klapuri
LP_cutoff_freq = 4; %Hz
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

    % calculate log magnitude energy spectrum with parameter mu = 1000
    % (idea introduced in [4], value for mu from [3])
    mu = 1000;

    % in the case of incomplete frames, pad out to AUDIO_WIN_LENGTH samples
    %X = fft(curr_frame, AUDIO_WIN_LENGTH);
    %log_power_spectrum = log(1 + mu*abs(X(1:AUDIO_WIN_LENGTH/2)))/log(1 + mu);

    % following steps in [4] we interpolate the log power spectrum by a factor of 2,
    % (restoring the length of the log spectrum vector to AUDIO_WIN_LENGTH)
    % and low pass filter with a butterworth filter: 6th order; cutoff 10Hz.
    %filtered_power_spectrum = filter(LP_num, LP_denom, upsample(log_power_spectrum, 2));

    % LPF doesn't seem to to do much...
    % use double length FFT and then throw away half the spectrum
    X = fft(curr_frame, 2*AUDIO_WIN_LENGTH);
    log_power_spectrum = log(1 + mu*abs(X(1:AUDIO_WIN_LENGTH)))/log(1 + mu);
    filtered_power_spectrum = log_power_spectrum;


    % half-wave rectified difference between this and prev frame's (interpolated)
    % log power spectrum
    unrectified_diff = filtered_power_spectrum - prev_log_power_spectrum;

    HWR_diff = unrectified_diff;
    negative_entries = find(HWR_diff(HWR_diff < 0));
    HWR_diff(negative_entries) = 0;

    if mod(k, 1000) == 0
        figure;
        subplot (3, 1, 1);
        plot(filtered_power_spectrum);
        title(sprintf('Filtered Power spectrum  %d', k));
        subplot(3, 1, 2);
        plot(unrectified_diff);
        title(sprintf('Unrectified Difference %d', k));
    end

    % save current power spectrum for next time
    prev_log_power_spectrum = filtered_power_spectrum;

    % weighted sum the result with the undifferentiated log_power_spectrum
    % makes analysis better according to [4]
    lambda = 0.8;
    spectral_flux = sum((1 - lambda)*filtered_power_spectrum + lambda*HWR_diff);

    % since this feature involves a differential, it won't make much sense for k = 1
    f(k) = spectral_flux;
end;

% export f to time instants to import into Sonic Visualiser
% f[k] represents the difference between frame k and frame k-1.
% However the frames overlap - at what point should this measurement be
% attributed to?
% in terms of time, we say that this value corresponds to the end of frame k-1.
% Then we put f(k) at (k-1)/FEATURE_RATE seconds

outfile = fopen(strcat('./', WAV_NAME, FEATURE_OUTPUT_SUFFIX, '.txt'), 'w+');
for k = 1:length(f)
    fprintf(outfile, '%f\t%f\n', (k-1)/FEATURE_RATE, f(k));
end;

%%%%%%%%%%%%%%%%%%%%
%% Processing of feature vector for periodicity
% In reality everything would be done in the previous loop,
% also tempo estimates would have to be made for past feature measurements
%%%%%%%%%%%%%%%%%%%%%%


% try upsampling feature?
%filtered_f = filter(LP_num, LP_denom, upsample(f, 2));

tempo_peak = zeros(NUM_FEATURE_FRAMES, 1);
tempo_confidence = zeros(NUM_FEATURE_FRAMES, 1);

feature_frame = @(n) f((n-1)*FEATURE_HOP_SIZE + 1: min(((n-1)*FEATURE_HOP_SIZE + FEATURE_WIN_LENGTH), end), 1);

ac_time_axis = (1:FEATURE_WIN_LENGTH)/FEATURE_RATE;

for k = 1:NUM_FEATURE_FRAMES
    curr_feature_frame = feature_frame(k); % no windowing so far

    % go through and calculate autocorrelations for each slice of FEATURE_LENGTH features
    % The following are equivalent operations (from [3])
    % A = xcorr(a)  - traditional (unnormalised) autocorrelation
    % A = ifft(abs(fft(a, 2*length(a)).^2)) - i.e. zero padding a to 2x its length,
    % then IDFT the magnitude spectrum
    %ac = xcorr(curr_feature_frame, 'coeff');

    feature_spectrum = fft(curr_feature_frame, 2*FEATURE_WIN_LENGTH);
    %ac = ifft(abs(feature_spectrum).^0.5);
    ac = ifft(abs(feature_spectrum).^0.5);

    normalised_ac = ac;

    if 0 %mod(k, 10) == 0
        figure; subplot(2, 1, 1);
        plot(ac_time_axis(1:length(curr_feature_frame)), curr_feature_frame);
        subplot(2, 1, 2);
        plot(ac_time_axis, normalised_ac(1:FEATURE_WIN_LENGTH));
        title(sprintf('Autocorrelation for k = %d', k));
    end
    % normalise to have area?? 1
    tempo_peaks = find(ac == max(ac));
    tempo_confidences = max(ac);

    % find peaks until confidence is within a certain ratio of the highest
    % ac(tempo_peak) = 0;
    % tempo_peak = find(ac == max(ac)) -> find second peak, make sure it's far enough away
end



% vim: tabstop=4 expandtab


