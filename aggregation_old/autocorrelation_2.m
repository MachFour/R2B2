% autocorrelation.m
% calculates the autocorrelation of a feature frame

function ac = autocorrelation(feature_frame)
    frame_length = length(feature_frame);

    % normalise to mean 0 and variance 1.
    mean = sum(feature_frame)/(frame_length);
    variance = sum(feature_frame.^2)/(frame_length-1);
    feature_frame = (feature_frame - mean)/sqrt(variance);
    
    ac = xcorr(feature_frame, 'none');
    % shift so that zero lag is at the first index
    ac = circshift(ac, frame_length);
    % discard redundant second half of samples and correct for shifting
    % bias (unbiased autocorrelation)
    ac = ac(1:frame_length)./(frame_length:-1:1)';
    
    % normalise so that first coefficient is 1, that way we can compare
    % between different autocorrelations
    ac = ac/ac(1);
    
    % The following are equivalent operations (from [3])
    % A = xcorr(a)  - traditional (unnormalised) autocorrelation
    % A = ifft(abs(fft(a, 2*length(a)).^2)) - i.e. zero padding a to 2x its length,
    % then IDFT the magnitude spectrum
    %feature_spectrum = fft(curr_feature_frame, 2*FEATURE_WIN_LENGTH);
    % square root suggested in "Streamlined Tempo Estimation..."
    %ac = ifft(abs(normalised_frame).^0.5);
end
