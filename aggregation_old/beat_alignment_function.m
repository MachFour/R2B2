% beat_alignment_function.m
% follows method in 'Context Dependent Beat tracking' to find the location 
% of beats in a frame of an onset detection function.
% beat spacing is given in samples.
function f = beat_alignment_function(frame, tempo_spacing)
    frame_length = length(frame);    
    % normalise frame data
    mean = sum(frame)/frame_length;
    variance = sum(frame.^2)/frame_length;
    
    % (we reverse frame to examine most recent onsets first)
    frame = (fliplr(frame) - mean)/sqrt(variance);
    SPIKE_WIDTH = 2;

    impulse_train = zeros(1, frame_length);
    % tempo is in lag_units (seconds)
    for w = 1:SPIKE_WIDTH
        impulse_train(w:tempo_spacing:end-tempo_spacing) = 1;       
    end
    % make impulses of decreasing height to weight more recent
    % onsets
   
    % exponential decay of weighting:
    % w(n) = a^(-a*n/N)
    % choose a = 2
    impulse_weighting = 2.^(-2/frame_length*(1:frame_length));
    impulse_train = impulse_train.*impulse_weighting;
    %figure; plot(impulse_train); hold on; plot(frame);
    %title(sprintf('Reversed feature frame and impulses for comb width %d', tempo_spacing));
    %xlabel('Samples');
    
    % only makes sense to consider shifts up to tempo_spacing, since
    % beyond that point the comb 'teeth' would fall off the end

    f = xcorr(frame, impulse_train, tempo_spacing);
    %figure; plot(f);
    %title('Beat alignment function');
    %xlabel('Comb shift (samples)');
    %ylabel('Strength');
    
    % return a vector of tempo_spacing length with zero lag as the first
    % sample
    f = circshift(f, [length(f), 0]);
    f = f(1:tempo_spacing);
end