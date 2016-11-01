%% This block runs Max's code and plots tempo vs. phase on a scatter plot
siarus = generate_data_set('test_songs/blue_rondo_ala_turk');

%%
clf
testFrame = 21;
plot_tempo_phase_estimates(siarus, testFrame);

%% This section is used to test the harmonic clustering algorithm
epsTempo = 0.1; % remember that these are relative distances!
allowed_harmonics = [1,2,3,4,5];
cluster = clusterByTempo(siarus, testFrame, allowed_harmonics, epsTempo); 

colours = {'green', 'red', 'blue', 'cyan','magenta','yellow'};

plotCluster(cluster, true); 

%% Now we deal with phase

% The following code in its current form assumes that the highest harmonic
% (once weighted) measures the true phase. This assumption is likely wrong
% and we'll need to come up with a way around it but it does give me a way
% to identify the strongest tempo hypothesis and test the method with which
% I do so.
allowed_phase_harmonics = [1,2,3,4];
eps_phase = 0.2;

[tempo, phase] = computeWinningTempo(cluster, allowed_harmonics, ...
    allowed_phase_harmonics, eps_phase);

%% Now lets try plotting tempo vs frame
frames = 1:length(siarus{1}); % get the number of frames
tempos = [];
for i=frames
    cluster = clusterByTempo(siarus, i, allowed_harmonics, epsTempo);
    [tempo, phase] = computeWinningTempo(cluster, allowed_harmonics, ...
        allowed_phase_harmonics, eps_phase);
    tempos = [tempos, tempo];
end

figure
scatter(frames, tempos, 'blue');