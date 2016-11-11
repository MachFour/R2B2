function [] = plot_tempo_phase_estimates(feature_estimates, frame)
    colours = {'green', 'red', 'blue', 'cyan','magenta','yellow'};
    NUM_FEATURE_CHANNELS = 5;
    
    glob_tempos = [];
    glob_phases = [];
    
    hold on
    for i = 1:NUM_FEATURE_CHANNELS
        tempos = [];
        phases = [];
        for k = 1:length(feature_estimates{i}{frame})
            tempos = [tempos, feature_estimates{i}{frame}{k}(1)];
            phases = [phases, feature_estimates{i}{frame}{k}(3)];
        end
        scatter(tempos, phases, colours{i});
        glob_tempos = [glob_tempos, tempos];
        glob_phases = [glob_phases, phases];
    end
    hold off
end