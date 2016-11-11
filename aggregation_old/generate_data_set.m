function [feature_estimates] = generate_data_set(wav_name)
    WAV_NAME = wav_name; % currently testing with 'someone'
    
    feature_tracker_2;
    
    feature_estimates = {};
    tempo_phase_estimates = {}; 

    tempos = {};
    phases = {};

    for i = 1:NUM_FEATURE_CHANNELS
        tempos  = feature_tempos{i};
        phases  = feature_phases{i};
        t_confs = feature_t_confs{i};
        p_confs = feature_p_confs{i};
        for j = 1:length(tempos)
            tempo_phase_estimates{j} = {};
            for k = 1:length(tempos{j})
                tuple = [tempos{j}(k), t_confs{j}(k),...
                    phases{j}(k), p_confs{j}(k)];
                tempo_phase_estimates{j}{k} = tuple;
            end
        end
        feature_estimates{i} = tempo_phase_estimates;
    end
end