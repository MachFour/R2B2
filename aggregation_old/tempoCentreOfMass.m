% note that this function needs the total_tempo_confidence to be
% pre-computed. this is easy as we usually need to call getTempoData
% beforehand in any case
function [centre_of_mass] = ...
    tempoCentreOfMass(cluster, total_tempo_confidence)

    weighted_tempos = []; % weighted by confidences

    for i=1:length(cluster)
        for j=1:length(cluster{i})
            curr_point = cluster{i}{j};
            % remember that curr_point(5) contains the harmonic number as
            % was computed during the clustering process
            weighted_tempos = [weighted_tempos,...
                curr_point(1)*curr_point(2)*curr_point(5)];
        end
    end

    centre_of_mass = sum(weighted_tempos)/total_tempo_confidence;
end