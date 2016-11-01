function [centre_of_mass] = ...
    tempoCentreOfMass(cluster, epsTempo, epsPhase)
    
    [minTempo, maxTempo, minPhase, maxPhase, sumTConf, sumPConf] = ...
        clusterData(cluster);

    weighted_tempos = []; % weighted by confidences
    weighted_phases = [];
    for i=1:length(cluster)
        for j=1:length(cluster{i})
            currPoint = cluster{i}{j};
            weighted_tempos = [weighted_tempos,...
                foldPoint(currPoint(1), maxTempo, epsTempo)*currPoint(2)];
            weighted_phases = [weighted_phases,...
                foldPoint(currPoint(3), maxPhase, epsPhase)*currPoint(4)];
        end
    end

    tempo_centroid = sum(weighted_tempos)/sumTConf;
    phase_centroid = sum(weighted_phases)/sumPConf;
end