function [tempo, phase] = computeWinningTempo(cluster, allowed_harmonics, ...
    allowed_phase_harmonics, eps_phase)
    
    sub_clusters = {};
    num_points = [];
    for i=1:length(allowed_harmonics)
        [sub_clusters{i}, num_points(i)] = ...
            computeSubCluster(cluster, allowed_harmonics(i));
    end

    % compute the centre of mass for the phase at the highest harmonic
    highest_harmonic = 1;
    for i=length(allowed_harmonics):-1:1
        if num_points(i) > 0
            highest_harmonic = i;
            break;
        end
    end

    % we are basically choosing the phase with which we are going to work with
    % from here on
    % this is where the assumption is used!
    weighted_phases = [];
    total_phase_confidence = 0;
    [nah1, nah2, total_tempo_confidence] = getTempoData(cluster);
    for i=1:length(sub_clusters{highest_harmonic})
        for j=1:length(sub_clusters{highest_harmonic}{i})
            curr_point = sub_clusters{highest_harmonic}{i}{j};
            weighted_phases = [weighted_phases curr_point(3)*curr_point(4)];
            total_phase_confidence = total_phase_confidence + curr_point(4);
        end
    end

    base = tempoCentreOfMass(cluster, total_tempo_confidence)...
        /allowed_harmonics(highest_harmonic);
    phase_centre_of_mass = sum(weighted_phases)/total_phase_confidence;

    %hold on
    %    line([base - 0.1, base + 0.1], ...
    %        [phase_centre_of_mass phase_centre_of_mass], 'Color', 'green')
    %hold off

    % now that we have the base phase and the highest harmonic, we test the
    % remaining harmonics for multiple clusters to collect evidence...

    % for all harmonics greater than the one we want to test, test the tempo
    % phase structure of the one we're testing against the other harmonics to
    % collect evidence...
    
    %tempo_evidences = {};
    %colours = {'red', 'green', 'blue', 'magenta'};

    maxEvidence = 0;
    evidenceA = 0;
    tempo = 0;
    phase = 0;
    %present_phase_harmonics = 0;

    for i=(highest_harmonic - 1):-1:1
        % sub_cluster = clusterWithFixedOffset(cluster, base_phase,...
        %   tempo_offset, eps_phase);
        curr_cluster = clusterWithFixedOffset(sub_clusters{i},...
            phase_centre_of_mass, base, allowed_phase_harmonics, eps_phase);
        evidenceA = evidenceA + collectEvidence(sub_clusters{highest_harmonic},...
            curr_cluster);
        
        % just compute the evidence of a single harmonic without worrying
        % about evidence from other harmonics
        evidenceB = collectEvidence(curr_cluster, {});
        
        if evidenceB > maxEvidence
            tempo = base*allowed_harmonics(highest_harmonic)/allowed_harmonics(i);
            phase = phase_centre_of_mass;
            maxEvidence = evidenceB;
        end
        % this is just plotting for testing purposes...
        %hold on
        %for m=1:length(curr_cluster)
        %    for j=1:length(curr_cluster{m})
        %        for k=1:length(allowed_phase_harmonics)
        %            if curr_cluster{m}{j}(5) == allowed_phase_harmonics(k)
        %                scatter(curr_cluster{m}{j}(1), curr_cluster{m}{j}(3), ...
        %                     colours{k});
        %                if present_phase_harmonics < k
        %                    present_phase_harmonics = k;
        %                end
        %            end
        %        end
        %    end
        %end

        %for k=1:present_phase_harmonics
        %    line([base*(i+1) - 0.1, base*(i+1) + 0.1], ...
        %        [phase_centre_of_mass + (k-1)*base,...
        %        phase_centre_of_mass + (k-1)*base],...
        %        'Color', 'green')
        %end
        %hold off
    end
    
    if evidenceA > maxEvidence
        tempo = base;
        phase = phase_centre_of_mass;
    end

    % repeat with remaining harmonics, stopping when we reach the fundamental
    % harmonic...
end