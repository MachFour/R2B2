% this function clusters by tempo using a harmonic metric such that points
% with close to 1/2, 1/3, 1/4 and 1/5 of the tempo are included in the same
% cluster. this function returns the winning cluster, which is the cluster
% that produces the highest cumulative confidence.
%
% note that there are a lot of efficiency fixes that can be made to this
% code. currently we are just brute force checking every point as the seed
% to a cluster when we could skip points with the same tempo as an
% example...
function [winning_cluster] = clusterByTempo(feature_estimates, ...
    frame, allowed_harmonics, eps_tempo)

    winning_cluster = {};
    top_score = 0;
    
    % strip away the frame layer for the purposes of clustering 
    data = {};
    for i=1:length(feature_estimates)
        data{i} = feature_estimates{i}{frame};
    end

    for i=1:length(data)
        for j=1:length(data{i})
            curr_cluster = gridQuery(data{i}{j}, data,...
                allowed_harmonics, eps_tempo);
            score = confidenceMetric(curr_cluster);
            if score > top_score
                winning_cluster = curr_cluster;
                top_score = score;
            end
        end
    end
end

% for a given seed point, this function computes a cluster given the 
% constraint of eps_tempo
function [cluster] = gridQuery(seed, data, allowed_harmonics, eps_tempo)
    cluster = {};
    cluster_from_feature = {};
    k = 1;
    
    for i=1:length(data) % iterate over all features
        for j=1:length(data{i}) % iterate over all estimates
            curr_point = data{i}{j};
            if  curr_point(1) <= seed(1)
                [fold, harmonic] = foldTempo(curr_point(1), seed(1), ...
                    allowed_harmonics, eps_tempo);
                if fold > 0 
                    curr_point(length(curr_point) + 1) = harmonic;
                    cluster_from_feature{k} = curr_point;
                    k = k + 1;
                end
            end
        end
        cluster{i} = cluster_from_feature;
        clusterPoints = {};
        k = 1;
    end  
end

% currently this just outputs the sum of both the phase and tempo
% confidences. nothing fancy but can be adjusted to produce different
% results
function [score] = confidenceMetric(cluster)
    score = 0;
    for i=1:length(cluster)
        phase_confidences = [];
        tempo_confidences = [];
        for k = 1:length(cluster{i})
            tempo_confidences = [tempo_confidences, cluster{i}{k}(2)];
            phase_confidences = [phase_confidences, cluster{i}{k}(4)];
        end
        score = score + sum(phase_confidences) + sum(tempo_confidences);
    end
end

function [fold, harmonic] = foldTempo(x, x_ref, allowed_harmonics, eps)
    fold = -1;
    for harmonic=allowed_harmonics
        if x*harmonic < x_ref + eps/(2*harmonic) && ...
                x*harmonic > x_ref - eps/(2*harmonic)
            fold = (x*harmonic);
            break;
        end
    end
end