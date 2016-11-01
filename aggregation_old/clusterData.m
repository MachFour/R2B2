% this function just computes useful values from the cluster data structure
function [min_tempo, max_tempo, total_tempo_confidence] = ...
    getTempoData(cluster)

    total_tempo_confidence = 0;
    max_tempo = 0;
    min_tempo = Inf;
    
    for i=1:length(cluster)
       for j=1:length(cluster{i})
           curr_point = cluster{i}{j};
           total_tempo_confidence = total_tempo_confidence + curr_point(2);
           if curr_point(1) < min_tempo
                min_tempo = curr_point(1);
           end
           if curr_point(1) > max_tempo
                max_tempo = curr_point(1);
           end
       end
    end
end