function [] = plotCluster(cluster, plot_lines)
    hold on
    for i=1:length(cluster)
       for j=1:length(cluster{i})
           currPoint = cluster{i}{j};
           scatter(currPoint(1), currPoint(3), 'black');
       end
    end

    [min_tempo, max_tempo, total_tempo_confidence] = getTempoData(cluster);

    base = tempoCentreOfMass(cluster, total_tempo_confidence);
    
    if plot_lines == true 
        for i=1:floor(max_tempo/min_tempo)
           line([base/i base/i], ...
               [0 1], 'Color', 'red');
        end
    end
    hold off
end