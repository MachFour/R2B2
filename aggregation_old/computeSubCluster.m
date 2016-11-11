% feed in the original cluster and the harmonic index and this function
% will return only the points that appear in cluster that are centred
% around the harmonic
function [sub_cluster, num_points] = computeSubCluster(cluster, harmonic)
    num_points = 0;
    sub_cluster = {};
    sub_cluster_from_feature = {};
    
    for i=1:length(cluster)
       k = 1;
       for j=1:length(cluster{i})
            if cluster{i}{j}(5) == harmonic
                sub_cluster_from_feature{k} = cluster{i}{j};
                num_points = num_points + 1;
                k = k + 1;
            end
       end
       sub_cluster{i} = sub_cluster_from_feature;
       sub_cluster_from_feature = {};
    end
end