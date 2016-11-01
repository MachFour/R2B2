% this algorithm clusters vertically
function [sub_cluster] = clusterWithFixedOffset(cluster, base_phase,...
    tempo_offset, allowed_phase_harmonics, eps_phase)
    
    sub_cluster = {};
    
    for i=1:length(cluster)
        sub_cluster{i} = {};
        m = 1;
        for j=1:length(cluster{i})
            curr_point = cluster{i}{j};
            for k=1:length(allowed_phase_harmonics)
                if curr_point(3) < (base_phase + (k-1)*tempo_offset) + eps_phase/2 ...
                        && curr_point(3) > (base_phase + (k-1)*tempo_offset)...
                        - eps_phase/2
                    curr_point(5) = allowed_phase_harmonics(k);
                    sub_cluster{i}{m} = curr_point;
                    m = m + 1;
                end
            end
        end
    end 
end