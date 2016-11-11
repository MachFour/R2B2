function [evidence_from_tempo_harmonic] = collectEvidence(base_cluster, ...
    test_cluster)
    
    evidence_base = 0;
    evidence_test = 0; 
    
    for i=1:length(base_cluster) 
        for j=1:length(base_cluster{i})
            evidence_base = evidence_base + base_cluster{i}{j}(3);
        end
    end
    
    for i=1:length(test_cluster) 
        for j=1:length(test_cluster{i})
            % we only add evidence from harmonics that are present and
            % ignore the phase fundamental
            if test_cluster{i}{j}(5) > 1
                evidence_test = evidence_test + test_cluster{i}{j}(3);
            end
        end
    end

    % this metric is currently arbitrary so we will need to have a good
    % think about the best way to define this...
    evidence_from_tempo_harmonic = evidence_base + evidence_test;
end