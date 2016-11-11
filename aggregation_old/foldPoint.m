function [dist] = foldPoint(x, x_ref, eps)
    dist = -1;
    allowedHarmonics = [1,2,3,4,5];
    
    for k=allowedHarmonics
        if x*k < x_ref + eps/(2*k) && x*k > x_ref - eps/(2*k)
            dist = (x*k);
            break;
        end
    end
end