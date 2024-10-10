function template = getTemplate(sig, points, loc, numpts, segLen)
    pt = zeros(numpts, segLen);
    for i=1:numpts
        tmp = sig(points(loc-i):points(loc-i+1)-1);
        pt(i, :) = spline(1:length(tmp), tmp, linspace(1, length(tmp), segLen));
    end 

    % return average 
    template = mean(pt, 1);
end