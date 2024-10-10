function template = getMultiTemplate(sig, points, loc, numpts, segLen)
    num_chs = size(sig, 2); 
    template = zeros(num_chs, segLen);
    for ch=1:num_chs    
        template(ch, :) = getTemplate(sig(:, ch), points, loc, numpts, segLen);
    end 
end 
