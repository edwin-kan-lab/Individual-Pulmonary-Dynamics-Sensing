function loc = getLocFromTime(locs, time_point, time_base, fs)    
    [~, loc] = min(abs(locs - (time_point - time_base) * fs));
end 