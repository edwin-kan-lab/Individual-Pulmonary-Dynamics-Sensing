% Simple utility function to get peak frequency of signal in region of
% interest 
function freq = getPeakHarmonic(sig, fRange, fs) 
    [P1, f] = ncsFFT(sig, fs, false, fRange);
    idx = f>fRange(1) & f<fRange(2); 
    [~, loc] = max(P1(idx));
    freqs = f(idx);
    freq = freqs(loc); 
end 