% Simple utility function to get breathing rate more robustly than peak
% harmonics
function freq = getBreathingRate(sig, fRange, fs) 
    [P1, f] = ncsFFT(sig, fs, false, fRange);
    [~, pks] = findpeaks(P1, 'MinPeakProminence', max(P1)/3);
    freqs =  f(pks); % Get top breathing harmonics
    idx = freqs>fRange(1) & freqs<fRange(2);
    freq = min(freqs(idx)); % Choose only 1 
end 