% Removes spikes from respiration signal, both at the smaller and larger
% scales, followed by filtering to remove DC trends. This is useful
% primarily for tidal volume extraction
% 
% Aakash Kapoor (ak2247)

function sig = respirationSpikeFilter(sig, fs, f_low, f_high)
    % Remove small spikes
    sig = medfilt1(sig, fs/2, 'omitnan');
    % Remove large spikes
    sig_diff = filloutliers(diff(sig), 'nearest', 'mean');
    sig = normalize(cumtrapz(sig_diff));
    % Filter 
    %sig = ncs_filt(sig, f_low, f_high, fs, 'bandpass');
    %sig = movmean(sig, fs/10); % OPTIONAL
end 