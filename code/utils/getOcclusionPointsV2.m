function [occlude_start, occlude_end] = getOcclusionPointsV2(sig, tFac)
    calib_win = 3; 
    hold_len = 3;
    tFac = 0.75; % Taking an initial threshold factor of 0.6
    % sig = sig.^2; % Experimental - Try square
    % Use the first few estimates as pre-occlusion levels 
    initialAvg = median(sig(1:calib_win)); 
    thresh = initialAvg * tFac; 
    % thresh = tFac; % because of disturbance, we do not take initial average 
    
    
    ctr = 0; 
    occlude_start = []; % Arrays to store occlusion start points 
    occlude_end = [];
    current_state = 1; % 0 = Occluded, 1 = Not Occluded
    
    for i=1:length(sig)
        if current_state == 0
            if sig(i) > thresh
                ctr = ctr + 1;
            end 
            if ctr >= hold_len
                % Occlusion has ended
                occlude_end = [occlude_end i];
                current_state = 1; 
                ctr = 0; 
            end 
        else 
            if sig(i) < thresh
                ctr = ctr + 1;
            end
            if ctr >= hold_len
                % Occlusion has started
                occlude_start = [occlude_start i];
                current_state = 0; 
                ctr = 0; 
            end 
        end 
    end 
end 