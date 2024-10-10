% This code converts airflow, F (L/s) from BIOPAC pneumotach transducer to
% volume, V (L). V is integration of F over time.
% May 03, 2019
% Pragya Sharma, ps847@cornell.edu
% 
% Updated: Mar 3, 2024
% Aakash Kapoor, ak2247@cornell.edu
% -------------------------------------------------------------------------

function [vol, airflow] = airflowToVolV1(airflow, fs, opts)
    % Airflow baseline correction
    meanDev = mean(airflow);
    airflow = airflow - meanDev;

    % Filtering for true zero-crossing detection
    if isfield(opts, 'filtType')
        if strcmp(opts.filtType, 'movmean')
            airflow = movmean(airflow, opts.winSize); 
        else
            % airflow = filterLpHp(airflow, fs, opts);
            airflow = ncs_filt(airflow, opts.fLow, opts.fHigh, fs, opts.filtType); 
        end 
    end 

    % Minimum time between two zero crossing points
    if ~isfield(opts, 'minTime')
        opts.minTime = 0.3; 
    end 

    % Now find volume by integrating for each cycle
    vol = zeros(length(airflow),1);
    t = (0:(length(airflow)-1))/fs;

    % z-cross does not give all the indices 
    if strcmp(opts.method, 'zcross')
        % Gives indices and pos/neg slope indication
        zcPosNeg = zeroCrossDet(airflow, fs, opts);

        % Starting integration from inspire beginning only. and making sure ends
        % with inspire beginning as well.
        if zcPosNeg(1,2) == -1
            zcPosNeg = zcPosNeg(2:end,:);
        end
        if zcPosNeg(end,2) == -1
            zcPosNeg = zcPosNeg(1:end-1,:);
        end
        
        for i = 1:2:(length(zcPosNeg))-2
            idxStart = zcPosNeg(i,1);
            idxEnd = zcPosNeg(i+2,1);
            
            vol(idxStart:idxEnd-1) = cumtrapz(t(idxStart:idxEnd-1), airflow(idxStart:idxEnd-1));
        end
    elseif strcmp(opts.method, 'peaks')
        % PREFERRED
        % Detect both positive and negative peaks
        % This fails when the end-inspiration peak is steeper than the
        % expiratory peak so is not as resilient a method 
        if ~isfield(opts, 'minPeakDist')
            opts.minPeakDist = int32(fs/4);
        end 
        % To mitigate the second positive peak from being chosen, smooth a
        % bit at the start
        af = movmean(airflow, floor(fs/20)); 
        [times, locs] = getFiducialPoints(af, 100); 
        
        for i=1:length(locs) 
            idxStart = times(locs(i));
            idxEnd = times(locs(i)+2); 
            
            vol(idxStart:idxEnd-1) = cumtrapz(t(idxStart:idxEnd-1), airflow(idxStart:idxEnd-1));
        end
    elseif strcmp(method, 'inpeaks')
        % Use peak detection only with inspiration peaks as they are the
        % only decently sized positive peaks 
        if ~isfield(opts, 'minPeakDist')
            opts.minPeakDist = int32(fs/4);
        end 
        % To mitigate the second positive peak from being chosen, smooth a
        % bit at the start
        af = movmean(airflow, floor(fs/20)); 
        dAir = normalize(diff(smooth(af)));
        
        % accounting for a 10% variability broadly 
        pkDist = 0.9 * round(fs/getPeakHarmonic(airflow, [0.1 0.4], fs)); 

        [~, inTimes] = findpeaks(dAir, 'MinPeakDistance', pkDist);
        
        for i=1:length(inTimes)-1
            vol(inTimes(i):inTimes(i+1)) = cumtrapz(airflow(inTimes(i):inTimes(i+1)));
        end 
    end 

    % At the end of this, sharp spikes in tidal volume will basically occur
    % only when spirometry is messed up. This will be used as an exclusion
    % criterion. 
end