% Generalized function which gives fiducial points for respiration based on
% airflow 
%
% Aakash Kapoor (ak2247)

function [times, locs] = getFiducialPoints(sig, fs)
    sig = movmean(sig, fs/10);
    dx = normalize(diff(sig));
    dIn = dx; dOut = dx; 
    for i=1:length(dx)
        if sig(i) < 0
            dIn(i) = 0; 
            dOut(i) = -dx(i);
        else
            dIn(i) = dx(i); 
            dOut(i) = 0;
        end 
    end 

    pkdist = 0.9 * round(fs/getBreathingRate(sig, [0.1 0.4], fs));

    [~, inTimes] = findpeaks(dIn, 'MinPeakDistance', pkdist);
    [~, exTimes] = findpeaks(dOut, 'MinPeakDistance', pkdist);
        
    % Starting integration from inspire beginning only and making sure 
    % ends with inspire beginning as well. Find volume by integrating 
    % for each cycle

    % Two-pointer approach
    ptr1 = 1; 
    ptr2 = 1; 
    cptr = 1; 
    strArr = char(zeros(length(inTimes) + length(exTimes), 1)); 
    times = zeros(length(inTimes) + length(exTimes), 1);
    
    while 1 
        if inTimes(ptr1) < exTimes(ptr2)
            strArr(cptr) = 'I';
            times(cptr) = inTimes(ptr1); 
            if ptr1~=length(inTimes)
                ptr1 = ptr1 + 1; 
                cptr = cptr + 1; 
            else
                strArr(cptr) = 'I';
                cptr = cptr + 1; 
                for i=ptr2:length(exTimes)
                    strArr(cptr) = 'O';
                    cptr = cptr + 1; 
                end 
                break; 
            end 
        else 
            strArr(cptr) = 'O';
            times(cptr) = exTimes(ptr2); 
            if ptr2~=length(exTimes)
                ptr2 = ptr2 + 1; 
                cptr = cptr + 1; 
            else
                strArr(cptr) = 'O';
                cptr = cptr + 1; 
                for i=ptr1:length(inTimes)
                    strArr(cptr) = 'I';
                    cptr = cptr + 1; 
                end 
                break; 
            end 
        end 
    end 
        
    strArr = convertCharsToStrings(strArr); 
    locs = strfind(strArr, 'IOI'); 
end 