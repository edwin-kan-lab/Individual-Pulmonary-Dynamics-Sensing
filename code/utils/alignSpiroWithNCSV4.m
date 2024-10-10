% Synchronize a given BIO trace with a given spirometer txt file which
% contains time stamps and annotations for events. Base time and offset are
% used to get total overall time axes as used by BIOPAC
% Spirometer used; Philips Respironics NM3 
% v4 uses v2 as base 
% Written by: Aakash Kapoor (ak2247)
% Updated: 15-Feb-2024

function [spiro_out, air_out, pao2_out, events] = ...
alignSpiroWithNCSV4(fName, baseTime, spiro_lead, tSteps, spiro_in, air_in, pao2_in)

    spiro_out = spiro_in;
    air_out = air_in;
    pao2_out = pao2_in;
    fileData = readlines(fName);
    
    opts = {};
    % movmean filter
    % opts.minPeakDist = 200; opts.method = 'peaks';
    % opts.filtType = 'movmean'; opts.winSize = 5;
    % ncs_filt
    opts.filtType = 'bandpass'; % removing DC will not impact Pk-Pk, only baseline
    opts.fHigh = 15; opts.fLow = 0.05;
    opts.method = 'peaks';

    for i=1:100
        if fileData(i).startsWith('Sample rate:')
            fs = str2num(extractAfter(fileData(i), ': '));
            % disp(['Spirometry is sampled at ', num2str(fs), 'Hz']);
        end 
        if fileData(i).startsWith('Creation date & time:')
            tmp = char(fileData(i));
            startDate = tmp(end-21:end-12);
            startTime = tmp(end-9:end);
            fileTime = datetime([startDate, ' ', startTime] , ...
                'InputFormat', 'dd.MM.yyyy HH:mm:ss.S');
            offset = seconds(fileTime - baseTime) - spiro_lead; 
        end 
        if fileData(i).startsWith('EVENTS:')
            eventIdx = i; 
        end 
        if fileData(i).startsWith('DATA:')
            dataIdx = i;
            % Extract events: from eventIdx to dataIdx - 2
            numEvents = dataIdx - eventIdx - 2;
            if numEvents ~=0
                for ind=1:numEvents 
                    tmp = strsplit(fileData(ind + eventIdx));
                    events(ind).time = tmp(3);
                    events(ind).label = [tmp{4:end}];
                end 
            else
                events = {};
            end 
            break;
        end 
    end

    clearvars fileData

    if ~exist('dataIdx', 'var')
        spiro_out = [];
        air_out = []; 
        events = {};
    else
        % Get all time stamps and data as (Time)(Flow)(PAO)(CO2)
        spiroData = table2array(readtable(fName, 'NumHeaderLines', dataIdx+1));
        spiroData(:, 1) = spiroData(:, 1) + offset; 
        
        [startA, endA, startB, endB, flag] = ...
                intersectArrays(spiroData(:, 1), tSteps, fs);

        if flag
            pao2_out(startB:endB) = spiroData(startA:endA, 3);
            % Obtain volume data by integrating airflow
            % Not the best algorithm. Gives worse results than cumtrapz for
            % good quality data but performs better in most cases
            [spiro_out(startB:endB), air_out(startB:endB)] = ...
                airflowToVolV1(spiroData(startA:endA, 2), fs, opts);
            
            % Locs are "fixed" by not being computed here as we can use the
            % same algorithm on airflow data at any time
        end
    end
end

function [startA, endA, startB, endB, flag] = ...
        intersectArrays(arrA, arrB, fs)

    startA = 1; endA = length(arrA);
    startB = 1; endB = length(arrB);
    
    if (arrA(startA) > arrB(endB) || arrB(startB) > arrA(endA))
        flag = false;
    else
        flag = true; 
        % Find common start position
        delStart = round((arrA(startA) - arrB(startB)) * fs);
        if delStart > 0
            startB = startB + delStart;
        else
            startA = startA - delStart;
        end
        % Find common end position
        delEnd = round((arrA(endA) - arrB(endB)) * fs);
        if delEnd > 0
            endA = endA - delEnd;
        else
            endB = endB + delEnd;
        end
    end 
end
