% NIH PigStudy - Respiration
% This file loads the data files corresponding to a given intervention
%
% After loading, the file will load the spirometry data and synchronize it
% with the collected NCS trace. 
%
% The file also loads the notes, saved as .txt files, and indexes them 
% based on their time during the study for later interpretation
%
% Synchronization occurs using a constant time axis that times all of these
% monitors against midnight of the day the study started (0000)

% Note: Ch1 & Ch3 is Pig's Right side, Ch2 & Ch4 is Pig's left side
% NcsData1 is Ch1 and Ch2, NcsData2 is Ch3 and Ch4
% ncs. data format:
% [T1R1 Amp; T1R1 Phs; T2R1 Amp; T2R1 Phs; T1R2 Amp; T1R2 Phs; T2R2 Amp; T2R2 Phs;
%  T3R3 Amp; T3R3 Phs; T4R3 Amp; T4R3 Phs; T3R4 Amp; T3R4 Phs; T4R4 Amp; T4R4 Phs]

%% Declare file to load
addpath(genpath(pwd)); 
show.respiration = 0; % If running in automated mode

if ~exist('runAllInd', 'var')
    clearvars;
    clc;
    close all;
    clear vars; 
    show.respiration = 1; 
    dFileDate = '02-06-24'; 
    toSave = true;
    interventionType = 'lungocclusion'; 
    % No lung occlusion for 03-07-24
end

if ~exist('preload_ncs', 'var') 
    preload_ncs = false; 
end 

if ~exist('toSave', 'var') 
    toSave = false; 
end 

% Number of seconds spirometry is ahead of NCS
% eg. spiro reads 10s while NCS reads 0s
spiro_lead = 0; 

% Define relevant signal bounds
switch dFileDate
    case '02-05-24'
        spiro_lead = 7.3; 
        switch interventionType
            case 'lungocclusion'
                start_times = [46652];
                end_times = [47828];
                calib_st = 46657;
                calib_en = 46684;
                start_true = 46697;  
                deflate_true = 47558;
                cpap_true = 47767;
                stage = {'Occlusion'};
            case 'tidalvol'
                start_times = [44344, 44371, 44421, 44524, 44606, ...
                    44658, 44744, 44784]; 
                end_times = [44371, 44419, 44517, 44565, 44652, ...
                    44741, 44782, 44842]; 
                stage = [150, 220, 290, 360, 430, 500, 570, 640];
        end
    case '02-06-24'
        spiro_lead = 9.2; 
        switch interventionType
            case 'lungocclusion'
                start_times = [46152];
                end_times = [47405];
                calib_st = 46164; 
                calib_en = 46192; 
                start_true = 46205;  
                deflate_true = 46907;
                cpap_true = 47270;
                stage = {'Occlusion'};
            case 'tidalvol'
                start_times = [44140, 44201, 44260, 44329, 44387, ...
                    44514, 44572, 44623]; 
                end_times = [44201, 44260, 44329, 44387, 44501, ...
                    44572, 44623, 44658]; 
                stage = [150, 220, 290, 360, 430, 500, 570, 640];
        end
    case '02-07-24'
        spiro_lead = 1.7;
        switch interventionType
            case 'lungocclusion'
                start_times = [46607];
                end_times = [47864];
                calib_st = 46607;
                calib_en = 46637; 
                start_true = 46671;  
                deflate_true = 47569; 
                cpap_true = 47832;
                stage = {'Occlusion'};
            case 'tidalvol'
                start_times = [44860, 44916, 44954, 45048, 45110, ...
                    45167, 45212, 45275]; 
                end_times = [44916, 44954, 45048, 45105, 45167, ...
                    45205, 45273, 45330]; 
                stage = [150, 220, 290, 360, 430, 500, 570, 640];
        end
    case '02-08-24'
        spiro_lead = -0.4;
        switch interventionType
            case 'lungocclusion'
                start_times = [49341];
                end_times = [50166];
                calib_st = 49341;
                calib_en = 49371;
                start_true = 49383;  
                deflate_true = 49858; 
                cpap_true = 50053;
                stage = {'Occlusion'};
            case 'tidalvol'
                start_times = [47731, 47772, 47802, 47840, 47889, ...
                    47961, 47994, 48055]; 
                end_times = [47769, 47802, 47840, 47889, 47934, ...
                    47994, 48055, 48090]; 
                stage = [150, 220, 290, 360, 430, 500, 570, 640];
        end
    case '03-04-24'
        spiro_lead = 7.1; 
        switch interventionType
            case 'lungocclusion'
                start_times = [49898];
                end_times = [52215];
                start_true = 49958;  
                calib_st = 49928;
                calib_en = 49947;
                deflate_true = 51744; 
                cpap_true = 51863;
                stage = {'Occlusion'};
            case 'tidalvol'
                start_times = [45872, 45925, 46008, 46096, 46172, ...
                    46240, 46312, 46382]; 
                end_times = [45916, 45988, 46080, 46160, 46230, ... 
                    46300, 46372, 46426]; 
                stage = [150, 220, 290, 360, 430, 500, 570, 640];
        end
    case '03-05-24'
        spiro_lead = -8.8;
        switch interventionType
            case 'lungocclusion'
                start_times = [44167];
                end_times = [45620];
                start_true = 44197;  
                calib_st = 44174; 
                calib_en = 44193; 
                deflate_true = 45224; 
                cpap_true = 45515;
                stage = {'Occlusion'};
            case 'tidalvol'
                start_times = [41667, 41700, 41735, 41775, 41819, ...
                    41859, 41923, 42012]; 
                end_times = [41695 41733, 41762, 41815, 41855, 41919, ...
                    42008, 42036]; 
                stage = [150, 220, 290, 360, 430, 500, 570, 640];
        end
    case '03-06-24'
        spiro_lead = 5; 
        switch interventionType
            case 'lungocclusion'
                start_times = [43001];
                end_times = [44867];
                start_true = 43132;  
                deflate_true = 44383; 
                calib_st = 43060; 
                calib_en = 43090;
                cpap_true = 44726;
                stage = {'Occlusion'};
            case 'tidalvol'
                start_times = [40127, 40190, 40231, 40258, 40299, ...
                    40335, 40374, 40413]; 
                end_times = [40188, 40227, 40255, 40295, 40331, ...
                    40370, 40412, 40486]; 
                stage = [150, 220, 290, 360, 430, 500, 570, 640];
        end
    case '03-07-24'
        spiro_lead = 6.7; 
        switch interventionType
            case 'tidalvol'
                start_times = [46599, 46617, 46656, 46722, 46783, ...
                    46836, 46888, 46953]; 
                end_times = [46609, 46646, 46716, 46780, 46832, ...
                    46885, 46937, 46985]; 
                stage = [150, 220, 290, 360, 430, 500, 570, 640];
        end
end

% Standardized name for file saving
interventionType_Orig = interventionType; 

% Fix for some files
switch dFileDate
    case '02-06-24'
        if strcmp(interventionType, 'tidalvol')
            interventionType = 'tidalvolume';
        end 
    case '02-07-24'
        if strcmp(interventionType, 'lungocclusion')
            interventionType = 'lungocclusion2';
        end
    case '02-08-24'
        if strcmp(interventionType, 'tidalvol')
            interventionType = 'tidalvol3';
        end
    case '03-04-24'
        if strcmp(interventionType, 'lungocclusion')
            interventionType = 'lungocclusion3';
        end 
    case '03-05-24'
        if strcmp(interventionType, 'tidalvol')
            interventionType = 'tidalvol2';
        end
    case '03-07-24'
        if strcmp(interventionType, 'lungocclusion')
            return;
        end 
end 

disp('File Specified');

%% Load files
dBasePath = 'E:\OneDrive - Cornell University\NIH Pig Study\'; 

% Code Speed-Up by preloading Saved Data file  
if preload_ncs
    load(fullfile(dBasePath, 'NIH Saved Data', [dFileDate ' ' ...
        interventionType_Orig 'V1.mat']));
    disp('Pre-loaded everything'); 
else
    d = dir(dBasePath);
    found = false;
    for i = 1:size(d, 1)
        if contains(d(i).name, dFileDate)
            dPath = fullfile(dBasePath, d(i).name, 'Data Files');
            d = dir(dPath);
            found = true;
        end
    end
    if ~found; error('Date not found'); end

    found = false;
    for i = 1:size(d, 1)
        if contains(d(i).name, 'NCS')
            if contains(d(i).name, '.mat')
                if contains(d(i).name, interventionType)
                    NcsFile = d(i).name;
                    found = true;
                end
            end
        end
    end

    if ~found; error('Intervention type not found'); end
    BioFile = strrep(NcsFile, 'NCS', 'BIO');
    RadioFile = strrep(NcsFile, 'NCS', 'Radio');

    disp('Loading');
    if isfile(fullfile(dPath, NcsFile))
        load(fullfile(dPath, NcsFile), 'NcsData1', 'NcsData2');
    else
        disp('Converting NCS');
        saveTDMStoMAT_NCSpig(dPath, strrep(NcsFile, '.mat', ''));
        load(fullfile(dPath, NcsFile), 'NcsData1', 'NcsData2');
    end

    % Gets rid of phase discontinuity at the beginning
    for i = 2:2:size(NcsData1, 2)
        NcsData1(1:20000, i) = median(NcsData1(20000:40000, i));
        NcsData2(1:20000, i) = median(NcsData2(20000:40000, i));
    end

    if isfile(fullfile(dPath, BioFile))
        load(fullfile(dPath, BioFile), 'BioData');
    else
        disp('Converting BIO');
        saveTDMStoMAT_BIOpig(dPath, strrep(BioFile, '.mat', ''));
        load(fullfile(dPath, BioFile), 'BioData');
    end

    if isfile(fullfile(dPath, RadioFile))
        load(fullfile(dPath, RadioFile), 'RadioData');
    else
        disp('Converting RADIO');
        saveTDMStoMAT_Radiopig(dPath, strrep(RadioFile, '.mat', ''));
        load(fullfile(dPath, RadioFile), 'RadioData');
    end

    radio.data = RadioData;
    radio.diff = diff(radio.data);

    load(fullfile(dBasePath, ['Pig Study ' dFileDate], ...
        [dFileDate '_powerlab_extracted_data_V3.mat']));

    % Extract date and time from filename
    date_.full = NcsFile(strfind(NcsFile, 'Saving')+6:end-4);
    date_.month = date_.full(1:2);
    date_.date = date_.full(3:4);
    date_.year = num2str(year(powerlab.record_info(1).data_start_dt));
    date_.hour = date_.full(6:7);
    date_.min = date_.full(8:9);
    date_.sec = date_.full(10:11);
    date_.filetime = labviewDateToDateTime(date_.full, date_.year);
    date_.baseline = datetime([date_.year '-' date_.month '-' ...
        date_.date ' 00:00:00'], 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
    date_.seconds_offset = seconds(date_.filetime - date_.baseline);

    disp('Loaded');
    clearvars i found d RadioData
 

    %% Load and sync spirometry 
    % There is a time lag so we need to fix this by making use of 
    % alignment in breath hold segments. This has been hard-coded 

    disp('Loading Spirometry');
    d = dir(dBasePath);
    found = false;
    for i = 1:size(d, 1)
        if contains(d(i).name, dFileDate)
            dPath = fullfile(dBasePath, d(i).name, 'Spirometry');
            d = dir(dPath);
            found = true;
            break;
        end
    end
    if ~found; error('Spirometry data not found'); end

    fSpiro = 100; 
    spiro.tSpiro = (1/fSpiro:1/fSpiro:length(NcsData1)/10000) + ...
                        date_.seconds_offset;
    spiro_data = zeros(length(spiro.tSpiro), 1);
    airflow_data = zeros(length(spiro.tSpiro), 1);
    pao2_data = zeros(length(spiro.tSpiro), 1);

    % Choose a baseline NCS to calibrate against

    for i=1:size(d, 1)
        if contains(d(i).name, '.txt')
            fName = fullfile(dPath, d(i).name);
            [spiro_data, airflow_data, pao2_data, events] = ... 
                alignSpiroWithNCSV4(fName, date_.baseline, ...
                    spiro_lead, spiro.tSpiro, spiro_data, ...
                    airflow_data, pao2_data);
        end 
    end

    spiro.volume = spiro_data; 
    spiro.airflow = airflow_data; 
    spiro.pao2 = pao2_data;
    spiro.fSpiro = fSpiro; 
    
    [times, locs] = getFiducialPoints(movmean(airflow_data, ...
        floor(fSpiro/20)), fSpiro);
    spiro.locs = times(locs); 

    % locs will automatically exclude BH locations so could be used for NCS
    % segmentation as well 

    if show.respiration
        figure('Position', [300 50 1200 900]); 
        
        ax1 = subplot(311); 
        plot(spiro.tSpiro, spiro.airflow); 
        title('Airflow (L/s)');
        
        ax2 = subplot(312); 
        plot(spiro.tSpiro, spiro.volume, ...
            '-s', 'MarkerIndices', spiro.locs, 'Color', '#4DBEEE', ...
            'MarkerFaceColor', '#77AC30', 'MarkerSize', 14, ...
            'LineWidth', 1.5);
        title('Tidal Volume (L)');
        
        ax3 = subplot(313); 
        plot(spiro.tSpiro, spiro.pao2);
        
        title('PaO2 (cm H20)');
        linkaxes([ax1, ax2, ax3], 'x'); 
    end 

    clearvars fSpiro found d dPath fName airflow_data pao2_data spiro_data
    disp('Spirometry Loaded');

    %% Load and index notes
    disp('Loading Notes');
    d = dir(dBasePath); 
    for i = 1:size(d, 1)
        if contains(d(i).name, dFileDate)
            dNotesPath = fullfile(dBasePath, d(i).name, 'Notes');
            d = dir(dNotesPath);
            break;
        end
    end

    inc = 1; notes = [];
    for i = 1:size(d, 1)
        if contains(d(i).name, '.txt')
            f = fileread(fullfile(dNotesPath, d(i).name));
            notes(inc).dateStr = d(i).name(end-14:end-4);
            notes(inc).datetime = ...
                labviewDateToDateTime(notes(inc).dateStr, date_.year);
            notes(inc).seconds = seconds(notes(inc).datetime - ...
                date_.baseline);
            notes(inc).text = f;
            inc = inc+1;
        end
    end

    clearvars i inc d f
    disp('Notes Loaded');

    %% Declare synchronous time axes
    disp('Syncing NCS and BIO');
    
    % Sync NCS and BIO
    p = 7.3306e-05; tOffset = 0.045;
    fNcsOrig = 10000; fBio = 1000/(1-p);
    tNcs = (1/fNcsOrig:1/fNcsOrig:length(NcsData1)/fNcsOrig) + ...
        date_.seconds_offset;
    bio.tBio = (1/fBio:1/fBio:length(BioData)/fBio) + tOffset + ...
        date_.seconds_offset;
    radio.tRadio = linspace(min(tNcs), max(tNcs), size(radio.data, 1)); 
    radio.fRadio = 1/(radio.tRadio(2)-radio.tRadio(1));
    
    % Store BioData
    bio.data = BioData; bio.fBio = fBio;
    clearvars p tOffset BioData fBio;
    
    %% NCS downsampling & filtering
    disp('Processing NCS');

    ncs.dsRate = 10; 
    ncs.fNcsDS = fNcsOrig/ncs.dsRate;
    
    ncs.ncsDS = zeros(length(NcsData1(:,1))/ncs.dsRate, ...
        size(NcsData1,2)*2);
    ncs.tNcsDS = tNcs(1:ncs.dsRate:end);
    
    for ch = 1:size(NcsData1, 2)
        % 900 MHz
        ncs.ncsDS(:, ch) = NcsData1(1:ncs.dsRate:end, ch);
        % 2 GHz
        ncs.ncsDS(:, ch+size(NcsData1,2)) = NcsData2(1:ncs.dsRate:end, ch);
    end
    ncs.chs = size(ncs.ncsDS,2);

    % NCS Filtering 
    ncs.ncsSpikeFilt = zeros(length(ncs.ncsDS)-1, 16);
    ncs.ncsRespFilt = zeros(length(ncs.ncsDS)-1, 16);
    for ch = 1:ncs.chs
        % Unwrap Phase
        if mod(ch, 2) == 0
            ncs.ncsDS(:, ch) = unwrap(ncs.ncsDS(:, ch) * pi/180);
        end 
        
        ncs.ncsSpikeFilt(:, ch) = ...
            respirationSpikeFilter(ncs.ncsDS(:, ch), ncs.fNcsDS, 0.05, 15); 
        ncs.ncsRespFilt(:, ch) = ncs_filt(ncs.ncsSpikeFilt(:, ch), ...
            0.05, 15, ncs.fNcsDS, 'bandpass');
    end
    
    % NCS CVI 
    if ~exist('opts', 'var')
        opts = {};
    end
    opts.do_cvi = false;
    
    if opts.do_cvi
        ncs.ncsCVI = zeros(length(ncs.ncsDS)-1, 16);
        cvi_opts = {}; 
        cvi_opts.fRange = [0.18 0.36]; 
        cvi_opts.optFcn = 'SNR';
        cvi_opts.magSteps = 4; 
        cvi_opts.phaseSteps = 150;
    
        optFcn = 'lin'; 
        for ch=1:2:16
            if ismember(ch, [3, 5, 11, 13])
                ncs.ncsCVI(:, ch) = ncs.ncsRespFilt(:, ch);
                ncs.ncsCVI(:, ch+1) = ncs.ncsRespFilt(:, ch+1);
            else
            [ncs.ncsCVI(:, ch), ncs.ncsCVI(:, ch+1)] = ...
                cvi(ncs.ncsRespFilt(:, ch), ncs.ncsRespFilt(:, ch+1), ...
                ncs.fNcsDS, cvi_opts);
            end 
        end 
    end 
    
    clearvars BP* HP* LP* NcsData* tNcs fNcsOrig foo
    disp('Completed Pre-Processing');
end 

%% Plot raw NCS and reference signals with notes
if show.respiration
    figure('Position', [300 50 1200 900]); clearvars ax;

    ax(1)=subplot('Position', [0.05 0.82 0.9 0.16]);
    yyaxis left;
    plot(ncs.tNcsDS, movmean(ncs.ncsDS(:, 1), 10), 'k', 'DisplayName', 'Amp', 'LineStyle', '-'); hold on; 
    yyaxis right;
    plot(ncs.tNcsDS, movmean(ncs.ncsDS(:, 2), 10), 'b', 'DisplayName', 'Phs', 'LineStyle', '-');
    title(['Ch1 (Right) 900 MHz RF ' dFileDate ' ' interventionType_Orig]); ax(1).XGrid = 'on';
    legend('Location', 'southeast', 'Orientation', 'horizontal'); xticklabels(ax(1), {});

    ax(2)=subplot('Position', [0.05 0.64 0.9 0.16]);
    yyaxis left;
    plot(ncs.tNcsDS, movmean(ncs.ncsDS(:, 7), 10), 'k', 'DisplayName', 'Amp', 'LineStyle', '-'); hold on; 
    yyaxis right;
    plot(ncs.tNcsDS, ncs.ncsDS(:, 8), 'b', 'DisplayName', 'Phs', 'LineStyle', '-');
    title('Ch2 (Left) 900 MHz RF'); ax(2).XGrid = 'on';
    legend('Location', 'southeast', 'Orientation', 'horizontal'); xticklabels(ax(2), {});
    
    ax(3)=subplot('Position', [0.05 0.46 0.9 0.16]); 
    yyaxis left;
    plot(ncs.tNcsDS, ncs.ncsDS(:, 9), 'k', 'DisplayName', 'Amp', 'LineStyle', '-'); hold on; 
    yyaxis right;
    plot(ncs.tNcsDS, ncs.ncsDS(:, 10), 'b', 'DisplayName', 'Phs', 'LineStyle', '-');
    title('Ch3 (Right) 2 GHz RF'); ax(3).XGrid = 'on';
    legend('Location', 'southeast', 'Orientation', 'horizontal'); xticklabels(ax(3), {});
    
    ax(4)=subplot('Position', [0.05 0.22 0.9 0.16]);
    yyaxis left;
    plot(ncs.tNcsDS, ncs.ncsDS(:, 15), 'k', 'DisplayName', 'Amp', 'LineStyle', '-'); hold on; 
    yyaxis right;
    plot(ncs.tNcsDS, ncs.ncsDS(:, 16), 'b', 'DisplayName', 'Phs', 'LineStyle', '-');
    title('Ch4 (Left) 2 GHz RF'); ax(4).XGrid = 'on';
    legend('Location', 'southeast', 'Orientation', 'horizontal'); xticklabels(ax(4), {});

    ax(5)=subplot('Position', [0.05 0.04 0.9 0.16]);
    yyaxis left; ax(5).YColor = 'k';
    
    plot(spiro.tSpiro, normalize(spiro.volume), 'k', 'DisplayName', 'Tidal Volume (L)'); hold on; 
    plot(spiro.tSpiro, normalize(spiro.pao2), 'g', 'DisplayName', 'PaO2 (cm H20)');
 
    yyaxis right; ax(5).YColor = 'b';

    title('Reference Features'); ax(5).XGrid = 'on';
    legend('Location', 'southeast', 'Orientation', 'horizontal');
    ax(5).XAxis.Exponent = 0; xtickformat(ax(5), '%.2f');
    for i = 1:size(notes, 2)
        if notes(i).seconds > min(ncs.tNcsDS) && notes(i).seconds < max(ncs.tNcsDS)
            text(ax(5), notes(i).seconds, 0.3, notes(i).text, 'Rotation', 60,...
                'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', 'r', 'FontWeight', 'bold');
        end
    end

    ax(6)=subplot('Position', [0.05 0.001 0.9 0.02]); 
    ax(6).XGrid = 'off'; ax(6).YGrid = 'off'; 
    plot(bio.tBio, zeros(size(bio.tBio))); xticklabels(ax(6), {}); yticklabels(ax(6), {});
    
    linkaxes(ax, 'x'); 
    xlim([min(ncs.tNcsDS)+5 max(ncs.tNcsDS)]);
    clearvars ax i
end

%% Plot Filtered NCS and reference signals with notes
if show.respiration
    figure('Position', [300 50 1200 900]); clearvars ax;

    ax(1)=subplot('Position', [0.05 0.82 0.9 0.16]);
    yyaxis left;
    plot(ncs.tNcsDS(1:end-1), ncs.ncsRespFilt(:, 1), 'k', 'DisplayName', 'Amp', 'LineStyle', '-'); hold on; 
    yyaxis right;
    plot(ncs.tNcsDS(1:end-1), ncs.ncsRespFilt(:, 2), 'b', 'DisplayName', 'Phs', 'LineStyle', '-');
    title(['Ch1 (Right) 900 MHz RF Filtered ' dFileDate ' ' interventionType_Orig]); ax(1).XGrid = 'on';
    legend('Location', 'southeast', 'Orientation', 'horizontal'); xticklabels(ax(1), {});

    ax(2)=subplot('Position', [0.05 0.64 0.9 0.16]);
    yyaxis left;
    plot(ncs.tNcsDS(1:end-1), ncs.ncsRespFilt(:, 7), 'k', 'DisplayName', 'Amp', 'LineStyle', '-'); hold on; 
    yyaxis right;
    plot(ncs.tNcsDS(1:end-1), ncs.ncsRespFilt(:, 8), 'b', 'DisplayName', 'Phs', 'LineStyle', '-');
    title('Ch2 (Left) 900 MHz RF Filtered '); ax(2).XGrid = 'on';
    legend('Location', 'southeast', 'Orientation', 'horizontal'); xticklabels(ax(2), {});
    
    ax(3)=subplot('Position', [0.05 0.46 0.9 0.16]); 
    yyaxis left;
    plot(ncs.tNcsDS(1:end-1), ncs.ncsRespFilt(:, 9), 'k', 'DisplayName', 'Amp', 'LineStyle', '-'); hold on; 
    yyaxis right;
    plot(ncs.tNcsDS(1:end-1), ncs.ncsRespFilt(:, 10), 'b', 'DisplayName', 'Phs', 'LineStyle', '-');
    title('Ch3 (Right) 2 GHz RF Filtered '); ax(3).XGrid = 'on';
    legend('Location', 'southeast', 'Orientation', 'horizontal'); xticklabels(ax(3), {});
    
    ax(4)=subplot('Position', [0.05 0.22 0.9 0.16]);
    yyaxis left;
    plot(ncs.tNcsDS(1:end-1), ncs.ncsRespFilt(:, 15), 'k', 'DisplayName', 'Amp', 'LineStyle', '-'); hold on; 
    yyaxis right;
    plot(ncs.tNcsDS(1:end-1), ncs.ncsRespFilt(:, 16), 'b', 'DisplayName', 'Phs', 'LineStyle', '-');
    title('Ch4 (Left) 2 GHz RF Filtered '); ax(4).XGrid = 'on';
    legend('Location', 'southeast', 'Orientation', 'horizontal'); xticklabels(ax(4), {});

    ax(5)=subplot('Position', [0.05 0.04 0.9 0.16]);
    yyaxis left; ax(5).YColor = 'k';
    
    plot(spiro.tSpiro, normalize(spiro.volume), 'k', 'DisplayName', 'Tidal Volume (L)'); hold on; 
    plot(spiro.tSpiro, normalize(spiro.pao2), 'g', 'DisplayName', 'PaO2 (cm H20)');
 
    yyaxis right; ax(5).YColor = 'b';

    title('Reference Features'); ax(5).XGrid = 'on';
    legend('Location', 'southeast', 'Orientation', 'horizontal');
    ax(5).XAxis.Exponent = 0; xtickformat(ax(5), '%.2f');
    for i = 1:size(notes, 2)
        if notes(i).seconds > min(ncs.tNcsDS) && notes(i).seconds < max(ncs.tNcsDS)
            text(ax(5), notes(i).seconds, 0.3, notes(i).text, 'Rotation', 60,...
                'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', 'r', 'FontWeight', 'bold');
        end
    end

    ax(6)=subplot('Position', [0.05 0.001 0.9 0.02]); 
    ax(6).XGrid = 'off'; ax(6).YGrid = 'off'; 
    plot(bio.tBio, zeros(size(bio.tBio))); xticklabels(ax(6), {}); yticklabels(ax(6), {});
    
    linkaxes(ax, 'x'); 
    xlim([min(ncs.tNcsDS)+5 max(ncs.tNcsDS)]);
    clearvars ax i
end

%% Plot CVI if exists
if show.respiration & opts.do_cvi
    figure('Position', [300 50 1200 900]); clearvars ax;

    ax(1)=subplot('Position', [0.05 0.82 0.9 0.16]);
    yyaxis left;
    plot(ncs.tNcsDS(1:end-1), ncs.ncsCVI(:, 1), 'k', 'DisplayName', 'Amp', 'LineStyle', '-'); hold on; 
    yyaxis right;
    plot(ncs.tNcsDS(1:end-1), ncs.ncsCVI(:, 2), 'b', 'DisplayName', 'Phs', 'LineStyle', '-');
    title(['Ch1 (Right) 900 MHz RF Filtered ' dFileDate ' ' interventionType_Orig]); ax(1).XGrid = 'on';
    legend('Location', 'southeast', 'Orientation', 'horizontal'); xticklabels(ax(1), {});

    ax(2)=subplot('Position', [0.05 0.64 0.9 0.16]);
    yyaxis left;
    plot(ncs.tNcsDS(1:end-1), ncs.ncsCVI(:, 7), 'k', 'DisplayName', 'Amp', 'LineStyle', '-'); hold on; 
    yyaxis right;
    plot(ncs.tNcsDS(1:end-1), ncs.ncsCVI(:, 8), 'b', 'DisplayName', 'Phs', 'LineStyle', '-');
    title('Ch2 (Left) 900 MHz RF Filtered '); ax(2).XGrid = 'on';
    legend('Location', 'southeast', 'Orientation', 'horizontal'); xticklabels(ax(2), {});
    
    ax(3)=subplot('Position', [0.05 0.46 0.9 0.16]); 
    yyaxis left;
    plot(ncs.tNcsDS(1:end-1), ncs.ncsCVI(:, 9), 'k', 'DisplayName', 'Amp', 'LineStyle', '-'); hold on; 
    yyaxis right;
    plot(ncs.tNcsDS(1:end-1), ncs.ncsCVI(:, 10), 'b', 'DisplayName', 'Phs', 'LineStyle', '-');
    title('Ch3 (Right) 2 GHz RF Filtered '); ax(3).XGrid = 'on';
    legend('Location', 'southeast', 'Orientation', 'horizontal'); xticklabels(ax(3), {});
    
    ax(4)=subplot('Position', [0.05 0.22 0.9 0.16]);
    yyaxis left;
    plot(ncs.tNcsDS(1:end-1), ncs.ncsCVI(:, 15), 'k', 'DisplayName', 'Amp', 'LineStyle', '-'); hold on; 
    yyaxis right;
    plot(ncs.tNcsDS(1:end-1), ncs.ncsCVI(:, 16), 'b', 'DisplayName', 'Phs', 'LineStyle', '-');
    title('Ch4 (Left) 2 GHz RF Filtered '); ax(4).XGrid = 'on';
    legend('Location', 'southeast', 'Orientation', 'horizontal'); xticklabels(ax(4), {});

    ax(5)=subplot('Position', [0.05 0.04 0.9 0.16]);
    yyaxis left; ax(5).YColor = 'k';
    
    plot(spiro.tSpiro, normalize(spiro.volume), 'k', 'DisplayName', 'Tidal Volume (L)'); hold on; 
    plot(spiro.tSpiro, normalize(spiro.pao2), 'g', 'DisplayName', 'PaO2 (cm H20)');
 
    yyaxis right; ax(5).YColor = 'b';

    title('Reference Features'); ax(5).XGrid = 'on';
    legend('Location', 'southeast', 'Orientation', 'horizontal');
    ax(5).XAxis.Exponent = 0; xtickformat(ax(5), '%.2f');
    for i = 1:size(notes, 2)
        if notes(i).seconds > min(ncs.tNcsDS) && notes(i).seconds < max(ncs.tNcsDS)
            text(ax(5), notes(i).seconds, 0.3, notes(i).text, 'Rotation', 60,...
                'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', 'r', 'FontWeight', 'bold');
        end
    end

    ax(6)=subplot('Position', [0.05 0.001 0.9 0.02]); 
    ax(6).XGrid = 'off'; ax(6).YGrid = 'off'; 
    plot(bio.tBio, zeros(size(bio.tBio))); xticklabels(ax(6), {}); yticklabels(ax(6), {});
    
    linkaxes(ax, 'x'); 
    xlim([min(ncs.tNcsDS)+5 max(ncs.tNcsDS)]);
    clearvars ax i
end

%% Volumetric Feature Extraction
% Extract features per-stage to prevent artifacts from excluded breaths 
opts.win = 3;
opts.complex = false; 
opts.seg = 1000; 
opts.num_comps = 1; 
opts.do_cvi = false; 

if opts.do_cvi
    opts.sel_channels = [1, 7, 9, 15];
else 
    opts.sel_channels = [1, 2, 7, 8, 9, 10, 15, 16];
end

ncs.locs = (spiro.locs - 1) * 10 + 1; 

% Each element per stage is num_seg x num_channels 
ncs.tidalVol = cell(length(stage), 1); 
ncs.airflow = cell(length(stage), 1); 
ncs.baseline = cell(length(stage), 1); 

% Each element per stage is num_vars x 1 
spiro.tidalVol = cell(length(stage), 1); 
spiro.airflow = cell(length(stage), 1);
spiro.baseline = cell(length(stage), 1); 

% Obtain Scoring Metrics for signal 
snr_scores = zeros(length(stage), length(opts.sel_channels)); 
lin_scores = zeros(length(stage), length(opts.sel_channels)); 

% Choose signal to process 
spiro_sig = spiro.volume; 
if opts.do_cvi
    ncs_sig = ncs.ncsCVI(:, opts.sel_channels);
else
    ncs_sig = ncs.ncsRespFilt(:, opts.sel_channels);
end

% Get PCA scores 
[coeff,score,latent,tsquared,explained,mu] = ...
    pca(ncs_sig);

for i=1:length(stage)
    % Get starting and ending values for each stage 
    start_loc = getLocFromTime(spiro.locs, start_times(i), ...
        date_.seconds_offset, 100);
    end_loc = getLocFromTime(spiro.locs, end_times(i), ...
        date_.seconds_offset, 100);
    
    opts.spiro_locs = spiro.locs(start_loc:end_loc);
    opts.ncs_locs = ncs.locs(start_loc:end_loc);
    
    [ncs_TVs, ncs_PCA, ncs_Base, spiro_TVs, snrs, lins] = ...
        getTidalVolumesV3(ncs_sig, spiro_sig, opts);
    
    snr_scores(i, :) = snrs; 
    lin_scores(i, :) = lins; 
    spiro.tidalVol{i} = spiro_TVs;
    if opts.complex
        opts.jump = 2;
    else 
        opts.jump = 1;
    end 
    ncs.tidalVol{i} = ncs_TVs(:, [1:opts.jump:8]);
    ncs.pcaTV{i} = ncs_PCA * ...
        explained(1:opts.num_comps)/sum(explained(1:opts.num_comps)); 
    ncs.baseline{i} = ncs_Base;
end 

calib_method = 'single'; % default
overall_start_loc = getLocFromTime(spiro.locs, start_times(1), date_.seconds_offset, 100);

clearvars ncs_TVs ncs_Base spiro_TVs spiro_Base cvi_TVs cvi_Base snrs lins

%% Run All Analysis Combinations 
% calibs = {'single', 'cubic', 'twopoint'};
calibs = {'single'}
combinations = {'all', 'topk'};
scorings = {'snr', 'lin', 'unif'}; 
occhs = {'lr', 'lr2', 'lrC'};

if strcmp(interventionType_Orig, 'lungocclusion')
    for sc=1:length(scorings)
        chs = occhs{sc};
        generateOcclusionResults(opts, overall_start_loc, spiro, ...
            end_times, date_, start_true, deflate_true, cpap_true, ncs, ...
            calib_st, calib_en, dFileDate, chs, toSave);
    end 
elseif strcmp(interventionType_Orig, 'tidalvol')
    for cind = 1:length(calibs)
        calib_method = calibs{cind};
        for sc = 1:length(scorings)
            scoring = scorings{sc};
            for cx = 1:length(combinations)
                comb = combinations{cx};     
                generateTidalVolumeResults(ncs, spiro, stage, comb, ...
                    calib_method, scoring, snr_scores, lin_scores, ...
                    dFileDate, toSave);
            end
        end 
    end 
end 
 

clearvars scorings combinations calibs scoring calib_method comb cind sc cx

%% Save necessary data 
version = '1';
if toSave
    clearvars NcsData1 NcsData2 ch i 
    disp(['Saving ' dFileDate ' ' interventionType_Orig]);
    runDate = datetime;
    save(fullfile(dBasePath, 'NIH Saved Data', [dFileDate ' ' interventionType_Orig 'V' version '.mat']));
    clearvars runDate
    disp(['Saved ' dFileDate ' ' interventionType_Orig]);
end