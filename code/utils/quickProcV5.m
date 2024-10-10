% quickProc.m
dPath = 'C:\Users\Kan Lab\Documents\Pig Study\Data\';
d = dir(dPath);
NCSFile = 'CasetestNCSSaving0124_164147.mat';
BIOFile = strrep(NCSFile, 'NCS', 'BIO');
RadioFile = strrep(NCSFile, 'NCS', 'Radio');

if isfile([dPath NCSFile])
    disp('Loading');
    load([dPath NCSFile], 'NcsData1', 'NcsData2');
else
    disp('Converting');
    saveTDMStoMAT_NCSpigV5(dPath, strrep(NCSFile, '.mat', ''));
    disp('Loading');
    load([dPath NCSFile], 'NcsData1', 'NcsData2');
end

if isfile([dPath BIOFile])
    disp('Loading');
    load([dPath BIOFile], 'BioData');
else
    disp('Converting');
    saveTDMStoMAT_BIOpig(dPath, strrep(BIOFile, '.mat', ''));
    disp('Loading');
    load([dPath BIOFile], 'BioData');
end

if isfile([dPath RadioFile])
    disp('Loading');
    load([dPath RadioFile], 'RadioData');
else
    disp('Converting');
    saveTDMStoMAT_RadiopigV5(dPath, strrep(RadioFile, '.mat', ''));
    disp('Loading');
    load([dPath RadioFile], 'RadioData');
end
disp('Complete');
%% Sync
ncs1Ch = 1; ncs2Ch = 1;

% Sync NCS and BIO
p = 7.3306e-05; tOffset = 0.045;
fNcs = 10000; fBio = 1000/(1-p);
tNcs = 1/fNcs:1/fNcs:length(NcsData1)/fNcs;
tBio = (1/fBio:1/fBio:length(BioData)/fBio)+tOffset;
tRadio = linspace(0, max(tNcs), size(RadioData, 1)); fRadio = 1/(tRadio(2)-tRadio(1));

% Trim
trim = [30 max(tNcs)-30];
ncs1 = NcsData1(round(trim(1)*fNcs):round(trim(2)*fNcs), :);
ncs2 = NcsData2(round(trim(1)*fNcs):round(trim(2)*fNcs), :);
bio = BioData(round(trim(1)*fBio):round(trim(2)*fBio), :);
radio = RadioData(round(trim(1)*fRadio):round(trim(2)*fRadio), :);
tNcs = tNcs(round(trim(1)*fNcs):round(trim(2)*fNcs));
tBio = tBio(round(trim(1)*fBio):round(trim(2)*fBio));
tRadio = tRadio(round(trim(1)*fRadio):round(trim(2)*fRadio));

%% Plot
figure;
ax(1)=subplot(311); 
yyaxis left; 
plot(tNcs(1:10:end), (ncs1(1:10:end, ncs1Ch)), 'k', 'DisplayName', 'Amp'); hold on; 
yyaxis right;
plot(tNcs(1:10:end), -(unwrap(ncs1(1:10:end, ncs1Ch+1))), 'b', 'DisplayName', 'Phs');
legend('Location', 'south', 'Orientation', 'horizontal');

ax(2)=subplot(312); 
yyaxis left; 
plot(tNcs(1:10:end), (ncs2(1:10:end, ncs2Ch)), 'k', 'DisplayName', 'Amp'); hold on; 
yyaxis right;
plot(tNcs(1:10:end), -(unwrap(ncs2(1:10:end, ncs2Ch+1))), 'b', 'DisplayName', 'Phs');
legend('Location', 'south', 'Orientation', 'horizontal');

ax(3)=subplot(313); 
plot(tBio, bio(:, 1), 'k', 'DisplayName', 'ECG'); hold on; 
% yyaxis left; plot(tBio, bio(:, 2), 'b', 'DisplayName', 'PA pressure'); hold on; 
% yyaxis right; plot(tBio, bio(:,3), 'Color', [1 0.5 0], 'DisplayName', 'AO pressure');
legend('Location', 'south', 'Orientation', 'horizontal');

linkaxes(ax, 'x');

%% FFT
% [P1, f] = ncsFFT(ncs(:, ncsCh), fNcs, true, [0.5 10]);

