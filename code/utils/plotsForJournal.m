%% Make Airflow Plots for showing airflow morphology
f = figure('Position', [800 500 800 600]); 
ax = gca(f);
box on;
sxx = 3; 
start_loc = getLocFromTime(spiro.locs, start_times(sxx), date_.seconds_offset, 100) + 1;
end_loc = getLocFromTime(spiro.locs, end_times(sxx), date_.seconds_offset, 100);

num_locs = end_loc - start_loc; 

end_loc = spiro.locs(end_loc);
tdiff = 100; 

template = zeros(1000, num_locs); 
c = 0.7;
for i=1:num_locs
    stId = spiro.locs(start_loc + i -1) - tdiff;
    enId = spiro.locs(start_loc + i) - tdiff; 

    tmp = spiro.airflow(stId:enId); 
    tmp = spline(1:length(tmp), tmp, linspace(1, length(tmp), 1000)); 
    plot(tmp, 'Color', [c c c]); hold on; 
    template(:, i) = tmp;
end 

mean_tmp = mean(template, 2); 
plot(mean_tmp, 'Color', '#0072BD', 'LineWidth', 3);
ylabel('Flow Rate (a.u.)', 'FontWeight','bold', 'FontSize', 18);
xlabel('Sampling Points', 'FontWeight','bold', 'FontSize', 18);
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
saveas(f, 'airflow', 'png');

%% Make BA
close all; 
clear all;
dates = {'02-05-24', '02-06-24', '02-07-24', '02-08-24', '03-04-24', ...
    '03-05-24', '03-06-24', '03-07-24'}; 
scoring = 'SNR';
calib_method = 'single'; 

predvols = []; 
tgtvols = []; 

for runAllInd = 1:length(dates)
    dFileDate = dates{runAllInd}; 
    load(strcat(dFileDate, '_TV_', scoring, '_', calib_method, '.mat'));
    
    predvols = [predvols; ncsVols];
    tgtvols = [tgtvols; spiroVols];
end
[val, idx] = sort(tgtvols, 'ascend'); 
tgtvols = tgtvols(idx); 
predvols = predvols(idx); 

ccoeff = corrcoef(tgtvols, predvols);
ccoeff = round(ccoeff(2), 2);
tidalVolMAPE = round(mape(predvols, tgtvols), 2);

% Correlation Plot
h = figure('Position', [800 500 700 500]); 
ax = gca(h);
box on;
hold on;
cpar = parula(3*length(tgtvols)); 
stcol = floor(1.2*length(tgtvols)); 
cpar = cpar(stcol:stcol+length(tgtvols)-1, :);

scatter(tgtvols * 1000, predvols * 1000, [], cpar, 'filled', 'Marker', 'o', ...
    'HandleVisibility', 'off');
ax = gca; ax.Box = 'on';
ax.XLim = [120 620];
ax.YLim = [60 730];
plot([ax.XLim(1) ax.XLim(2)], [ax.XLim(1) ax.XLim(2)], 'LineWidth', 2, ...
    'DisplayName', 'y = x', 'Color', [0.49 0.18 0.56], 'LineStyle', '--');
%plot(1000 * tgtvols, 1000 * tgtvols, 'LineWidth', 2, 'DisplayName', ...
%        'y = x', 'HandleVisibility', 'off')
xlabel('Spirometry Tidal Volume Estimate (mL)', 'FontSize', ...
        17, 'FontWeight', 'bold');
ylabel('NFRF Tidal Volume Estimate (mL)', 'FontSize', 17, ...
        'FontWeight', 'bold');
text(498, 232, ['r^2 = ' num2str(ccoeff)], 'FontSize', 16, ...
        'FontWeight', 'bold', 'Color', [0.49 0.18 0.56]);
text(463, 157, ['MAPE = ' num2str(tidalVolMAPE)], 'FontSize', ...
        16, 'FontWeight', 'bold', 'Color', [0.25 0.53 0.04])
legend('FontSize', 16, 'Location', 'northwest', 'NumColumns', 2);
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;

% BA Plots
h2 = figure('Position', [800 500 800 500], 'Visible', 'on'); 
box on;
hold on;
 
cpar = parula(3*length(tgtvols)); 
stcol = floor(1.2*length(tgtvols)); 
cpar = cpar(stcol:stcol+length(tgtvols)-1, :);
ax = gca; ax.Box = 'on';
all_samples = 1000*(tgtvols - predvols); 
all_mean = mean(all_samples); all_std = std(all_samples);

scatter(500 * (tgtvols + predvols), 1000*(tgtvols - predvols), [], cpar, ...
            'filled', 'Marker', 'o', 'HandleVisibility', 'off');
ax.XLim = [100 670];
plot([ax.XLim(1) ax.XLim(2)], [all_mean all_mean], 'k', ...
        'LineWidth', 2, 'Color', [0.49 0.18 0.56], 'HandleVisibility','off', ...
        'Color', [0.64 0.08 0.18]);  
plot([ax.XLim(1) ax.XLim(2)], [all_mean-1.96 * all_std all_mean-1.96 * all_std], ...
        'k', 'LineWidth', 2, 'LineStyle', '--', 'HandleVisibility','off', ...
        'Color', [0.49 0.18 0.56]);
plot([ax.XLim(1) ax.XLim(2)], [all_mean+1.96 * all_std all_mean+1.96 * all_std], ...
        'k', 'LineWidth', 2, 'LineStyle', '--', 'HandleVisibility','off', ...
        'Color', [0.64 0.08 0.18]);

text(120, 145, ['Mean + 1.96 SD: ' num2str(round(all_mean + 1.96 * all_std, 2)) ' mL'], ...
        'FontSize', 15, 'FontWeight', 'bold', 'Color', [0.64 0.08 0.18]);
text(529, 3, ['Mean: ' num2str(round(all_mean, 2)) ' mL'], ...
        'FontSize', 15, 'FontWeight', 'bold', 'Color', [0.49 0.18 0.56]);
text(120, -119, ['Mean - 1.96 SD: ' num2str(round(all_mean - 1.96 * all_std, 2)) ' mL'], ...
        'FontSize', 15, 'FontWeight', 'bold', 'Color', [0.64 0.08 0.18]);
%legend(ax, 'location', 'southeast', 'FontSize', 15, 'NumColumns', 3);
ylabel('Volume Difference (mL)', 'FontWeight','bold', 'FontSize', 17);
xlabel('Mean Volume (mL)', 'FontWeight','bold', 'FontSize', 17);
ax.YLim = [all_mean - 3.5 * all_std all_mean + 3 * all_std]; 
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;

%% Make Tidal Volume - Spiro Plots
% Considering file 02-05 tidal vol loaded using main script 
f = figure('Position', [800 500 700 600]); 
box on;
sxx = 3; 
start_loc = 120;
end_loc = 125;

num_locs = end_loc - start_loc; 

end_loc = spiro.locs(end_loc);
tdiff = 100; 

template = zeros(1000, num_locs); 
ncs_template = zeros(1000, num_locs); 

c = 0.7;
ax1 = subplot(211);
for i=1:num_locs
    stId = spiro.locs(start_loc + i -1) - tdiff;
    enId = spiro.locs(start_loc + i) - tdiff; 
 
    tmp = spiro.volume(stId:enId);
    % Process to make it tidal volume 

    tmp = spline(1:length(tmp), tmp, linspace(1, length(tmp), 1000)); 
    plot(tmp, 'Color', [c c c]); hold on; 
    template(:, i) = tmp;
end 

mean_tmp = mean(template, 2); 
plot(mean_tmp, 'Color', '#800080', 'LineWidth', 3);
ylabel('Tidal Volume (L)', 'FontWeight','bold', 'FontSize', 15);
xlabel('Sampling Points', 'FontWeight','bold', 'FontSize', 15);
title('Spirometer Measurement', 'FontWeight','bold', 'FontSize', 16);
ylim([min(mean_tmp)-0.02 max(mean_tmp)+0.03])
%xticklabels([])
ax1.XAxis.FontSize = 14;
ax1.YAxis.FontSize = 14;

tdiff2 = 1000;
ax2 = subplot(212); 
for i=1:num_locs
    stId = ncs.locs(start_loc + i -1) - tdiff2;
    enId = ncs.locs(start_loc + i) - tdiff2; 
 
    tmp = ncs.ncsRespFilt(stId:enId, 8);
    % Process to make it tidal volume 

    tmp = spline(1:length(tmp), tmp, linspace(1, length(tmp), 1000)); 
    plot(tmp, 'Color', [c c c]); hold on; 
    ncs_template(:, i) = tmp;
end 

mean_ncs_tmp = mean(ncs_template, 2); 
plot(mean_ncs_tmp, 'Color', '#15B01A', 'LineWidth', 3);
ylabel('Amplitude (a.u.)', 'FontWeight','bold', 'FontSize', 15);
xlabel('Sampling Points', 'FontWeight','bold', 'FontSize', 15);
%xticklabels([]);
title('NFRF Measurement', 'FontWeight','bold', 'FontSize', 16);
ylim([min(mean_ncs_tmp)-0.2 max(mean_ncs_tmp)+0.25]);
ax2.XAxis.FontSize = 14;
ax2.YAxis.FontSize = 14;
saveas(f, 'ncs_spiro', 'png');
