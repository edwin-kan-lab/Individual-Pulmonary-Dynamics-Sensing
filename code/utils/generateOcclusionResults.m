function generateOcclusionResults(opts, overall_start_loc, spiro, end_times, ...
    date_, start_true, deflate_true, cpap_true, ncs, calib_st, calib_en, ...
    ncs_sig, dFileDate, scoring)
    
    toSave = true; 
    ncs_sig = ncs.ncsRespFilt;
    
    % Time start and end for time series signal plots 
    l1_loc = overall_start_loc + opts.win - 1;
    l2_loc = getLocFromTime(spiro.locs, end_times(end), date_.seconds_offset, 100);
    l1 = ncs.locs(l1_loc); % Corrected
    l2 = ncs.locs(l2_loc);
    
    gt_start_loc = getLocFromTime(spiro.locs, start_true, ...
        date_.seconds_offset, 100) - l1_loc;
    gt_deflate_loc = getLocFromTime(spiro.locs, deflate_true, ...
        date_.seconds_offset, 100) - l1_loc;
    gt_cpap_loc = getLocFromTime(spiro.locs, cpap_true, ...
        date_.seconds_offset, 100) - l1_loc;
    
    ncsTVs = vertcat(ncs.tidalVol{:}); % breaths x channels 
    
    % Calibrate - Single Point 
    for i=1:8
        cfac = mean(spiro.tidalVol{1}(1:10))/mean(ncsTVs(1:10, i)); 
        ncsTVs(:, i) = ncsTVs(:, i) * cfac; 
    end 
    
    % Compute scores using a given calibration time 
    scores = zeros(8, 1); 
    cst = getLocFromTime(spiro.locs, calib_st, date_.seconds_offset, 100);
    cen = getLocFromTime(spiro.locs, calib_en, date_.seconds_offset, 100);
    for ch=1:8 
        scores(ch) = getSNROcclusion(ncs_sig(ncs.locs(cst): ncs.locs(cen), ch), ...
                1000, [0.18 0.4]);
    end 
    
    if strcmp(scoring, 'lr')
        % Only use 900 MHz
        left_tv = (ncsTVs(:, 3) * scores(3) + ncsTVs(:, 4) * scores(4))/...
            (scores(3) + scores(4));
        right_tv = (ncsTVs(:, 1) * scores(1) + ncsTVs(:, 2) * scores(2))/...
            (scores(1) + scores(2));
        left_ncs = (normalize(ncs_sig(:, 7), 'range') * scores(3) + ...
            normalize(ncs_sig(:, 8), 'range') * scores(4))/(scores(3) + scores(4));
        right_ncs = (normalize(ncs_sig(:, 1), 'range') * scores(1) + ...
            normalize(ncs_sig(:, 2), 'range') * scores(2))/(scores(1) + scores(2));
    elseif strcmp(scoring, 'lr2')
        % Only use 2 GHz
        left_tv = (ncsTVs(:, 7) * scores(7) + ncsTVs(:, 8) * scores(8))/...
            (scores(7) + scores(8));
        right_tv = (ncsTVs(:, 5) * scores(5) + ncsTVs(:, 6) * scores(6))/...
            (scores(5) + scores(6));
        left_ncs = (normalize(ncs_sig(:, 15), 'range') * scores(7) + ...
            normalize(ncs_sig(:, 16), 'range') * scores(8))/(scores(7) + scores(8));
        right_ncs = (normalize(ncs_sig(:, 9), 'range') * scores(5) + ...
            normalize(ncs_sig(:, 10), 'range') * scores(6))/(scores(5) + scores(6));
    else
        % Use both 
        left_tv = (ncsTVs(:, 7) * scores(7) + ncsTVs(:, 8) * scores(8)...
            + ncsTVs(:, 3) * scores(3) + ncsTVs(:, 4) * scores(4))/(scores(7) ...
            + scores(8) + scores(3) + scores(4));
        right_tv = (ncsTVs(:, 1) * scores(1) + ncsTVs(:, 2) * scores(2)...
            + ncsTVs(:, 5) * scores(5) + ncsTVs(:, 6) * scores(6))/(scores(1) ...
            + scores(2) + scores(5) + scores(6));
        left_ncs = (normalize(ncs_sig(:, 15), 'range') * scores(7) + ...
            normalize(ncs_sig(:, 16), 'range') * scores(8) + ...
            normalize(ncs_sig(:, 7), 'range') * scores(3) + ...
            normalize(ncs_sig(:, 8), 'range') * scores(4))/(scores(7) + scores(8) + scores(3) + scores(4));
        right_ncs = (normalize(ncs_sig(:, 9), 'range') * scores(5) + ...
            normalize(ncs_sig(:, 10), 'range') * scores(6) + ...
            normalize(ncs_sig(:, 1), 'range') * scores(1) + ...
            normalize(ncs_sig(:, 2), 'range') * scores(2))/(scores(1) + scores(2) + scores(5) + scores(6));
    end 
    
    % These lengths will be L2_Loc - L1_Loc
    lr = left_tv./right_tv;
    
    % Detect Occlusion using thresholding on moving average 
    [a, b] = getOcclusionPointsV3(lr); 

    % Error Handling: No point detected
    if length(a) == 0
        a = 1; 
    end 
    if length(b) == 0
        b = length(lr_comb);
    end 
    
    % Figure 1
    h = figure('Position', [600 300 850 1000]); 
    box on; 
    t = tiledlayout(2,1,'TileSpacing', 'compact');
    cpar = parula(100);
    ax1 = nexttile;
    plot(left_tv, 'LineWidth', 2, 'Color', cpar(45, :), ...
        'DisplayName', 'Left Lung TV'); hold on;
    plot(right_tv, 'LineWidth', 2, 'Color', cpar(65, :), ...
        'DisplayName', 'Right Lung TV');
    title('NFRF Tidal Volume Estimate (mL)', 'FontSize', 14, ...
        'FontWeight', 'bold');
    legend('location', 'north', 'FontSize', 13);
    % Fixed - Use V2 points
    xline(a, '-', {'Predicted Start'}, ...
        'LineWidth', 2, 'Color', [0.33 0.07 0.39], 'HandleVisibility','off', ...
        'FontSize', 14, 'FontWeight', 'bold');
    xline(b, '-', {'Predicted End'}, ...
        'LineWidth', 2, 'Color', 'k', 'HandleVisibility','off', ...
        'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Tidal Volume (mL)', 'FontSize', 14, 'FontWeight', 'bold');

    ax2 = nexttile; 
    plot(lr, 'LineWidth', 2, 'Color', [0.85 0.33 0.10]); 
    title('NFRF Sided Tidal Volume Ratio', 'FontSize', 14, 'FontWeight', 'bold');
    xline(gt_start_loc, '-', {'Blocker Inflated'}, ...
        'LineWidth', 2, 'Color', [0.33 0.07 0.39], ...
        'HandleVisibility','off', 'FontSize', 14, 'FontWeight', 'bold');
    xline(gt_deflate_loc, '-', {'Blocker Deflated'}, ...
        'LineWidth', 2, 'Color', [0.15 0.54 0.05], ...
        'HandleVisibility','off', 'FontSize', 14, 'FontWeight', 'bold');
    xline(gt_cpap_loc, '-', {' '}, ...
        'LineWidth', 2, 'Color', [0 0 0], ...
        'HandleVisibility','off', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Volume Ratio', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Breaths', 'FontSize', 14, 'FontWeight', 'bold');
    linkaxes([ax1, ax2], 'x');
    xlim([1 length(lr)]);
    
    if toSave
        saveas(h, strcat(dFileDate, '_', scoring, '_OcclusionLR'), 'png');
    end 

    % Figure 2 
    h = figure('Position', [500 100 1300 1050]); 
    box on;
    
    % Make time series look the same 
    mu_left = mean(left_ncs(ncs.locs(cst):ncs.locs(cen)));
    mu_right = mean(right_ncs(ncs.locs(cst):ncs.locs(cen)));
    pkpk_left = peak2peak(left_ncs(ncs.locs(cst):ncs.locs(cen)));
    pkpk_right = peak2peak(right_ncs(ncs.locs(cst):ncs.locs(cen)));
    right_ncs = (right_ncs - mu_right)./pkpk_right; 
    left_ncs = (left_ncs - mu_left)./pkpk_left; 

    times = ncs.tNcsDS(l1:l2);
    timeAx = ncs.tNcsDS(1:end-1) - ncs.tNcsDS(1); 
    cpar = parula(100);
    ax1 = subplot('Position', [0.1 0.72 0.85 0.24]);

    plot(timeAx, left_ncs, 'LineWidth', 1.5, 'Color', cpar(45, :), ...
        'DisplayName', 'Left Lung TV'); hold on;
    title('NFRF Left Lung Signal', 'FontSize', 17, ...
        'FontWeight', 'bold');
    
    ht = ax1.YLim(2)-ax1.YLim(1);
    base = ax1.YLim(1);

    xline(timeAx(ncs.locs(a(1)+l1_loc)), '-', {}, ...
        'LineWidth', 2, 'Color', [0.59 0.13 0.21], ...
        'HandleVisibility','off', 'FontSize', 16, 'FontWeight', 'bold');
    ll = min(length(a), length(b));
    xline(timeAx(ncs.locs(b(ll)+l1_loc)), '-', {}, ...
        'LineWidth', 2, 'Color', [0.59 0.13 0.21], ...
        'HandleVisibility','off', 'FontSize', 16, 'FontWeight', 'bold');
    
    ylim1 = [base base+ht];
    ylabel('Amplitude (a.u.)', 'FontSize', 16, 'FontWeight', 'bold');
    xlabel('Time (s)', 'FontSize', 16, 'FontWeight', 'bold');
    
    ax1.XAxis.FontSize = 15;
    ax1.YAxis.FontSize = 15;
    
    ax2 = subplot('Position', [0.1 0.39 0.85 0.24]);
    ax2.FontSize = 16; 
    plot(timeAx, right_ncs, 'LineWidth', 1.5, 'Color', cpar(65, :), ...
        'DisplayName', 'Right Lung TV'); hold on;
    title('NFRF Right Lung Signal', 'FontSize', 17, ...
        'FontWeight', 'bold');
    ht = ax2.YLim(2)-ax2.YLim(1);
    base = ax2.YLim(1);
    xlabel('Time (s)', 'FontSize', 16, 'FontWeight', 'bold');

    xline(timeAx(ncs.locs(a(1)+l1_loc)), '-', {}, ...
        'LineWidth', 2, 'Color', [0.59 0.13 0.21], ...
        'HandleVisibility','off', 'FontSize', 16, 'FontWeight', 'bold');
    ll = min(length(a), length(b));
    xline(timeAx(ncs.locs(b(ll)+l1_loc)), '-', {}, ...
        'LineWidth', 2, 'Color', [0.59 0.13 0.21], ...
        'HandleVisibility','off', 'FontSize', 16, 'FontWeight', 'bold');
    
    ylim2 = [base base+ht];
    ylabel('Amplitude (a.u.)', 'FontSize', 16, 'FontWeight', 'bold');
    
    linkaxes([ax1, ax2], 'xy');
    
    ax2.XLim = [timeAx(l1) timeAx(l2)];
    ax2.YLim = (ylim1 + ylim2)./2;
    
    ax2.XAxis.FontSize = 15;
    ax2.YAxis.FontSize = 15;

    ax3 = subplot('Position', [0.1 0.07 0.85 0.23]);
    plot(lr, 'LineWidth', 2, 'Color', [0.85 0.33 0.10]); 
    title('NFRF Sided Tidal Volume Ratio', 'FontSize', 17, 'FontWeight', 'bold');
    xline(gt_start_loc, '-', {'Blocker', 'Inflated'}, ...
        'LineWidth', 2, 'Color', [0.33 0.07 0.39], ...
        'HandleVisibility','off', 'FontSize', 16, 'FontWeight', 'bold');
    xline(gt_deflate_loc, '-', {'Blocker', 'Deflated'}, ...
        'LineWidth', 2, 'Color', [0.15 0.54 0.05], ...
        'HandleVisibility','off', 'FontSize', 16, 'FontWeight', 'bold');
    xline(gt_cpap_loc, '-', {'CPAP'}, ...
        'LineWidth', 2, 'Color', [0 0 0], ...
        'HandleVisibility','off', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('Volume Ratio', 'FontSize', 16, 'FontWeight', 'bold');
    xlabel('Breaths', 'FontSize', 16, 'FontWeight', 'bold');
    xlim([1 length(lr)]);
    
    ax3.XAxis.FontSize = 15;
    ax3.YAxis.FontSize = 15;
    
    if toSave
        saveas(h, strcat(dFileDate, '_', scoring, '_OcclusionOverall'), 'png');
        % Save for further analysis 
        disp('Saving Occlusion Data');
        save(strcat(dFileDate, '_', scoring, '_Occlusion.mat'), 'a', 'b', ...
            'gt_start_loc', 'gt_deflate_loc', 'gt_cpap_loc', 'start_true', ...
            'deflate_true', 'cpap_true'); 
    end 
end