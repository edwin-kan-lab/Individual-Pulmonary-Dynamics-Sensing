function generateTidalVolumeResults(ncs, spiro, stage, combination, ...
    calib_method, scoring, snr_scores, lin_scores, dFileDate, toSave)
     
    ncs.tidalVolScored = cell(length(stage), 1); 
    markers = ['o', 's', '^', 'd', 'p', 'h'];

    % Calibration 
    if strcmp(calib_method, 'twopoint')
        % Linear fit with two stages - 290ml and 360ml
        for ch=1:8
            x = vertcat(ncs.tidalVol{3}(:, ch), ncs.tidalVol{4}(:, ch));
            y = vertcat(spiro.tidalVol{3}, spiro.tidalVol{4});
            coefficients = polyfit(x, y, 1);

            for i=1:length(stage)
                ncs.tidalVol{i}(:, ch) = polyval(coefficients, ncs.tidalVol{i}(:, ch));
            end 
        end 

        % PCA 
        x = vertcat(ncs.pcaTV{1}, ncs.pcaTV{2}, ncs.pcaTV{3});
        y = vertcat(spiro.tidalVol{1}, spiro.tidalVol{2}, spiro.tidalVol{3});
        coefficients = polyfit(x, y, 2);
        
        for i=1:length(stage)
            ncs.pcaTV{i} = polyval(coefficients, ncs.pcaTV{i});
        end 
    else 
        % Simple magnitude scaling - TV 360ml
        for ch=1:8
            cfac = mean(spiro.tidalVol{4})/mean(ncs.tidalVol{4}(:, ch));
            for i=1:length(stage)
                if strcmp(calib_method, 'cubic')
                    % vhat = mean(ncs.tidalVol{i}(:, ch));
                    vhat = mean(ncs.baseline{i}(:, ch));
                    % v = mean(ncs.tidalVol{4}(:, ch));
                    v = mean(ncs.baseline{4}(:, ch));
                    corrfac = 1 + (vhat/v - 1)^1/3;
                    cfac = cfac/corrfac; 
                end 
                ncs.tidalVol{i}(:, ch) = ncs.tidalVol{i}(:, ch) * cfac;
            end 
        end 

        % PCA 
        cfac = median(spiro.tidalVol{4})/median(ncs.pcaTV{4});
        for i=1:length(stage)
            ncs.pcaTV{i} = ncs.pcaTV{i} * cfac; 
        end 
    end 
    
    % Get the top-k channels 
    if strcmp(combination, 'topk')
        [s, idx] = sort(mean(snr_scores, 1), 'descend');
        sel_ch = idx(1:4); % IMPORTANT: This selects the K in top-K. Best idea is to use 3 or 4
        snr_scores = snr_scores(:, sel_ch); 
        lin_scores = lin_scores(:, sel_ch);         
    else
        sel_ch = 1:8;
    end 
    
    % Normalized quality scores
    %snr_scores = (snr_scores - min(snr_scores, [], 2))./(max(snr_scores, [], 2) - min(snr_scores, [], 2));
    %lin_scores = (lin_scores - min(lin_scores, [], 2))./(max(lin_scores, [], 2) - min(lin_scores, [], 2));

    unif_scores = ones(size(snr_scores));

    for i=1:length(stage)
        ncs.tidalVol{i} = ncs.tidalVol{i}(:, sel_ch);
        rows = size(ncs.tidalVol{i}, 1);
        scoredVols = zeros(rows, 1); 
        for j=1:rows
            if strcmp(scoring, 'snr')
                scoredVols(j) = dot(ncs.tidalVol{i}(j, :), ...
                    snr_scores(i, :))/sum(snr_scores(i, :));
            elseif strcmp(scoring, 'lin')
                scoredVols(j) = dot(ncs.tidalVol{i}(j, :), ...
                    lin_scores(i, :))/sum(lin_scores(i, :));
            else
                % uniform
                scoredVols(j) = dot(ncs.tidalVol{i}(j, :), ...
                    unif_scores(i, :))/sum(unif_scores(i, :));
            end 
        end 
        ncs.tidalVolScored{i} = scoredVols; 
    end 
    
    % Score tidal volumes 
    ncsVols = vertcat(ncs.tidalVolScored{:});
    pcaVols = vertcat(ncs.pcaTV{:});
    spiroVols = vertcat(spiro.tidalVol{:});

    ccoeff = corrcoef(spiroVols, ncsVols);
    ccoeff = ccoeff(2);
    pcaCorr = corrcoef(spiroVols, pcaVols);
    pcaCorr = pcaCorr(2); 

    tidalVolMAPE = mape(ncsVols, spiroVols);
    pcaMAPE = mape(pcaVols, spiroVols); 

    % Correlation Plot
    h = figure('Position', [800 500 600 500], 'Visible', 'off'); 
    box on;
    hold on;
    x = [];
    y = []; 
    z = []; 

    cpar = parula(2*length(stage)); 
    
    % Make pretty plots
    for st=1:length(stage)
        scatter(spiro.tidalVol{st} * 1000, ncs.tidalVolScored{st} * 1000, ...
            'DisplayName', ['TV: ' num2str(stage(st)) ' mL'], ...
            'MarkerFaceColor',  cpar(2*st-1, :), 'MarkerEdgeColor', 'k', ...
            'SizeData', 100, 'Marker', markers(mod(st, length(markers)) + 1));

        x = [x; spiro.tidalVol{st}];
        y = [y; ncs.tidalVolScored{st}];
    end 
    ax = gca; ax.Box = 'on';
    plot([ax.XLim(1) ax.XLim(2)], [ax.XLim(1) ax.XLim(2)], ...
        'LineWidth', 2, 'DisplayName', 'y = x', 'HandleVisibility', 'off')
    %plot(1000 * spiroVols, 1000 * spiroVols, 'LineWidth', 2, 'DisplayName', ...
    %    'y = x', 'HandleVisibility', 'off')
    xlabel('Spirometry Tidal Volume Estimate (mL)', 'FontSize', ...
        14, 'FontWeight', 'bold');
    ylabel('NFRF Tidal Volume Estimate (mL)', 'FontSize', 14, ...
        'FontWeight', 'bold');
    text(700 * max(spiroVols), 170, ['r^2 = ' num2str(ccoeff)], 'FontSize', 14, ...
        'FontWeight', 'bold');
    text(700 * max(spiroVols), 140, ['MAPE = ' num2str(tidalVolMAPE)], 'FontSize', ...
        14, 'FontWeight', 'bold')
    legend('FontSize', 14, 'Location', 'northwest', 'NumColumns', 2);
    if toSave
        saveas(h, strcat(dFileDate, '_', scoring, '_', calib_method, ...
            '_TV'), 'png');
    end 
    
    h2 = figure('Position', [800 500 600 500], 'Visible', 'off'); 
    box on;
    hold on;
    x = [];
    y = []; 
    z = []; 

    cpar = parula(2*length(stage)); 
    
    % PCA Plots
    for st=1:length(stage)
        scatter(spiro.tidalVol{st} * 1000, ncs.pcaTV{st} * 1000, ...
            'DisplayName', ['TV: ' num2str(stage(st)) ' mL'], ...
            'MarkerFaceColor',  cpar(2*st-1, :), 'MarkerEdgeColor', 'k', ...
            'SizeData', 100, 'Marker', markers(mod(st, length(markers)) + 1));

        x = [x; spiro.tidalVol{st}];
        z = [z; ncs.pcaTV{st}];
    end 

    ax = gca; ax.Box = 'on';
    plot([ax.XLim(1) ax.XLim(2)], [ax.XLim(1) ax.XLim(2)], ...
        'LineWidth', 2, 'DisplayName', 'y = x', 'HandleVisibility', 'off')
    %plot(1000 * spiroVols, 1000 * spiroVols, 'LineWidth', 2, 'DisplayName', ...
    %    'y = x', 'HandleVisibility', 'off')
    xlabel('Spirometry Tidal Volume Estimate (mL)', 'FontSize', ...
        14, 'FontWeight', 'bold');
    ylabel('NFRF Tidal Volume Estimate (mL)', 'FontSize', 14, ...
        'FontWeight', 'bold');
    text(700 * max(spiroVols), 170, ['r^2 = ' num2str(pcaCorr)], 'FontSize', 14, ...
        'FontWeight', 'bold');
    text(700 * max(spiroVols), 140, ['MAPE = ' num2str(pcaMAPE)], 'FontSize', ...
        14, 'FontWeight', 'bold')
    legend('FontSize', 14, 'Location', 'northwest', 'NumColumns', 2);
    
    if toSave
        saveas(h2, strcat(dFileDate, '_', calib_method, ...
            '_TV_PCA'), 'png');
    end 

    % BA Plot 
    h3 = figure('Position', [800 500 700 600], 'Visible', 'off'); 
    box on;
    hold on;
    x = [];
    y = []; 
    cpar = parula(2 * length(stage)); 
    
    % Make pretty plots
    for st=1:length(stage)
        scatter(500 * (spiro.tidalVol{st} + ncs.tidalVolScored{st}), ...
            1000*(spiro.tidalVol{st} - ncs.tidalVolScored{st}), ...
            'DisplayName', ['TV: ' num2str(stage(st)) ' mL'], ...
            'MarkerFaceColor',  cpar(2*st-1, :), 'MarkerEdgeColor', 'k', ...
            'SizeData', 100, 'Marker', markers(mod(st, length(markers)) + 1));

        x = [x; spiro.tidalVol{st}];
        y = [y; ncs.tidalVolScored{st}];
    end 
    
    ax = gca; ax.Box = 'on';
    all_samples = 1000*(spiroVols - ncsVols); 
    all_mean = mean(all_samples); all_std = std(all_samples);

    plot([ax.XLim(1) ax.XLim(2)], [all_mean all_mean], 'k', ...
        'LineWidth', 2, 'HandleVisibility','off');  
    plot([ax.XLim(1) ax.XLim(2)], [all_mean-1.96 * all_std all_mean-1.96 * all_std], ...
        'k', 'LineWidth', 2, 'LineStyle', '--', 'HandleVisibility','off');
    plot([ax.XLim(1) ax.XLim(2)], [all_mean+1.96 * all_std all_mean+1.96 * all_std], ...
        'k', 'LineWidth', 2, 'LineStyle', '--', 'HandleVisibility','off');

    text(120, all_mean + 1.96 * all_std + 4, ...
        ['Mean + 1.96 SD: ' num2str(round(all_mean + 1.96 * all_std, 2)) ' mL'], ...
        'FontSize', 14, 'FontWeight', 'bold');
    text(120, all_mean - 6, ['Mean: ' num2str(round(all_mean, 2)) ' mL'], ...
        'FontSize', 14, 'FontWeight', 'bold');
    text(120, all_mean - 1.96 * all_std + 5, ...
        ['Mean - 1.96 SD: ' num2str(round(all_mean - 1.96 * all_std, 2)) ' mL'], ...
        'FontSize', 14, 'FontWeight', 'bold');
    legend(ax, 'location', 'southeast', 'FontSize', 15, 'NumColumns', 3);
    ylabel('Volume Difference (mL)', 'FontWeight','bold', 'FontSize', 17);
    xlabel('Mean Volume (mL)', 'FontWeight','bold', 'FontSize', 17);
    ax.YLim = [all_mean - 3.5 * all_std all_mean + 3 * all_std]; 
    
    if toSave
        saveas(h3, strcat(dFileDate, '_', scoring, '_', calib_method, ...
            '_TV_BA'), 'png');
    end
    
    % PCA - BA
    h4 = figure('Position', [800 500 700 600], 'Visible', 'off'); 
    box on;
    hold on;
    x = [];
    z = []; 
    cpar = parula(2*length(stage)); 
    
    % Make pretty plots
    for st=1:length(stage)
        scatter(500 * (spiro.tidalVol{st} + ncs.pcaTV{st}), ...
            1000*(spiro.tidalVol{st} - ncs.pcaTV{st}), ...
            'DisplayName', ['TV: ' num2str(stage(st)) ' mL'], ...
            'MarkerFaceColor',  cpar(2*st-1, :), 'MarkerEdgeColor', 'k', ...
            'SizeData', 100, 'Marker', markers(mod(st, length(markers)) + 1));

        x = [x; spiro.tidalVol{st}];
        z = [z; ncs.pcaTV{st}];
    end 
    
    ax = gca; ax.Box = 'on';
    all_samples = 1000*(spiroVols - pcaVols); 
    all_mean = mean(all_samples); all_std = std(all_samples);

    plot([ax.XLim(1) ax.XLim(2)], [all_mean all_mean], 'k', ...
        'LineWidth', 2, 'HandleVisibility','off');  
    plot([ax.XLim(1) ax.XLim(2)], [all_mean-1.96 * all_std all_mean-1.96 * all_std], ...
        'k', 'LineWidth', 2, 'LineStyle', '--', 'HandleVisibility','off');
    plot([ax.XLim(1) ax.XLim(2)], [all_mean+1.96 * all_std all_mean+1.96 * all_std], ...
        'k', 'LineWidth', 2, 'LineStyle', '--', 'HandleVisibility','off');

    text(120, all_mean + 1.96 * all_std + 4, ...
        ['Mean + 1.96 SD: ' num2str(round(all_mean + 1.96 * all_std, 2)) ' mL'], ...
        'FontSize', 14, 'FontWeight', 'bold');
    text(120, all_mean - 6, ['Mean: ' num2str(round(all_mean, 2)) ' mL'], ...
        'FontSize', 14, 'FontWeight', 'bold');
    text(120, all_mean - 1.96 * all_std + 5, ...
        ['Mean - 1.96 SD: ' num2str(round(all_mean - 1.96 * all_std, 2)) ' mL'], ...
        'FontSize', 14, 'FontWeight', 'bold');
    legend(ax, 'location', 'southeast', 'FontSize', 15, 'NumColumns', 3);
    ylabel('Volume Difference (mL)', 'FontWeight','bold', 'FontSize', 17);
    xlabel('Mean Volume (mL)', 'FontWeight','bold', 'FontSize', 17);
    ax.YLim = [all_mean - 3.5 * all_std all_mean + 3 * all_std]; 
    
    if toSave
        saveas(h4, strcat(dFileDate, '_', calib_method, ...
        '_TV_BA_PCA'), 'png');

        % Save for further analysis 
        disp(['Saving TV for ' calib_method ' and ' scoring]);
        save(strcat(dFileDate, '_TV_', scoring, '_', calib_method, '.mat'), 'ncsVols', 'pcaVols', 'spiroVols');
    end 
end