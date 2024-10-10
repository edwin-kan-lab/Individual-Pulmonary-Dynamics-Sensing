function [ncsVols, pcaVol, ncsBase, spiroVols, snrs, lins] = ...
    getTidalVolumesV3(ncs_sig, spiro_sig, opts)
    
    % Find locs for stage
    spiro_locs = opts.spiro_locs; 
    ncs_locs = opts.ncs_locs; 
    
    if opts.complex
        jump = 2;
    else
        jump = 1;
    end 
    
    % For PCA
    num_comps = opts.num_comps; 
    [coeff,score,latent,tsquared,explained,mu] = pca(ncs_sig);
    
    pca_sig = score(:, 1:num_comps); 

    snrs = zeros(1, length(opts.sel_channels));
    lins = zeros(1, length(opts.sel_channels));

    for ch=1:jump:length(opts.sel_channels)
        snrs(ch) = getSNR(ncs_sig(ncs_locs(1): ncs_locs(end), ch), ...
                1000, [0.18 0.35]);
        lins(ch) = getLinearity(ncs_sig(ncs_locs(1): ncs_locs(end), ch), ...
                 1000, [0.18 0.35]);
    end

    if length(spiro_locs) <= opts.win
        % In this case, just compute for each breath since we want error
        % bars 
        spiroVols = zeros(length(spiro_locs)-1, 1);
        spiroBase = zeros(length(spiro_locs)-1, 1);
        pcaVol = zeros(length(spiro_locs)-1, num_comps);
        ncsVols = zeros(length(spiro_locs)-1, length(opts.sel_channels));
        ncsBase = zeros(length(spiro_locs)-1, length(opts.sel_channels));
        
        for i=1:length(spiro_locs)-1
            sig = spiro_sig(spiro_locs(i): spiro_locs(i+1));
            spiroVols(i) = peak2peak(sig);
            
            % Needed for complex peak-peak
            if opts.complex
                [~, st] = max(sig); [~, en] = min(sig);
            end 

            % PCA Vol 
            for nc=1:num_comps
                sig = pca_sig(ncs_locs(i): ncs_locs(i+1), nc);
                pcaVol(i, nc) = peak2peak(sig);
            end

            for ch=1:jump:length(opts.sel_channels)
                if opts.complex
                    compl = ncs_sig(ncs_locs(i): ncs_locs(i+1), ch).*exp(ncs_sig(ncs_locs(i): ncs_locs(i+1), ch+1)*1i);
                    compl = movmean(compl, 10); 
                    ncsVols(i, ch) = abs(compl(st) - compl(en));
                    ncsVols(i, ch + 1) = 0;
                else
                    ncsVols(i, ch) = peak2peak(ncs_sig(ncs_locs(i): ncs_locs(i+1), ch));
                end

                ncsBase(i, ch) = mean(ncs_sig(ncs_locs(i): ncs_locs(i+1), ch));
                % Needed for complex Pk-Pk 
                if opts.complex
                    ncsBase(i, ch + 1) = mean(ncsTemp(ch + 1, :));
                end 
            end 
        end 
    else
        num_breaths = length(spiro_locs) - opts.win + 1;
        spiroVols = zeros(num_breaths-1, 1);
        spiroBase = zeros(num_breaths-1, 1);
        pcaVol = zeros(num_breaths-1, num_comps);
        ncsVols = zeros(num_breaths-1, length(opts.sel_channels));
        ncsBase = zeros(num_breaths-1, length(opts.sel_channels));

        for i=1:num_breaths-1
            sig = getTemplate(spiro_sig, spiro_locs, i+opts.win, ...
                opts.win, opts.seg);
            spiroVols(i) = peak2peak(sig);
            spiroBase(i) = mean(sig); 
            
            % Needed for complex peak-peak
            if opts.complex
                [~, st] = max(sig); [~, en] = min(sig);
            end 
            
            % PCA Vol 
            sig = getMultiTemplate(pca_sig, ncs_locs, i+opts.win, ...
                    opts.win, opts.seg);
            for nc=1:num_comps                
                pcaVol(i, nc) = peak2peak(sig(nc, :));
            end 
            
            
            ncsTemp = getMultiTemplate(ncs_sig, ncs_locs, i+opts.win, ...
                    opts.win, opts.seg);
            
            for ch=1:jump:length(opts.sel_channels)
                if opts.complex
                    compl = ncsTemp(ch, :).*exp(ncsTemp(ch+1, :)*1i);
                    compl = movmean(compl, 10); 
                    ncsVols(i, ch) = abs(compl(st) - compl(en));
                    ncsVols(i, ch + 1) = 0;
                else
                    ncsVols(i, ch) = peak2peak(ncsTemp(ch, :));
                end
                   
                ncsBase(i, ch) = mean(ncsTemp(ch, :));
                % Needed for complex Pk-Pk 
                if opts.complex
                    ncsBase(i, ch + 1) = mean(ncsTemp(ch + 1, :));
                end 
            end 
        end
    end 
end 