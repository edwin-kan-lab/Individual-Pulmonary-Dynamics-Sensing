function [mag, ph] = cvi(mag, ph, fs, opts)
    ncs_comp = mag.*exp(ph*1i); 
    [I_map, Q_map] = vector_injection_v3(ncs_comp, fs, opts.magSteps, ...
        opts.phaseSteps, opts.fRange, opts.optFcn, opts);
    magVec = linspace(0, 1, opts.magSteps); 
    phVec = linspace(0, 2*pi, opts.phaseSteps);

    % Magnitude
    [~, magInc] = max(max(I_map, [], 1));
    [~, phInc] = max(max(I_map, [], 2));
    mag_cross_cor_inject_phasor = magVec(magInc)*exp(phVec(phInc)*1i);
    mag_cross_cor_ncs_complex = ncs_comp + mag_cross_cor_inject_phasor;
    mag_cross_cor_ncs = abs(mag_cross_cor_ncs_complex);
    mag = rescale(mag_cross_cor_ncs, -1, 1);

    % Phase
    [~, magInc] = max(max(Q_map, [], 1));
    [~, phInc] = max(max(Q_map, [], 2));
    ph_cross_cor_inject_phasor = magVec(magInc)*exp(phVec(phInc)*1i);
    ph_cross_cor_ncs_complex = ncs_comp + ph_cross_cor_inject_phasor;
    ph_cross_cor_ncs = unwrap(angle(ph_cross_cor_ncs_complex));
    ph = rescale(ph_cross_cor_ncs, -1, 1);
end 