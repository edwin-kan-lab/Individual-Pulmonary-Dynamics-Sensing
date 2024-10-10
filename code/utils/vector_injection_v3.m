function [mag_map, ph_map] = ...
    vector_injection_v3(ncs_complex, Fs, mag_step, ph_step, fRange, optFcn, opts)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code is used for performing vector injection with a given objective
% function to optimize a certain signal feature or regularize signal 
% morphology.
%
% % INPUTS:
% ncs_complex:   Complex quadrature signal array.
% Fs:       I/Q magnitude/phase signal sampling rate.
% mag_step: Injected vector magnitude scan steps.
% ph_step:  Injected vector phase scan steps.
% fRange:   Target subject motion frequency range.
% tempType: Template type: 'square'or 'sine'
% optFcn: 'SNR': SNR optimization, 'lin':1st/2nd harm optimization
%         'temp': template similarity optimization
%
% opts: options for each objective function that define behavior
% SNR: .nfRange: noise floor range
% temp: .tempType: template type | .similarity: xcorr or MSE
%       .winSecs: seconds searched | .absolute: t/f for absolute value
%
% Jianlin Zhou | Cornell University | jz899@cornell.edu
% Edited: Thomas Conroy: 5/5/2022
% Edited for greater parameter flexibility and choice of optimization
% functions rather than always running all
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Defining default options
if strcmp(optFcn, 'SNR')
    if ~isfield(opts,'nfRange')
        opts.nfRange = [10 50];
        fprintf('No nfRange specified, using [10 50].\n');
    end
elseif strcmp(optFcn, 'temp')
    if ~isfield(opts,'tempType')
        opts.tempType = 'square';
        fprintf('No tempType specified, using square.\n');
    end
    if ~isfield(opts,'similarity')
        opts.similarity = 'xcorr';
        fprintf('No similarity specified, using xcorr.\n');
    end
    if ~isfield(opts,'winSecs')
        opts.winSecs = 'xcorr';
        fprintf('No winSecs specified, using 1 second.\n');
    end
    if ~isfield(opts,'absolute')
        opts.absolute = true;
        fprintf('No absolute specified, using absolute value.\n');
    end
end

mag_scan = linspace(0,0.9,mag_step);
ph_scan = linspace(0,2*pi,ph_step);

for i_mag = 1 : length(mag_scan)
    for i_ph = 1 : length(ph_scan)

        ncs_complex_injection_step = mag_scan(i_mag)*exp(ph_scan(i_ph)*1i);
        ncs_complex_imitation_orig_scan = ncs_complex + ncs_complex_injection_step;

        % Magnitude analysis
        complex_sig_mag = detrend(abs(ncs_complex_imitation_orig_scan));
        % Phase analysis
        complex_sig_phase = detrend(unwrap(angle(ncs_complex_imitation_orig_scan)));
        
        if strcmp(optFcn, 'SNR')
            % Objective Function 1 Fundamental tone SNR
            mag_tmp = func_obj_fundamentalTone(complex_sig_mag,Fs,fRange,opts.nfRange);
            ph_tmp = func_obj_fundamentalTone(complex_sig_phase,Fs,fRange,opts.nfRange);
        elseif strcmp(optFcn, 'lin')
            % Objective Function 2 Fundamental tone to 2nd harmonic ratio
            mag_tmp = func_obj_1st_2nd_ratio(complex_sig_mag,Fs,fRange);
            ph_tmp = func_obj_1st_2nd_ratio(complex_sig_phase,Fs,fRange);
        elseif strcmp(optFcn, 'temp')
            % Objective Function 3 Template similarity
            mag_tmp = func_obj_Template(complex_sig_mag, Fs, fRange, opts);            
            ph_tmp = func_obj_Template(complex_sig_phase, Fs, fRange, opts);
        elseif strcmp(optFcn, 'cons')
            mag_tmp = func_obj_consistency(complex_sig_mag, Fs);            
            ph_tmp = func_obj_consistency(complex_sig_phase, Fs);
        elseif strcmp(optFcn, 'SNR_lin')
            mag_tmp = func_obj_SNRLin(complex_sig_mag,Fs,fRange,opts.nfRange);
            ph_tmp = func_obj_SNRLin(complex_sig_phase,Fs,fRange,opts.nfRange);
        elseif strcmp(optFcn, 'SIR')
            mag_tmp = func_obj_SIR(complex_sig_mag, Fs, fRange);
            ph_tmp = func_obj_SIR(complex_sig_phase, Fs, fRange);
        end
       
        % Generate vector injection scanning result map        
        mag_map(i_ph, i_mag) = mag_tmp;

        % Generate vector injection scanning result map 
        ph_map(i_ph, i_mag) = ph_tmp;

 %       disp(['Vector Injection Scanning Progress: ', ...
 %           num2str(((i_mag-1)*ph_step+i_ph)/(mag_step*ph_step)*100),'%'])

    end
end

    
end

function [score] = func_obj_SNRLin(ncs_sig, Fs, fRange, nfRange)
    % Version 1
    %snr_score = func_obj_fundamentalTone(sig, fs, fRange, nfRange);
    %lin_score = func_obj_1st_2nd_ratio(sig, fs, fRange);
    %score = snr_score + lin_score;

    % Version 2
    L_sig = length(ncs_sig);             % Length of signal
    Y_sig = fft(detrend(ncs_sig));
    P2 = abs(Y_sig/L_sig); L = round(L_sig/2);
    P1 = P2(1:round(L_sig/2)+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L_sig/2))/L_sig;
    % F
    mag_1_peak = max(P1(fRange(2) > f & f > fRange(1))); 
    %for i=1:4
    %    % Add all harmonic powers upto a specified number (4 here)
    %    mag_1_peak = mag_1_peak + max(P1(i*fRange(2) > f & f > i*fRange(1)));
    %end 

    amp_noise_freq  = median(P1(nfRange(2)>f&f>nfRange(1)));
    score = 10*log10(mag_1_peak/amp_noise_freq);
    score = score + 2*func_obj_1st_2nd_ratio(ncs_sig, Fs, fRange); 
    % hyper-parameter of relative linearity weighting can be tested 
end 

function [score] = func_obj_SIR(ncs_sig, Fs, fRange)
% Objective function 1 fundamental tone SNR
    L_sig = length(ncs_sig);             % Length of signal
    Y_sig = fft(detrend(ncs_sig));
    P2 = abs(Y_sig/L_sig); L = round(L_sig/2);
    P1 = P2(1:round(L_sig/2)+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L_sig/2))/L_sig;
    % F
    mag_1_peak = max(P1(fRange(2)>f&f>fRange(1)));
    amp_noise_freq  = median(P1(fRange(2)>f&f>fRange(1)));
    score = 10*log10(mag_1_peak/amp_noise_freq);

    % Add linearity as well
    score = score + func_obj_1st_2nd_ratio(ncs_sig, Fs, fRange);
end

function [SNR_sig_fundamental] = func_obj_fundamentalTone(ncs_sig,Fs,fRange,nfRange)
% Objective function 1 fundamental tone SNR
    L_sig = length(ncs_sig);             % Length of signal
    Y_sig = fft(detrend(ncs_sig));
    P2 = abs(Y_sig/L_sig); L = round(L_sig/2);
    P1 = P2(1:round(L_sig/2)+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L_sig/2))/L_sig;
    % F
    mag_1_peak = max(P1(fRange(2)>f&f>fRange(1)));

    amp_noise_freq  = median(P1(nfRange(2)>f&f>nfRange(1)));
    SNR_sig_fundamental = 10*log10(mag_1_peak/amp_noise_freq);

end

function inv_mse_signal = func_obj_consistency(ncs_sig, Fs)
    winLength = Fs;
    mse_all = [];
    for i = winLength:winLength:length(ncs_sig)-winLength
        mse_all(end+1) = sum(( zscore(ncs_sig(i+1:i+winLength)) - zscore(ncs_sig(i-winLength+1:i)) ).^2);
    end
    inv_mse_signal = 1/sum(mse_all);
end


function [SNR_sig_1st_2nd_ratio] = func_obj_1st_2nd_ratio(ncs_sig,Fs,fRange)
% Objective function 2 fundamental tone to 2nd harmonic ratio
    L_sig = length(ncs_sig);             % Length of signal
    Y_sig = fft(detrend(ncs_sig));
    P2 = abs(Y_sig/L_sig); 
    P1 = P2(1:round(L_sig/2)+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L_sig/2))/L_sig;
    % F
    mag_1_peak = max(P1(fRange(2)>f&f>fRange(1)));
    mag_2_peak = max(P1(fRange(2)*2>f&f>fRange(1)*2));
    SNR_sig_1st_2nd_ratio = 10*log10(mag_1_peak /mag_2_peak );

end

function cross_cor = func_obj_Template(ncs_sig, Fs, fRange, opts)
% Objective function 3 correlation to a waveform template
    templateComp = loadTemplate(opts.tempType, ncs_sig, Fs, fRange);
    
    if strcmp(opts.similarity, 'xcorr')
        corr_arr = xcorr(rescale(ncs_sig, -1, 1), templateComp', opts.winSecs*Fs, 'normalized');
    elseif strcmp(opts.similarity, 'sse')
        corr_arr = sse_win(rescale(ncs_sig, -1, 1), rescale(templateComp, -1, 1), opts.winSecs*Fs);
    end

    if opts.absolute
        cross_cor = max(abs(corr_arr));
    else
        cross_cor = max(corr_arr);
    end
end


