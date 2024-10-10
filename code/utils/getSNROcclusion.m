function score = getSNROcclusion(ncs_sig, Fs, fRange)
    L_sig = length(ncs_sig);             % Length of signal
    Y_sig = fft(detrend(ncs_sig));
    P2 = abs(Y_sig/L_sig); 
    P1 = P2(1:round(L_sig/2)+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L_sig/2))/L_sig;
    % F
    mag_1_peak = max(P1(fRange(2)>f&f>fRange(1)));
    % maximum interference is from motion rather than heart. 
    interf1 = max(P1(fRange(1)>f&f>0));
    % interf1 = median(P1(3>f&f>0.5));
    
    score = mag_1_peak/interf1;

    % This makes it mathemtically unstable 
    if 0 
    for i=2:4
        mag_pk = max(P1(fRange(2)*i>f&f>fRange(1)*i));
        %inter_mn = median(P1(fRange(2)*i>f&f>fRange(1)*i));
        inter_mn = interf1; 
        score = score + 10*log10(mag_pk/inter_mn);
    end 
    score = score/4; % alternative method of computation
    end 
end