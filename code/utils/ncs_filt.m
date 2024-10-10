 function ncsAmp_Filt = ncs_filt(iqAmp,f_low,f_high,Fs,filtType)

   %% Filtering    
    orderHP = 20; % 20 original
    
    if strcmp(filtType,'lowpass') 
        % Nonlinear phase filter but better filtering characteristics
        filtLP = fdesign.lowpass('N,F3db',orderHP,f_high,Fs); 
        % kaiserwin is much faster than equiripple without much degradation (?)
        HdLP = design(filtLP,'butter'); 
        ncsAmp_Filt = filtfilt(HdLP.sosMatrix,HdLP.ScaleValues,iqAmp);
    elseif strcmp(filtType,'highpass')
        % Nonlinear phase filter but better filtering characteristics
        filtHP = fdesign.highpass('N,F3db',orderHP,f_low,Fs); 
        % 'butter' with 'N,F3db' specifications 
        HdHP = design(filtHP,'butter'); 
        ncsAmp_Filt = filtfilt(HdHP.sosMatrix,HdHP.ScaleValues,iqAmp);      
    elseif strcmp(filtType,'bandpass')
        % Nonlinear phase filter but better filtering characteristics
        filtHP = fdesign.highpass('N,F3db',orderHP,f_low,Fs); 
        % 'butter' with 'N,F3db' specifications 
        HdHP = design(filtHP,'butter'); 
        % Nonlinear phase filter but better filtering characteristics
        filtLP = fdesign.lowpass('N,F3db',orderHP,f_high,Fs); 
        % kaiserwin is much faster than equiripple without much degradation (?)
        HdLP = design(filtLP,'butter'); 
        ncsAmp_Filt = filtfilt(HdHP.sosMatrix,HdHP.ScaleValues,iqAmp);
        ncsAmp_Filt = filtfilt(HdLP.sosMatrix,HdLP.ScaleValues,ncsAmp_Filt);  
    elseif strcmp(filtType, 'bandstop')
        % Nonlinear phase filter but better filtering characteristics
        filtHP = fdesign.highpass('N,F3db',orderHP,f_high,Fs); 
        % 'butter' with 'N,F3db' specifications 
        HdHP = design(filtHP,'butter'); 
        % Nonlinear phase filter but better filtering characteristics
        filtLP = fdesign.lowpass('N,F3db',orderHP,f_low,Fs); 
        % kaiserwin is much faster than equiripple without much degradation (?)
        HdLP = design(filtLP,'butter'); 
        ncsAmp_Filt = filtfilt(HdHP.sosMatrix,HdHP.ScaleValues,iqAmp);
        ncsAmp_Filt2 = filtfilt(HdLP.sosMatrix,HdLP.ScaleValues,ncsAmp_Filt);  
        ncsAmp_Filt = ncsAmp_Filt2 + ncsAmp_Filt;
    else
        ncsAmp_Filt = detrend(iqAmp);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end