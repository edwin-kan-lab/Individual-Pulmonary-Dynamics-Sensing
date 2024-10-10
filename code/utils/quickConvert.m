% quickConvert.m
dPath = ['E:\OneDrive - Cornell University\NIH Pig Study\Pig Study 03-07-24\Data Files\'];
d = dir(dPath);


%%
for i = 1:length(d)
    if contains(d(i).name, 'NCS') && contains(d(i).name, '.tdms')
        NCSFile = d(i).name;
        BIOFile = strrep(NCSFile, 'NCS', 'BIO');
        RadioFile = strrep(NCSFile, 'NCS', 'Radio');
        
        try
            saveTDMStoMAT_NCSpigV5(dPath, strrep(NCSFile, '.tdms', ''));
            saveTDMStoMAT_BIOpig(dPath, strrep(BIOFile, '.tdms', ''));
            saveTDMStoMAT_RadiopigV5(dPath, strrep(RadioFile, '.tdms', ''));
        catch 
            disp(['Unable to convert ', NCSFile, '. The file may be corrupted or have zero size']); 
        end
    end
end