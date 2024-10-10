% Run all experiments
clc; 
addpath(genpath(pwd))
dates = {'02-05-24', '02-06-24', '02-07-24', '02-08-24', '03-04-24', ...
    '03-05-24', '03-06-24', '03-07-24'}; 
toSave = true; 

for runAllInd = 1:length(dates)
    clearvars -except dates files* runAllInd toSave
    dFileDate = dates{runAllInd}; 
    
    files = {'tidalvol', 'lungocclusion'}; 

    for filesInc = 1:length(files)
        clearvars -except dates files* calibs scorings runAllInd ...
            dFileDate preload_ncs toSave
        interventionType = files{filesInc};
        disp(['Run all processing ' dFileDate ' ' interventionType]);
        
        study_main
    end
end

% Ensure everything is closed 
close all;
clear all; 