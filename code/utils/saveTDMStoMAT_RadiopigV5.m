% Originally Pragya Sharma, 12 July 2018
% Edited by Thomas Conroy, 01 January 2020: 4 channel with syncronized
% BIOPAC data
% tbc38@cornell.edu

function [convertedData,RadioData,subjectInfo] = saveTDMStoMAT_RadiopigV5(dataPath, fileName)

% SAVETDMSTOMAT specifies data path, and file name to be converted from
% TDMS to mat, and saveFile name is the name of the saved mat file
% It uses convertTDMS.m to convert TDMS to mat.

% dataPath = 'C:\Users\tbc38\Documents\Pig Study\TrialPigData\PigData\';
% fileName = 'CaseBIOSaving0308_102942';
saveFileName =  fileName;

[convertedData, RadioData, subjectInfo] = callConvert(dataPath,fileName);
% [convertedData, NcsData, subjectInfo] = callConvert(dataPath,fileName);

saveFileName = [dataPath,saveFileName];%,'.mat'];

save(saveFileName, 'RadioData', 'subjectInfo');
% save(saveFileName, 'NcsData', 'subjectInfo');

end

function [convertedData, Data, subjectInfo] = callConvert(dataPath,fileName)
% CALLCONVERT takes data path and file name as input. It calls 

filePathName = [dataPath,fileName,'.tdms'];
[convertedData,~,~,~,~] = convertTDMS(false,filePathName);

Data = [ convertedData.Data.MeasuredData(4).Data,...
            ];
Data = [Data(1:12:end) ...
        Data(2:12:end) ...
        Data(3:12:end) ...
        Data(4:12:end) ...
        Data(5:12:end) ...
        Data(6:12:end) ...
        Data(7:12:end) ...
        Data(8:12:end) ...
        Data(9:12:end) ...
        Data(10:12:end) ...
        Data(11:12:end) ...
        Data(12:12:end) ...
        ];
subjectInfo = convertedData.Data.Root.Property; 
end

