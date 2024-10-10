% Originally Pragya Sharma, 12 July 2018
% Edited by Thomas Conroy, 01 January 2020: 4 channel with syncronized
% BIOPAC data
% tbc38@cornell.edu

function [convertedData,NcsData,subjectInfo] = saveTDMStoMAT_NCSpigV5(dataPath, fileName)

% SAVETDMSTOMAT specifies data path, and file name to be converted from
% TDMS to mat, and saveFile name is the name of the saved mat file
% It uses convertTDMS.m to convert TDMS to mat.

% dataPath = 'C:\Users\tbc38\Documents\Pig Study\TrialPigData\PigData\';
% fileName = 'CaseNCSSaving0308_102942';
saveFileName =  fileName;

[convertedData, NcsData, subjectInfo] = callConvert(dataPath,fileName);
% [convertedData, NcsData, subjectInfo] = callConvert(dataPath,fileName);

NcsData1 = NcsData(:,1:8);
NcsData2 = NcsData(:,9:16);

saveFileName = [dataPath,saveFileName];%,'.mat'];

save(saveFileName, 'NcsData1', 'NcsData2', 'subjectInfo');
% save(saveFileName, 'NcsData', 'subjectInfo');

end

function [convertedData, NcsData, subjectInfo] = callConvert(dataPath,fileName)
% CALLCONVERT takes data path and file name as input. It calls 

filePathName = [dataPath,fileName,'.tdms'];
[convertedData,~,~,~,~] = convertTDMS(false,filePathName);

NcsData = [ convertedData.Data.MeasuredData(4).Data,...
            convertedData.Data.MeasuredData(5).Data,...
            convertedData.Data.MeasuredData(6).Data,...
            convertedData.Data.MeasuredData(7).Data,...
            convertedData.Data.MeasuredData(8).Data, ...
            convertedData.Data.MeasuredData(9).Data,...
            convertedData.Data.MeasuredData(10).Data,...
            convertedData.Data.MeasuredData(11).Data,...
            convertedData.Data.MeasuredData(12).Data,...
            convertedData.Data.MeasuredData(13).Data,...
            convertedData.Data.MeasuredData(14).Data,...
            convertedData.Data.MeasuredData(15).Data,...
            convertedData.Data.MeasuredData(16).Data, ...
            convertedData.Data.MeasuredData(17).Data,...
            convertedData.Data.MeasuredData(18).Data,...
            convertedData.Data.MeasuredData(19).Data,...
            ];

subjectInfo = convertedData.Data.Root.Property; 
end

