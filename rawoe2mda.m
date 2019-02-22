function [] = rawoe2mda(directory)

if nargin == 0
    directory = pwd;
end

cd(directory);

rawFiles = dir('*CH*.continuous');
chansPerTrode = length(rawFiles);
allData = cell(chansPerTrode,1);

for ii=1:chansPerTrode
   [data,timestamps] = load_open_ephys_data_faster(rawFiles(ii).name);
   allData{ii} = data;
end

totalTime = max(timestamps)-min(timestamps)+10/1000;
startTime = min(timestamps);

fullData = zeros(length(timestamps),chansPerTrode);

for ii=1:chansPerTrode
   fullData(:,ii) = allData{ii}; 
end
clear allData;

writemda(fullData',sprintf('raw.full.mda'),'float64');

Fs = round(1/(timestamps(2)-timestamps(1)));
save('RecordingInfo.mat','totalTime','startTime','Fs','timestamps','chansPerTrode');
end