% useage of this is the following:
%   run the rest of the mountainsort pipeline
%   then, you should have a directory filled with -mounsort.mat files
%   go into that directory and run this code
%   the output will be a single file with all of the data from the
%   experiment, including the LFP ... the LFP obviously adds a ton of data
%   to the file, so if you have a ton of neurons, you may have to remove
%   the LFP or think about storing the data individually

files = dir('*-mounsort.mat');

allUnits = 0;
for ii=1:length(files)
   load(files(ii).name,'totalUnits');
   allUnits = allUnits+totalUnits;
end

clear totalUnits;

allData = struct('fileName',cell(allUnits,1),'recordDate',cell(allUnits,1),...
    'recordTime',cell(allUnits,1),'expDay',cell(allUnits,1),'SpikeWaveForm',cell(allUnits,1),...
    'spikeTimes',cell(allUnits,1),'meanFR',cell(allUnits,1),'eventTimes',cell(allUnits,1),...
    'eventIDs',cell(allUnits,1),'LFP',cell(allUnits,1),'trodalness',...
    cell(allUnits,1),'gender',cell(allUnits,1),'species',cell(allUnits,1));

preAmpGain = 1;

count = 0;
for ii=1:length(files)
   load(files(ii).name); 
   for jj=1:totalUnits
      count = count+1;
      allData(count).fileName = files(ii).name(1:end-13);
      allData(count).recordDate = DateTime;
      allData(count).recordTime = totalTime;
      
      index = regexp(files(ii).name,'ay');
      
      allData(count).expDay = str2double(files(ii).name(index+2));
      
      tempwave = cell(3,1);
      for kk=1:3
          tempwave{kk} = meanwaves{jj,kk};
      end
      
      allData(count).SpikeWaveForm = tempwave;
      allData(count).spikeTimes = allts{jj};
      allData(count).meanFR = length(allts{jj})/totalTime;
      allData(count).eventTimes = timeStamps;
      allData(count).eventIDs = svStrobed;
      
      Chans = find(~cellfun(@isempty,allad));numChans = length(Chans);
      if numChans==0
          
      elseif length(allad{Chans(1)}) < 10
          
      else
          tempData = allad{jj};
          voltage = 1000.*((tempData).*SlowPeakV)./(0.5*...
              (2^SlowADResBits)*adgains(jj)*preAmpGain);
          
          n = 2;
          notch = 60/(adfreqs(jj)/2);
          bw = notch/n;
          [b,a] = iirnotch(notch,bw);
          allData(count).LFP = filtfilt(b,a,voltage);
      end
      
      allData(count).trodalness = chansPerTrode;
      allData(count).gender = 'Male';
      allData(count).species = 'C57BL/6';
   end
end

save('RachelLiz_CompileUnitData.mat','allData');