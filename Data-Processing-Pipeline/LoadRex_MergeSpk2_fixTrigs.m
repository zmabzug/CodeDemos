function [REXData, varargout] =LoadRex_MergeSpk2_fixTrigs(rexName)

% Given the name of some Rex data files (the A and E files,
% without the ‘A’ or ‘E’ on the end), attempts a conversion of this data
% into a Matlab file that contains all spike, code, saccade, and other
% data.
% Loads data exported from Spike2 and adds/replaces spike times, unit IDs and waveforms
% Optional outputs: [data, Sp2Trigs] = ...
%                   [data, Sp2Trigs, nClus] =           

%% Get directories sorted out
currDir = pwd;
slashes = strfind(rexName, '/');
backslashes = strfind(rexName, '\');
lastFolderChar = max([slashes backslashes]);
folder = rexName(1:lastFolderChar);
filename = rexName(lastFolderChar+1:end);
cd(folder)

%% If LMA file...
ind = strfind(filename, '_A0');
LMA = 0;
Sp2filename = filename;
if ~isempty(ind)
    LMA = 1;
    filename = filename(1:ind-1);
end

%% Load Rex data
REXData = mrdr('-a','-d', filename);

%% get data from Spike2

if exist([folder Sp2filename 's.mat'])
    load([folder Sp2filename 's.mat']);
    load([folder filename 't.mat']); % use regular filename for triggers always
elseif exist([folder 'Spike2Exports/' Sp2filename 's.mat'])
    load([folder './Spike2Exports/' Sp2filename 's.mat'])
    load([folder './Spike2Exports/' filename 't.mat'])
else
    error('No Spike2 file found.')
end
    
% find which channel contains the Spike data, and which the Sync triggers
varlist=who; %list variables
varlist=varlist(~cellfun(@isempty,strfind(varlist,filename))); %restrict to the ones that start with the file name (the ones just loaded)
try
    chName = cellfun(@(x) eval([x '.title']), varlist,'UniformOutput',false);
catch
    for varN=1:size(varlist,1)
        chName{varN}=eval([varlist{varN} '.title']);
    end
end
chName(strcmp(chName, 'trig')) = cellstr('trigger'); % rename trigs to triggers
chName(strcmp(chName, '')) = cellstr('trigger'); % rename blanks to triggers
eval(['Spk2Data = ' cell2mat(varlist(cellfun(@isempty,strfind(chName,'trigger')))) ';']);
eval(['Spk2Trig = ' cell2mat(varlist(~cellfun(@isempty,strfind(chName,'trigger')))) ';']);
nClus = unique(Spk2Data.codes(:,1));

% find which triggers are the trial starts
SyncTrigs=int32(round((Spk2Trig.times-Spk2Trig.times(1))*1000));
% Extract triggers from REXData
allTrigs = [];
for n = 1:length(REXData)
    codes = [REXData(n).Events.Code];
    times = [REXData(n).Events.Time];
    
    % This klooge is not needed for files recorded after 9/1/2016
    % (previously, a second trigger code was not dropped on error trials
    % and the error code had to be used as a surrogate
%     trigsReal = times(codes==1502)';
%     trigsErr = times(codes==17385)' + 2;
%     trigs = sort([trigsReal; trigsErr]);

    trigs = times(codes==1502)';
    trigsPerTrial(n,1) = length(trigs);
    allTrigs = [allTrigs; trigs];
end
allTrigs = allTrigs - allTrigs(1);

% Compare differences in triggers!
minIter = min(length(allTrigs), length(SyncTrigs));
maxIter = max(length(allTrigs), length(SyncTrigs));
RexTrigs = double(allTrigs);
Sp2Trigs = double(SyncTrigs);
AdjSp2Trigs = double(SyncTrigs);
rexCounter = 1;
sp2Counter = 1;
for k = 1:maxIter;
    if rexCounter>size(RexTrigs)
        break
    end
    rawOffset(k) = RexTrigs(rexCounter) - AdjSp2Trigs(sp2Counter);
    if abs(rawOffset(k))<=7 % correct minor deviations in Adjusted only
        AdjSp2Trigs(sp2Counter:end) = AdjSp2Trigs(sp2Counter:end) + rawOffset(k);
        rexCounter = rexCounter + 1;
        sp2Counter = sp2Counter + 1;
    elseif abs(rawOffset(k))>5 % bad triggering
        if rawOffset(k)>0 % extra trigger(s) in Spike2
            AdjSp2Trigs(sp2Counter) = NaN;
            fprintf('Extra Spike2 trigger found: %i\n', sp2Counter)
            sp2Counter = sp2Counter + 1; % keep rexCounter to try and find a new match for it
        end
        if rawOffset(k)<0 % missing trigger(s) in Spike2 -> should never happen -> deliver error
            error('There are triggers missing from Spike2, which this code is not currently equipped to handle.')
        end
    end
end
if sp2Counter ~= rexCounter
    warning('%i extra triggers found and removed from %st.mat.\n', sp2Counter-rexCounter, filename)
end
newSp2Trigs = Sp2Trigs(~isnan(AdjSp2Trigs));
newAdjSp2Trigs = AdjSp2Trigs(~isnan(AdjSp2Trigs));

%plot(RexTrigs - newSp2Trigs)
SyncTrigs = newSp2Trigs; % rewriting!

% TrialTimes=sort([[REXData.aStartTime] [REXData.aEndTime]]);
TrialTimes=sort([REXData.tStartTime]);
TrialTimes=TrialTimes-TrialTimes(1);
% find trigger pattern (e.g., Start / Stop)
trigSeqStep=mode(diff(find(ismember(floor(SyncTrigs/10),floor(TrialTimes/10)))));
firstTrig=find(ismember(floor(SyncTrigs/10),floor(TrialTimes/10)),1);
startTrigs=Spk2Trig.times(firstTrig:trigSeqStep:end);
endTrigs=Spk2Trig.times(firstTrig+1:trigSeqStep:end);

% if length(startTrigs)<size(TrialTimes,2) %misssing some trials
%     maxDiv=min(diff(TrialTimes));
%     for trialSeek=1:round(maxDiv/50-1)
%         startTrigs=sort(unique([startTrigs; ...
%         find(ismember(floor(SyncTrigs/int32(50*trialSeek)),...
%         floor(TrialTimes/int32(50*trialSeek))))]));
%     end
% end

for  trialNumber = 1:min([length(startTrigs) REXData(end).trialNumber])
    REXData(trialNumber).Units=Spk2Data.codes(Spk2Data.times>=startTrigs(trialNumber) &...
        Spk2Data.times<=endTrigs(trialNumber),1);
    REXData(trialNumber).SpikeTimes=round((Spk2Data.times(Spk2Data.times>=startTrigs(trialNumber) &...
        Spk2Data.times<=endTrigs(trialNumber))-startTrigs(trialNumber))*1000);
    REXData(trialNumber).Waveforms=Spk2Data.values(Spk2Data.times>=startTrigs(trialNumber) &...
        Spk2Data.times<=endTrigs(trialNumber),:);
end

varargout{1} = SyncTrigs;
varargout{2} = nClus;

cd(currDir);
helpDebug = 0;

end

