clear; clc;

% select a data file to look at
filename = 'M004L3P1_0601';

% sent trigger times: from REX, stored in allTrigs
% received trigger times: from Spike2, stored in Spk2Trigs
% when things are working correctly, one trigger is sent at the start of each trial, and one trigger is sent at the end of each trial

% extract trigger times from Spike2, extract behavioral data from REX
[data, Spk2Trigs] = LoadRex_MergeSpk2(filename);

% or use this function instead to verify that phantom triggers were removed from Spike2
%[data, Spk2Trigs] = LoadRex_MergeSpk2_fixTrigs(filename);

%% Extract trigger times from REX data
% M007/20160828a and earlier: there is no second 1502 on error trials
allTrigs = [];
for n = 1:length(data)
    codes = [data(n).Events.Code];
    times = [data(n).Events.Time];
    % these are the "codes" that correspond to times when triggers are sent
    trigsReal = times(codes==1502)';
    trigsErr = times(codes==17385)' + 2; % klooge
    trigs = sort([trigsReal; trigsErr]);
    trigsPerTrial(n,1) = length(trigs);
    allTrigs = [allTrigs; trigs];
end

% normalize
allTrigs = double(allTrigs - allTrigs(1));

%% Compare to Spike2 triggers
nTrigsREX = length(allTrigs);
nTrigsSp2 = length(Spk2Trigs);
% in case there aren't the same number of triggers, make vectors same length for easy plotting
minTrigs = min(nTrigsREX, nTrigsSp2);
allTrigsPlot = allTrigs(1:minTrigs)./1000;
Spk2TrigsPlot = double(Spk2Trigs(1:minTrigs))./1000;

% plot trigger times against each other, and unity line for comparison
subplot(1,3,1)
plot(Spk2TrigsPlot, allTrigsPlot, 'r.', 'MarkerSize', 10)
set(gca,'FontSize',16)
xlabel('Spike2 trigger times (s)')
ylabel('REX trigger times (s)')
title(sprintf('REX triggers: %i, Spk2 triggers: %i', nTrigsREX, nTrigsSp2));
hold on
plot(Spk2TrigsPlot, Spk2TrigsPlot, 'k--'); % unity line
hold off

% plot histograms of deviations between sent and received triggers
subplot(2,3,2)
histogram(allTrigsPlot-Spk2TrigsPlot)
xlabel('Deviation in trigger times (REX-Sp2) (s)')
ylabel('Count')
set(gca,'FontSize',16)

% plot deviations as a function of sent trigger times, to look at how deviation changes over time
subplot(2,3,5)
plot(allTrigsPlot, allTrigsPlot-Spk2TrigsPlot)
ylabel('Deviation in trigger times (REX-Sp2) (s)')
xlabel('REX trigger times')
set(gca,'FontSize',16)

% look at histogram of gaps between received start-of-trial triggers
diffSpk2 = diff(Spk2TrigsPlot);
subplot(2,3,3)
histogram(diffSpk2(1:2:end));
xlabel('During trial')
set(gca,'FontSize',16)
% look at histogram of gaps between received end-of-trial triggers
subplot(2,3,6)
histogram(diffSpk2(2:2:end));
title('Gap between Spk2 triggers')
xlabel('Between trial')
set(gca,'FontSize',16)

