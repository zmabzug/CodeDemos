clear; clc;

filename = 'M004L3P1_0601';
%filename = 'M005L3P1_0621';
%filename = 'M006L4A0_1481';

[data, Spk2Trigs] = LoadRex_MergeSpk2_fixTrigs(filename);
%[data, Spk2Trigs] = LoadRex_MergeSpk2(filename);

%% Extract trigger times from REX data
% M007/20160828a and earlier: there is no second 1502 on error trials
allTrigs = [];
for n = 1:length(data)
    codes = [data(n).Events.Code];
    times = [data(n).Events.Time];
    trigsReal = times(codes==1502)';
    trigsErr = times(codes==17385)' + 2; % because
    trigs = sort([trigsReal; trigsErr]);
    trigsPerTrial(n,1) = length(trigs);
    allTrigs = [allTrigs; trigs];
end

allTrigs = double(allTrigs - allTrigs(1));

%% Compare to Spike2 triggers
nTrigsREX = length(allTrigs);
nTrigsSp2 = length(Spk2Trigs);
% in case there aren't the same number of triggers
minTrigs = min(nTrigsREX, nTrigsSp2);
allTrigsPlot = allTrigs(1:minTrigs)./1000;
Spk2TrigsPlot = double(Spk2Trigs(1:minTrigs))./1000;

subplot(1,3,1)
plot(Spk2TrigsPlot, allTrigsPlot, 'r.', 'MarkerSize', 10)
set(gca,'FontSize',16)
xlabel('Spike2 trigger times (s)')
ylabel('REX trigger times (s)')
title(sprintf('REX triggers: %i, Spk2 triggers: %i', nTrigsREX, nTrigsSp2));
hold on
plot(Spk2TrigsPlot, Spk2TrigsPlot, 'k--'); % unity line
hold off

subplot(2,3,2)
histogram(allTrigsPlot-Spk2TrigsPlot)
xlabel('Deviation in trigger times (REX-Sp2) (s)')
ylabel('Count')
set(gca,'FontSize',16)

subplot(2,3,5)
plot(allTrigsPlot, allTrigsPlot-Spk2TrigsPlot)
ylabel('Deviation in trigger times (REX-Sp2) (s)')
xlabel('REX trigger times')
set(gca,'FontSize',16)

diffSpk2 = diff(Spk2TrigsPlot);
subplot(2,3,3)
histogram(diffSpk2(1:2:end));
xlabel('During trial')
set(gca,'FontSize',16)
subplot(2,3,6)
histogram(diffSpk2(2:2:end));
title('Gap between Spk2 triggers')
xlabel('Between trial')
set(gca,'FontSize',16)

% %% Investigate further using trials processed by rex_process
% if exist([filename '.mat'], 'file')
%    load([filename '.mat'])
%    trigCodesPerTrial = sum(allcodes==1502 | allcodes==17385, 2);
%    % Find if there are any trials (other than the last one) with ~=2
%    weirdTrials = find(trigCodesPerTrial(1:end-1)~=2);    
% end
