clear; clear global; clc;

% define data directories
files = dir('fits-*.mat');
nModels = length(files);
names = {files.name}';
subjectFilesOld = dir('../../Psychtoolbox/Results/Subjects/processed/*.mat');
nSubjectsOld = length(subjectFilesOld);
subjectFilesNew = dir('../../Psychtoolbox/JoystickResults/Subjects/*.mat');
subjectFiles = [subjectFilesOld; subjectFilesNew];
nSubjectsNew = length(subjectFilesNew);

nSubjects = length(subjectFiles);

% load in models and relevant parameters
for fid = 1:nModels
    load(files(fid).name);
    allFits{fid} = results;
    BIC(fid,:) = results.BIC;
    NLL(fid,:) = results.NLL;
    nParams(fid) = results.nParams;
    maxPerf(fid,:) = results.bestFits(:,end-4)';
    clear results
end

%% extract biases to determine chance
for sid = 1:nSubjects
    if sid<=nSubjectsOld
        data_dir = '../../Psychtoolbox/Results/Subjects/processed/';
        load([data_dir subjectFiles(sid).name]);
        oldData = 1;
        r1pref = mean(strcmp(data.ruleID, 'R1'));
        ruleBias(sid) = max(r1pref, 1-r1pref);
    else
        data_dir = '../../Psychtoolbox/JoystickResults/Subjects/';
        load([data_dir subjectFiles(sid).name]);
        if sid==nSubjectsOld+1 % the first new subject
            data = data(1:500,:);
        end
        oldData = 0;
        ruleBias(sid) = max(histcounts(categorical(data.chosenInd), 'Normalization', 'probability'));
    end
    nTrials(sid) = size(data,1);
end
key = [zeros(1,nSubjectsOld) ones(1,nSubjectsNew)];
%%
% Remove any subjects with no data
BIC(:,all(isnan(BIC))) = [];
BICold = BIC(:, 1:nSubjectsOld);
BICnew = BIC(:, nSubjectsOld+1:end);
% Look at differences in BIC
dBIC = BIC - repmat(min(BIC), nModels, 1);

% do Bayesian model selection from Friston
[alpha,exp_r,xp,pxp,bor] = spm_BMS((-1/2).*BIC', [], [], [], [], 0.25*ones(1,nModels));
[alphaO,exp_rO,xpO,pxpO,borO] = spm_BMS((-1/2).*BICold', [], [], [], [], 0.25*ones(1,nModels));
[alphaN,exp_rN,xpN,pxpN,borN] = spm_BMS((-1/2).*BICnew', [], [], [], [], 0.25*ones(1,nModels));
[sortedBIC, ind] = sortrows([1-exp_r' dBIC]); % use 1-exp_r to get sorting in right direction
sortedBIC = sortedBIC(:, 2:end);
sortedNames = {files(ind).name}';
[indX5, indY5] = find(sortedBIC<5);
[indX0, indY0] = find(sortedBIC==0);

%% plot Bayesian model selection results
figure(1)
clf(1)
for n=1:2
    subplot(1,2,n)
    if n==1
        exp_plot = exp_rO;
        dBIC_plot = dBIC(:, 1:nSubjectsOld);
    else
        exp_plot = exp_rN;
        dBIC_plot = dBIC(:, nSubjectsOld+1:end);
    end
    [sortedBIC, ind] = sortrows([1-exp_plot' dBIC_plot]); % use 1-exp_r to get sorting in right direction
    sortedBIC = sortedBIC(:, 2:end);
    sortedNames = {files(ind).name}';
    [indX5, indY5] = find(sortedBIC<5);
    [indX0, indY0] = find(sortedBIC==0);
    h = imagesc(min(sortedBIC, 20));
    colormap(flipud(jet));
    if n==1
        cb = colorbar;
        cb.YTickLabel{end} = '>20';
    end
    shading flat;
    xlabel('Subject #')
    yyaxis left;
    axL = gca;
    axL.YTick = 1:nModels;
    axL.YTickLabel = strrep(sortedNames, '_', '\_');
    title('\DeltaBIC')
    axL.FontSize = 16;
    set(h,'AlphaData',~isnan(sortedBIC))
    
    yyaxis right;
    plot(indY5, nModels-indX5+0.5, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b')
    hold on
    plot(indY0, nModels-indX0+0.5, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g')
    hold off
    ylim([0 nModels])
    axR = gca;
    axR.YTick = 0.5:nModels;
    axR.YTickLabel = flipud(nParams(ind)');
    ylabel('nParams')
end
%% plot prediction improvement over chance
figure(2)
clf(2)
predImp = maxPerf-repmat(ruleBias,nModels,1);
sortedImp = sortrows([1-exp_r' predImp]); % sort predictions in the same way as dBIC
sortedImp = round(100*sortedImp(:, 2:end), 2);
predMax = nanmax(sortedImp(~strcmp(sortedNames, 'fits_humans_null.mat'),:));
predMaxMat = repmat(predMax, nModels, 1);
predFluc = sortedImp-predMaxMat;
predFluc(strcmp(sortedNames, 'fits_humans_null.mat'),:) = NaN;
predFluc(:, predMax==0) = NaN;
[indX5P, indY5P] = find(predFluc>-1);
[indX0P, indY0P] = find(predFluc==0);
h = imagesc(predFluc);
colormap(jet);
cb = colorbar;
shading flat;
xlabel('Subject #')
yyaxis left;
axL = gca;
axL.YTick = 1:nModels;
axL.YTickLabel = strrep(sortedNames, '_', '\_');
title('Prediction Improvement')
axL.FontSize = 16;
set(h,'AlphaData',~isnan(predFluc))

yyaxis right;
plot(indY5P, nModels-indX5P+0.5, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b')
hold on
plot(indY0P, nModels-indX0P+0.5, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g')
hold off
ylim([0 nModels])
axR = gca;
axR.YTick = 0.5:nModels;
axR.YTickLabel = flipud(nParams(ind)');
ylabel('nParams')

%% look at best model for each subject irrespective of BMS results
figure(3)
clf(3)
for n = 1:nSubjects
    b(n) = bar(n, predMax(n));
    hold on
    if n<=nSubjectsOld
        b(n).FaceColor = 'green';
    elseif n>nSubjectsOld
        b(n).FaceColor = 'cyan';
    end
end

ax = b(1).Parent;
ax.XTick = 1:nSubjects;
set(gca, 'FontSize', 16);
ylabel('Prediction Improvement (%)')
xlabel('Subject Number')
xlim([0.4 nSubjects+0.6])
%% best models for different subsets of subjects
figure(4)
clf(4)

newBIC = BIC(:, 1:nSubjectsOld);
joyBIC = BIC(:, nSubjectsOld+1:end);
[anew,exp_new,newxp,~,~] = spm_BMS((-1/2).*newBIC', [], [], [], [], 0.25*ones(1,nModels));
[ajoy,exp_joy,joyxp,~,~] = spm_BMS((-1/2).*joyBIC', [], [], [], [], 0.25*ones(1,nModels));
allxp = [newxp; joyxp];
allxp = [exp_new; exp_joy];
[~,c] = find(allxp>0.05);
modelsOfInterest = unique(c);
allxpRel = allxp(:, modelsOfInterest);
namesRel = names(modelsOfInterest);
barh([1-sum(allxpRel,2) allxpRel])
set(gca, 'YTickLabel', {'Old-Easy Subjects'; 'Joystick'})
set(gca, 'FontSize', 16);
xlabel('Model probability for a random subject')
legend(['All others'; strrep(namesRel, '_', '\_')], 'Location', 'best')
title('Bayesian Model Selection')

% try high performers only
[ahigh,exp_high,highxp,~,~] = spm_BMS((-1/2).*BIC(:,predMax>20 & 1:nSubjects>nSubjectsOld)', [], [], [], [], 0.25*ones(1,nModels));
%% pick best levels for each factor
myNames = names(3:end);
typeKey = ~cellfun('isempty', strfind(myNames, 'full')) + ... % full
    ~cellfun('isempty', strfind(myNames, '-sameAlpha')) + ... % full-sameAlpha
    ~cellfun('isempty', strfind(myNames, 'Only'))*3 + ... % rewOnly
    ~cellfun('isempty', strfind(myNames, 'size')); % sizeOnly
lapseKey = cellfun('isempty', strfind(myNames, 'noLapse'));
kappaKey = cellfun('isempty', strfind(myNames, 'holdKappa'));
biasKey = ~cellfun('isempty', strfind(myNames, 'useBias'));

% type
thisAlpha = anew(3:end);
alphaType = [sum(thisAlpha(typeKey==1)) sum(thisAlpha(typeKey==2)) sum(thisAlpha(typeKey==3)) sum(thisAlpha(typeKey==4))];
xpType = spm_dirichlet_exceedance(alphaType,1e6);

alphaLapse = [sum(thisAlpha(lapseKey==0)) sum(thisAlpha(lapseKey==1))];
xpLapse = spm_dirichlet_exceedance(alphaLapse,1e6);

alphaKappa = [sum(thisAlpha(kappaKey==0)) sum(thisAlpha(kappaKey==1))];
xpKappa = spm_dirichlet_exceedance(alphaKappa,1e6);

alphaBias = [sum(thisAlpha(biasKey==0)) sum(thisAlpha(biasKey==1))];
xpBias = spm_dirichlet_exceedance(alphaBias,1e6);

%% comparing lapse rates and model fits
figure(5)
clf(5)
biased = ruleBias>.5;

colors = repmat([0 0 0], nSubjects, 1);
colors(:,1) = biased';

lossPerTrial = NLL./repmat(nTrials, nModels, 1);
[lpt, best] = min(lossPerTrial, [], 1);
subplot(1,3,1)
h = scatter(allFits{17}.bestFits(nSubjectsOld+1:end,5), predMax(nSubjectsOld+1:end), 100, colors(nSubjectsOld+1:end,:), 'filled');
hold on
plot([0 max(h.XData)*1.1], [20 20], 'b--')
hold off
xlim([0 max(h.XData)*1.1])
ylabel('Prediction Improvement for Best Model (%)')
xlabel('Lapse parameter from full model')
set(gca, 'FontSize', 16);
pbaspect([1 1 1]);

subplot(1,3,2)
h = scatter(lpt(nSubjectsOld+1:end), predMax(nSubjectsOld+1:end), 100, colors(nSubjectsOld+1:end,:), 'filled');
hold on
plot([0 max(h.XData)*1.1], [20 20], 'b--')
hold off
xlim([0 max(h.XData)*1.1])
xlabel('NLL normalized by # of trials')
set(gca, 'FontSize', 16);
pbaspect([1 1 1]);

subplot(1,3,3)
varFits = std(NLL,[],1).^2;
h = scatter(log10(varFits(nSubjectsOld+1:end)), predMax(nSubjectsOld+1:end), 100, colors(nSubjectsOld+1:end,:), 'filled');
hold on
plot([0 max(h.XData)*1.1], [20 20], 'b--')
hold off
xlim([0 max(h.XData)*1.1])
xlim([0 6])
xlabel('log_{10}(Variance) of NLL across Models')
set(gca, 'FontSize', 16);
pbaspect([1 1 1]);