function [] = simpleInferenceWrapper(varargin)
% Model Options:
% 'full', 'sizeOnly', 'rewardOnly', 'dumbSize', 'dumbBias'
% '-sameAlpha'
% '-holdKappa'
% '-noLapse'
% '-useBias'

if nargin==0
    clear; clear global; clc;
    modelName = 'full';
    nRuns = 4;
    fillNewSubjects = 0;
else
    modelName = varargin{1};
    nRuns = 4;
    fillNewSubjects = 0;
end

if nargin>1
    nRuns = varargin{2};
end

if nargin>2
    fillNewSubjects = varargin{3};
end

subjectFilesOld = dir('../../Psychtoolbox/Results/Subjects/processed/*.mat');
nSubjectsOld = length(subjectFilesOld);
subjectFilesNew = dir('../../Psychtoolbox/JoystickResults/Subjects/*.mat');
subjectFiles = [subjectFilesOld; subjectFilesNew];
nSubjectsNew = length(subjectFilesNew);

nSubjects = length(subjectFiles);

zackData = ~cellfun('isempty', strfind({subjectFiles.name}', 'ZA'));

%% select model
% 'full', 'sizeOnly', 'rewardOnly', 'full-sameAlpha', 'dumbSize'

[filename, LB, UB, params] = getModelParams(modelName);
nParams = length(LB);

%%
if fillNewSubjects
    load('fits-dumbBias0.mat')
    nOldDone = sum(isnan(results.bestFits(:,2)));
    nNewDone = size(results.bestFits,1) - nOldDone;
    clear results
    nExtraOld = nSubjectsOld - nOldDone;
    nExtraNew = nSubjectsNew - nNewDone;
    key = zeros(nSubjects,1);
    key(nOldDone+1:nSubjectsOld) = 1;
    key(nNewDone+1:nSubjectsNew) = 1;
    %
    basefilename = filename(1:end-5);
    oldFiles = dir(['./' basefilename '0.mat']);
    filename = oldFiles.name; % we want to rewrite
    load(oldFiles.name)
    nAlreadyDone = size(results.allFits, 1);
    nRuns = size(results.allFits, 2); % rewrite so it can fit in the existing matrix
end

%% Preallocation!
nTrials = nan(nSubjects,1);
NLL = nan(nSubjects,1);
BIC = nan(nSubjects,1);
if params.bias
    bestFits = nan(nSubjects, nParams+7); % because of the extra for new data
else
    bestFits = nan(nSubjects, nParams+6);
end

if fillNewSubjects % fill in old data
    NLL(~key) = results.NLL;
    BIC(~key) = results.BIC;
    bestFits(~key, :) = results.bestFits;
    nTrials(~key) = results.nTrials;
    fits = nan(nSubjects, nRuns, size(results.bestFits,2));
    fits(~key, :, :) = results.allFits;
end

%% Run subject loop
tick = 1; % for making sure that the extra bias parameter is added once
for sid = 1:nSubjects
    if fillNewSubjects & ~key(sid)
        continue
    end
    
    if sid<=nSubjectsOld
        data_dir = '../../Psychtoolbox/Results/Subjects/processed/';
        load([data_dir subjectFiles(sid).name]);
        oldData = 1;
        rule = strcmp(data.ruleID, 'R1')+1;
    else
        data_dir = '../../Psychtoolbox/JoystickResults/Subjects/';
        load([data_dir subjectFiles(sid).name]);
        if sid==19
            data = data(1:500,:);
        end
        oldData = 0;
        rule = data.chosenInd;
        if params.bias && tick % only run this once!
            tick = 0;
            % need to add extra bias parameter by doubling the last bias bound
            LB = [LB LB(end)];
            UB = [UB UB(end)];
            nParams = length(LB);
        end
    end
    nTrials(sid) = size(data,1);
    
    % alpha, beta, kappa, lambda, lapse
    sid
    parfor k = 1:nRuns
        inits = LB + (UB-LB).*rand(1,nParams);
        opts = optimset('Display', 'none', 'MaxFunEvals', 30000, 'MaxIter', 15000, 'TolX', 1e-5);
        if strcmp(params.type, 'dumb')
            [x, fval, exitflag, output] = fmincon(@(coef) dumbInference(coef, data, oldData, params, 0), inits, [], [], [], [], LB, UB, [], opts);
        else
            [x, fval, exitflag, output] = fmincon(@(coef) simpleInference(coef, data, oldData, params, 0), inits, [], [], [], [], LB, UB, [], opts);
        end
        if params.bias && oldData % we need to add an extra spot so things fit nicely
            fits(sid,k,:) = [x NaN fval 0 0 0 output.iterations exitflag];
        else
            fits(sid,k,:) = [x fval 0 0 0 output.iterations exitflag];
        end
    end
    %%
    for k = 1:nRuns
        clear -global predictions
        global predictions
        if strcmp(params.type, 'dumb')
            dumbInference(fits(sid,k,:), data, oldData, params, 1);
        else
            simpleInference(fits(sid,k,:), data, oldData, params, 1);
        end
        accuracy = mean(rule==predictions);
        stayTrials = find(diff(rule)==0)+1;
        stayAcc = mean(rule(stayTrials)==predictions(stayTrials));
        switchTrials = find(diff(rule))+1;
        switchAcc = mean(rule(switchTrials)==predictions(switchTrials));
        fits(sid, k, end-4) = accuracy;
        fits(sid, k, end-3) = stayAcc;
        fits(sid, k, end-2) = switchAcc;
    end
    
    thisSubj = squeeze(fits(sid,:,:));
    bestRun = find(thisSubj(:,end-5) == min(thisSubj(:,end-5)));
    bestRun = bestRun(1);
    bestFits(sid,:) = squeeze(fits(sid,bestRun,:))';
    NLL(sid) = fits(sid,bestRun,end-5);
    BIC(sid) = 2*NLL(sid) + nParams*log(nTrials(sid)); % BIC is already calculated using the correct number of parameters when using bias
end

results.allFits = fits;
results.bestFits = bestFits;
results.nParams = nParams;
if params.bias % just to make it clear that it's different across groups
    results.nParams = nParams - 0.5;
end
results.nRuns = nRuns;
results.nTrials = nTrials;
results.NLL = NLL;
results.BIC = BIC;

save(filename, 'results');