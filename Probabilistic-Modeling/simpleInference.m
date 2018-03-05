function [err] = simpleInference(coef, data, oldData, params, globalOutputs)
%% deciper params
c = 1;
useSize = 0;
useRew = 0;

% determine which parameters to fit based on params structure
if ~isempty(strfind(params.type, 'both')) % double alpha and kappa if using both
    useSize = 1;
    useRew = 1;
    
    alpha = coef(c); c=c+1;
    if params.sameAlpha % only use one alpha if forcing to be the same (aka alpha and beta are same)
        beta = coef(c-1);
    else
        beta = coef(c); c=c+1;
    end
    if params.kappaFree
        kappa = coef(c); c=c+1;
        lambda = coef(c); c=c+1;
    else
        kappa = 1;
        lambda = 1;
    end
elseif strcmp(params.type, 'size')
    useSize = 1;
    alpha = coef(c); c=c+1;
    
    if params.kappaFree
        kappa = coef(c); c=c+1;
    else
        kappa = 1;
    end
elseif strcmp(params.type, 'reward')
    useRew = 1;
    beta = coef(c); c=c+1;
    if params.kappaFree
        lambda = coef(c); c=c+1;
    else
        lambda = 1;
    end
end

if params.lapse
    lapse = coef(c); c=c+1;
else
    lapse = 1e-9; % to try to prevent zero probabilities
end

%%%% always last assignment %%%%
if params.bias
    bias1 = coef(c); c=c+1;
    if ~oldData
        bias2 = coef(c); c=c+1;
    end
end

%% ensure compatibility across data sources
if oldData % for compatibility with data from old task
    nRules = 2;
    rule = strcmp(data.ruleID, 'R1')+1;
    % assume size=red, brightness=green, and ignore blue
else
    nRules = 3;
    rule = data.chosenInd;
end
nTrials = size(data,1);
r = [0:nTrials-1];
zeroTrials = zeros(nTrials,1);

%% calculate probabilities based only on size information
if useSize
    if oldData
        redSize = abs(data.sizeLeft - data.sizeRight);
        greenSize = abs(data.contLeft - data.contRight);
        scale = 72.3; % determined
        greenSize = greenSize*scale;
    else
        redSize = data.redWidth;
        greenSize = data.greenWidth;
        blueSize = data.blueWidth;
    end
    
    % implement learning and smoothing
    alphaVec = (alpha.*ones(length(nTrials))).^r;
    redSmoothed = tsmovavg([zeroTrials; redSize]', 'w', fliplr(alphaVec));
    redSmoothed = redSmoothed(nTrials:end-1) + redSize';
    greenSmoothed = tsmovavg([zeroTrials; greenSize]', 'w', fliplr(alphaVec));
    greenSmoothed = greenSmoothed(nTrials:end-1) + greenSize';
    
    % calculate running estimate of mean
    redMean = redSmoothed./cumsum(alphaVec);
    greenMean = greenSmoothed./cumsum(alphaVec);
    
    if ~oldData
        blueSmoothed = tsmovavg([zeroTrials; blueSize]', 'w', fliplr(alphaVec));
        blueSmoothed = blueSmoothed(nTrials:end-1) + blueSize';
        blueMean = blueSmoothed./cumsum(alphaVec);
    end
    
    % assume constant variance but changing mean
    redVar = cumsum((redSize'-redMean).^2)./(r+1);
    greenVar = cumsum((greenSize'-greenMean).^2)./(r+1);
    if ~oldData
        blueVar = cumsum((blueSize'-blueMean).^2)./(r+1);
    end
    
    % calculate probabilities that one choice is more beneficial than others
    % p(r>g & r>b) = p(r>g)*p(g>b) + p(r>b)*p(b>g)
    if oldData
        rVSg = 1-normcdf(0, redMean-greenMean, sqrt(redVar + greenVar));
        pRed = rVSg;
        pGreen = 1-rVSg;
        allProbs = [pRed' pGreen'];
    else
        rVSg = 1-normcdf(0, redMean-greenMean, sqrt(redVar + greenVar));
        rVSb = 1-normcdf(0, redMean-blueMean, sqrt(redVar + blueVar));
        gVSb = 1-normcdf(0, greenMean-blueMean, sqrt(greenVar + blueVar));
        pRed = rVSg.*gVSb + rVSb.*(1-gVSb);
        pGreen = (1-rVSg).*rVSb + gVSb.*(1-rVSb);
        pBlue = (1-rVSb).*rVSg + (1-gVSb).*(1-rVSg);
        allProbs = [pRed' pGreen' pBlue'];
    end
    
    % control explore/exploit tradeoff (kappa)
    allProbs = allProbs./repmat(sum(allProbs,2),1,nRules);
    allProbs(allProbs==0) = min(allProbs(allProbs>0)); % get rid of any potential zero probabilities
    allProbs = allProbs.^kappa;
    allProbs(isinf(allProbs)) = realmax/1e2; % get rid of any infinites
    sizeProbs = allProbs./repmat(sum(allProbs,2),1,nRules);
    clear allProbs
else
    sizeProbs = ones(nTrials, nRules);
end
%% calculate probabilities based only on reward
if useRew
    betaVec = (beta.*ones(length(nTrials))).^r;
    
    if oldData
        redAlpha = strcmp(data.ruleID, 'R0').*(data.reward==1);
        redBeta = strcmp(data.ruleID, 'R0').*(data.reward==0);
        greenAlpha = strcmp(data.ruleID, 'R1').*(data.reward==1);
        greenBeta = strcmp(data.ruleID, 'R1').*(data.reward==0);
    else
        redAlpha = strcmp(data.chosenTarget, 'red').*(data.perf==1);
        redBeta = strcmp(data.chosenTarget, 'red').*(data.perf==0);
        greenAlpha = strcmp(data.chosenTarget, 'green').*(data.perf==1);
        greenBeta = strcmp(data.chosenTarget, 'green').*(data.perf==0);
        blueAlpha = strcmp(data.chosenTarget, 'blue').*(data.perf==1);
        blueBeta = strcmp(data.chosenTarget, 'blue').*(data.perf==0);
    end
    
    % use beta distributions to estimate dynamics reward rates
    % learning and smoothing
    redAlphaSmoothed = tsmovavg([zeroTrials; redAlpha]', 'w', fliplr(betaVec));
    redAlphaSmoothed = redAlphaSmoothed(nTrials:end-1) + redAlpha' + beta;
    redBetaSmoothed = tsmovavg([zeroTrials; redBeta]', 'w', fliplr(betaVec));
    redBetaSmoothed = redBetaSmoothed(nTrials:end-1) + redBeta' + beta;
    redEP = exp(psi(redAlphaSmoothed) - psi(redAlphaSmoothed+redBetaSmoothed));
    
    greenAlphaSmoothed = tsmovavg([zeroTrials; greenAlpha]', 'w', fliplr(betaVec));
    greenAlphaSmoothed = greenAlphaSmoothed(nTrials:end-1) + greenAlpha' + beta;
    greenBetaSmoothed = tsmovavg([zeroTrials; greenBeta]', 'w', fliplr(betaVec));
    greenBetaSmoothed = greenBetaSmoothed(nTrials:end-1) + greenBeta' + beta;
    greenEP = exp(psi(greenAlphaSmoothed) - psi(greenAlphaSmoothed+greenBetaSmoothed));
    
    if ~oldData
        blueAlphaSmoothed = tsmovavg([zeroTrials; blueAlpha]', 'w', fliplr(betaVec));
        blueAlphaSmoothed = blueAlphaSmoothed(nTrials:end-1) + blueAlpha' + beta;
        blueBetaSmoothed = tsmovavg([zeroTrials; blueBeta]', 'w', fliplr(betaVec));
        blueBetaSmoothed = blueBetaSmoothed(nTrials:end-1) + blueBeta' + beta;
        blueEP = exp(psi(blueAlphaSmoothed) - psi(blueAlphaSmoothed+blueBetaSmoothed));
        allProbs = [redEP' greenEP' blueEP'];
    else
        allProbs = [redEP' greenEP'];
    end
    
    % control explore/exploit tradeoff (lambda)    
    allProbs = allProbs./repmat(sum(allProbs,2),1,nRules);
    allProbs(allProbs==0) = min(allProbs(allProbs>0)); % get rid of any potential zero probabilities
    allProbs = allProbs.^lambda;
    allProbs(isinf(allProbs)) = realmax/1e2; % get rid of any infinites
    rewardProbs = allProbs./repmat(sum(allProbs,2),1,nRules);
    clear allProbs
else
    rewardProbs = ones(nTrials, nRules);
end

%% determine bias
if params.bias
    if oldData
        biasProbVec = [bias1 1-bias1];
    else
        biasProbVec = [bias1 bias2 max(1-bias1-bias2, 0)];
        biasProbVec = biasProbVec./sum(biasProbVec);
    end
    biasProbs = repmat(biasProbVec, nTrials, 1);
else
    biasProbs = ones(nTrials, nRules);
end

%% integrate visual + outcome + bias + lapse
allProbs = sizeProbs.*rewardProbs.*biasProbs;
allProbs = [ones(1,nRules) ; allProbs(1:end-1, :)];
allProbs = allProbs./repmat(sum(allProbs,2),1,nRules);
allProbs = (1-lapse).*allProbs + lapse.*ones(nTrials,nRules)/nRules;

% calculate NLL to be minimized
pCr = diag(allProbs(:,rule));
err = -sum(log(pCr));

if isnan(err) || isinf(err)
    err = 1e9;
end

[~,preds] = max(allProbs, [], 2);

% output variables if necessary
if globalOutputs==1
    global predictions
    predictions = preds;
elseif globalOutputs==2
    global sizeP
    global rewP
    sizeP = sizeProbs;
    rewP = rewardProbs;
end
