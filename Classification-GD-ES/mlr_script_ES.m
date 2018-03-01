% multinomial logistic regression
clear;

%whosedata = 'mine'; % or 'paul'
whosedata = 'paul';
type = 'normal'; % or 'multinomial'
validation = 'leave-one-out'; % or 'thirds'
%validation = 'thirds';

if strcmp(whosedata, 'mine')
    cd('~/Desktop/zackdata/Rasters/Sac1');
    files = dir('S*');
    %files = files(1);
elseif strcmp(whosedata, 'paul')
    cd('~/Desktop/pauldata/zackRasters');
    load('SEF_alignrefix_highvslow.mat')
    files = allrasts;
else break
end

for i=1:length(files) %make sure length(files) gives right number
    
    % Clear variables
    clearvars -except files i oops oddstab resultsT resultsR results whosedata type validation
    
    % Splitpoint for loading data
    if strcmp(whosedata, 'mine')
        eval(['load ' files(i).name]);
        allrastmat = [rasts{1}; rasts{2}; rasts{3}; rasts{4}];
    elseif strcmp(whosedata, 'paul')
        allrastmat = [files{i,1}; files{i,2}];
    end
    
    % Set parameters
    binsize = 100;
    nT = size(allrastmat,1); % # of trials
    nL = size(allrastmat,2); % # of timepoints
    
    % Sort rasters into bins
    nbins = floor(nL/binsize);
    for b = 1:nbins;
        binned(:,b) = sum(allrastmat(:,(b-1)*binsize+1:b*binsize),2);
    end
    
    % Remove NaN trials from data
    int = isnan(binned(:,1));
    full = binned(~int,:);
    nT2 = size(full,1);
    
    % Splitpoint for generating truths
    if strcmp(whosedata, 'mine')
        % Make keys (truths) for trialtype and rule
        if strcmp(type, 'multinomial') % 1:SSR0, 2:SSR1, 3:INR0, 4:INR1
            key = [1*ones(size(rasts{1},1), 1); 2*ones(size(rasts{2},1), 1); 3*ones(size(rasts{3},1), 1); 4*ones(size(rasts{4},1), 1)];
            key = key(~int,:);
        else % normal binary (logistic) comparisons
            keyTT = 1+[zeros(size(rasts{1},1), 1); zeros(size(rasts{2},1), 1); ones(size(rasts{3},1), 1); ones(size(rasts{4},1), 1)];
            keyRR = 1+[zeros(size(rasts{1},1), 1); ones(size(rasts{2},1), 1); zeros(size(rasts{3},1), 1); ones(size(rasts{4},1), 1)];
            keyTT = keyTT(~int,:);
            keyRR = keyRR(~int,:);
            key = [keyTT keyRR];
        end
    elseif strcmp(whosedata, 'paul')
        % Make keys for bet
        key = 1 + [zeros(size(files{i,1},1), 1); ones(size(files{i,2},1), 1)];
        key = key(~int,:);
    end

    % Splitpoint for running ES
    if strcmp(whosedata, 'mine')
        if strcmp(type, 'multinomial') % 1:SSR0, 2:SSR1, 3:INR0, 4:INR1
            % DO NOTHING RIGHT NOW
        else % normal binary (logistic) comparisons
            pINS = sum(keyTT-1)/length(keyTT);
            pSS = 1 - pINS;
            pR1 = sum(keyRR-1)/length(keyRR);
            pR0 = 1 - pR1;
            % Run ES on trial type
            if strcmp(validation, 'leave-one-out')
                [prob_valid pred_valid Ptrain Ptest] = logitES(full, keyTT-1);
                resultsT(i,:) = [mean(pred_valid) mean(Ptrain) mean(Ptest)];
                % Use binomial test to see if prediction accuracy is better than chance
                [phatT pciT] = binofit(sum(pred_valid), length(pred_valid), 0.05);
                hT = pciT(1) > (max(pINS,pSS));
            elseif strcmp(validation, 'thirds')
                [perf_validT PtrainT PtestT] = logitESthirds(full, keyTT-1);
                [hT pT] = ttest(perf_validT, max(pINS,pSS), 'tail', 'right');
            end
            if hT
                oops(i,1) = 1;
            else oops(i,1) = 0;
            end
            % Run ES on rules
            if strcmp(validation, 'leave-one-out')            
                [prob_valid pred_valid Ptrain Ptest] = logitES(full, keyRR-1);
                resultsR(i,:) = [mean(pred_valid) mean(Ptrain) mean(Ptest)];
                % Use binomial test to see if prediction accuracy is better than chance    
                [phatR pciR] = binofit(sum(pred_valid), length(pred_valid), 0.05);
                hR = pciR(1) > (max(pR1,pR0));
            elseif strcmp(validation, 'thirds')
                [perf_validR PtrainR PtestR] = logitESthirds(full, keyRR-1);
                [hR pR] = ttest(perf_validR, max(pR1,pR0), 'tail', 'right');
            end
            if hR
                oops(i,2) = 1;
            else oops(i,2) = 0;
            end
        end
    elseif strcmp(whosedata, 'paul')
        pL = sum(key-1)/length(key);
        pH = 1 - pL;
        % Run ES on bets
        if strcmp(validation, 'leave-one-out')
            [prob_valid pred_valid Ptrain Ptest] = logitES(full, key-1);
            results(i,:) = [mean(pred_valid) mean(Ptrain) mean(Ptest)];
            % Use binomial test to see if prediction accuracy is better than chance
            [phat pci] = binofit(sum(pred_valid), length(pred_valid), 0.05);
            h = pci(1) > (max(pL,pH));
        elseif strcmp(validation, 'thirds')
            [perf_valid Ptrain Ptest] = logitESthirds(full, key-1);
            [h p] = ttest(perf_valid, max(pL,pH), 'tail', 'right');   
        end
        if h
            oops(i,1) = 1;
        else oops(i,1) = 0;
        end
    end
    
%     idxT = sub2ind(size(pihatT), (1:length(pihatT))', keyTT);
%     idxR = sub2ind(size(pihatR), (1:length(pihatR))', keyRR);
%     
%     % Do log-odds test
%     oddsT = pihatT(idxT)./((sum(keyTT-1).*(keyTT-1) + sum(2-keyTT).*(2-keyTT) - 1)/(length(keyTT)-1));
%     oddsR = pihatR(idxR)./((sum(keyRR-1).*(keyRR-1) + sum(2-keyRR).*(2-keyRR) - 1)/(length(keyRR)-1));
%     logoddsT = log2(oddsT);
%     logoddsR = log2(oddsR);
%     eloT = mean(logoddsT);
%     eloR = mean(logoddsR);
%     pT = 1/(1+2^-eloT);
%     pR = 1/(1+2^-eloR);
%     oddstab(i,:) = [mean(oddsT) eloT pT oops(i,1) mean(oddsR) eloR pR oops(i,2)];
    
    [i sum(oops,1)]
end

