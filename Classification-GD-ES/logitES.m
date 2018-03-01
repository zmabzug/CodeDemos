function [prob_valid pred_valid Ptrain Ptest] = logitES(binrast, truth)
% Implementing early stopping for logistic regression. Input binned spike
% counts and truths for trial conditions.


% binrast =binrast - ones(size(binrast,1),1)*mean(binrast);
% CC=cov(binrast);
% binrast = binrast*pinv(sqrtm(CC));

[nTrials nBins] = size(binrast);
begin = 1; % usually
%begin = nTrials; % do last trial for testing

for n = begin:nTrials;
    % Leave one out
    index = boolean(ones(nTrials,1));
    index(n) = false;
    smallbin = binrast(index,:);
    smallkey = truth(index,:); % key should be 0s and 1s
    
    % Initialize weights to zero and add constant to bins
    gma = 0.03;
    W = zeros(nBins+1,1)'; 
    smallbin = [ones(nTrials-1,1) smallbin];
    max_iter = 2000;
    
    % Sort remaining trials into training and testing
    idx = randperm(nTrials-1);
%     idx1 = find(smallkey==1);
%     idx0 = find(smallkey==0);
%     idx1 = permute(idx1,randperm(length(idx1)));
%     idx0 = permute(idx0,randperm(length(idx0)));
%     halfTrials1 = ceil(length(idx1)/2);
%     halfTrials0 = ceil(length(idx0)/2);
% 
%     train = smallbin(idx0(1:halfTrials0),:);
%     train = [train;smallbin(idx1(1:halfTrials1),:);];
%     trainkey = [zeros(halfTrials0,1);ones(halfTrials1,1)];
%     
%     test = smallbin(idx0(halfTrials0+1:end),:);
%     test = [test;smallbin(idx1(halfTrials1+1:end),:);];
%     testkey = [zeros(length(idx0)-halfTrials0,1);ones(length(idx1)-halfTrials1,1)];
    
    halfTrials = ceil(nTrials/2);
    train = smallbin(idx(1:halfTrials),:);
    trainkey = smallkey(idx(1:halfTrials));
    test = smallbin(idx(halfTrials+1:end),:);
    testkey = smallkey(idx(halfTrials+1:end));

%    pinvCstatic = inv(diag(diag(-train'*train/4))); % non-conjugate
%    (constant) gradient

    % Iterate
    for iter = 1:max_iter;
        % For test data: calculate log-prob
        wdx_test = test*W';
%        lp_test = sum(testkey.*wdx_test - log(1 + exp(wdx_test))); 
        lp_test = 1+log2(exp(1))*mean(testkey.*wdx_test - log(1 + exp(wdx_test))); % mutual information (in bits)

        % For training data: calculate log-prob, gradient, double gradient
        wdx = train*W';
        pygx = 1./(1+exp(-wdx));
        lp_train = 1+log2(exp(1))*mean(trainkey.*wdx - log(1 + exp(wdx)));
        Dlp_train = (trainkey - pygx)'*train;
        DDlp_train = -train'*(repmat(pygx.*(1-pygx), 1, nBins+1).*train);

        % Update weights
        Wtemp = W - (gma.*pinv(DDlp_train)*Dlp_train')';
%        Wtemp = W - (gma.*pinvCstatic*Dlp_train')'; % non-conjugate
%        gradient

        % For test data: calculate new log-prob
        wdx_test = test*Wtemp';
        lp_test_NEW = 1+log2(exp(1))*mean(testkey.*wdx_test - log(1 + exp(wdx_test)));
        
        %Output log-probabilities and %correct for train & test
        logprobtrain(iter) = lp_train;
        logprobtest(iter) = lp_test;
        prob_train = exp(trainkey.*(train*W'))./(1 + exp(train*W')); % P(Y=y|x)
        % prob_train = exp(1.*(train*W'))./(1 + exp(train*W')); % P(Y=1|x)
        prob_test = exp(testkey.*(test*W'))./(1 + exp(test*W'));
        pred_train = round(prob_train);
        pred_test = round(prob_test);
        perf_train(iter) = mean(pred_train);
        perf_test(iter) = mean(pred_test);
        
        % Check for exit criteria
        if lp_test_NEW < lp_test
            break
        else W = Wtemp;
        end   
%         W = Wtemp;
        
        % End if max_iter is reached
        if iter == max_iter;
            disp('Max iter reached')
        end
    end
    
    % Plotting diagnostics (only for one left-out trial)
    if begin==nTrials;
        subplot(1, 2, 1);
        plot(1:50, logprobtrain(1:50), 'k-o', 1:50, logprobtest(1:50), 'r-o');
        xlim([0.5 50.5]);
        legend('Training', 'Testing');
        xlabel('Iteration #');
        ylabel('Mutual information (bits)');

        subplot(1, 2, 2);
        plot(1:50, perf_train(1:50), 'k-o', 1:50, perf_test(1:50), 'r-o');
        xlim([0.5 50.5]);
        legend('Training', 'Testing');
        xlabel('Iteration #');
        ylabel('Classification Performance');
        ylim([0 1]);
    end

    % Validation on left-out trial
    prob_valid(n) = exp(truth(n).*([1 binrast(n,:)]*W'))./(1 + exp([1 binrast(n,:)]*W')); % P(Y=y|x)
    pred_valid(n) = prob_valid(n)>0.5;
    
    % Also output final performance on testing and training
    Ptrain(n) = perf_train(end);
    Ptest(n) = perf_test(end);
    
end
