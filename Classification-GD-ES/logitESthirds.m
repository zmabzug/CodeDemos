function [perf_valid Ptrain Ptest] = logitESthirds(binrast, truth)
% Implementing early stopping for logistic regression. Input binned spike
% counts and truths for trial conditions.

[nTrials nBins] = size(binrast);
numreps = 20;

for n = 1:numreps;    
    % Initialize weights to zero and add constant to bins
    gma = 0.03;
    W = randn(nBins+1,1)';
    smallbin = [ones(nTrials,1) binrast];
    max_iter = 2000;
    
    % Randomly remove 1/3 of trials for validation
    idx = randperm(nTrials);
    thirdTrials = floor(nTrials/3);
    valid = smallbin(idx(1:thirdTrials), :);
    test = smallbin(idx(thirdTrials+1:2*thirdTrials), :);
    train = smallbin(idx(2*thirdTrials+1:end), :);

    validkey = truth(idx(1:thirdTrials));
    testkey = truth(idx(thirdTrials+1:2*thirdTrials));
    trainkey = truth(idx(2*thirdTrials+1:end));

    % Iterate
    for iter = 1:max_iter
        % For test data: calculate log-prob
        wdx_test = test*W';
        lp_test = 1+log2(exp(1))*mean(testkey.*wdx_test - log(1 + exp(wdx_test))); % mutual information (in bits)

        % For training data: calculate log-prob, gradient, double gradient
        wdx = train*W';
        pygx = 1./(1+exp(-wdx));
        lp_train = 1+log2(exp(1))*mean(trainkey.*wdx - log(1 + exp(wdx)));
        Dlp_train = (trainkey - pygx)'*train;
        DDlp_train = -train'*(repmat(pygx.*(1-pygx), 1, nBins+1).*train);

        % Update weights
        Wtemp = W - (gma.*pinv(DDlp_train)*Dlp_train')';

        % For test data: calculate new log-prob
        wdx_test = test*Wtemp';
        lp_test_NEW = 1+log2(exp(1))*mean(testkey.*wdx_test - log(1 + exp(wdx_test)));
        
        %Output log-probabilities and %correct for train & test
        logprobtrain(iter) = lp_train;
        logprobtest(iter) = lp_test;
        prob_train = exp(trainkey.*(train*W'))./(1 + exp(train*W')); % P(Y=y|x)
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
        
        % End if max_iter is reached
        if iter == max_iter
            disp('Max iter reached')
        end
    end

    % Validation on left-out trial
    prob_valid = exp(validkey.*(valid*W'))./(1 + exp(valid*W')); % P(Y=y|x)
    pred_valid = prob_valid>0.5;
    perf_valid(n) = mean(pred_valid);
    
    % Also output final performance on testing and training
    Ptrain(n) = perf_train(end);
    Ptest(n) = perf_test(end);
    
end

end
