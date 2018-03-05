clear; clear global; clc;

nRuns = 8;
modelList = {
    'full-noLapse';
    'sizeOnly-noLapse';
    'rewardOnly-noLapse';
    'full-sameAlpha-noLapse';
    };
% modelList = {'dumbBias'};

fillNewSubjects = 1;

runAllModels = 1;
if runAllModels
    m1 = {'full'; 'sizeOnly'; 'rewardOnly'; 'full-sameAlpha'};
    m2 = {[]; '-holdKappa'};
    m3 = {[]; '-noLapse'};
    m4 = {[]; '-useBias'};
    for m = 1:length(m1)*length(m2)*length(m3)*length(m4)
        m1ind = ceil(m./(length(m2)*length(m3)*length(m4)));
        m2ind = mod(ceil(m./1),2)+1;
        m3ind = mod(ceil(m./2),2)+1;
        m4ind = mod(ceil(m./4),2)+1;
        modelList{m} = [m1{m1ind} m2{m2ind} m3{m3ind} m4{m4ind}];
    end
    modelList{m+1} = 'dumbSize';
    modelList{m+2} = 'dumbBias';
end

nModelsToRun = length(modelList);

for mid = 1:nModelsToRun
    try
        simpleInferenceWrapper(modelList{mid}, nRuns, fillNewSubjects)
    catch
        keyboard
    end
end