function [] = rex_process_Sp2(human, monkey, varargin)
% rex_process_Sp2(human, monkey, filename[, filename, filename..., autoRewrite])
% e.g. rex_process_Sp2('Zack', 'Muzzy', 'M004L3P1_0601', 'M005L3P1_0621', 1)
% 
% rex_process_Sp2 will check to see if a file has already been processed,
% and will ask if you'd like to re-write any already-processed data. Add in
% autoRewrite = 1 to automatically re-write any processed data files
% without asking for confirmation [autoRewrite = 0 by default].

% here, we can accomodate different root directories (depending on the user and their computer)
if strcmp(human, 'Zack')
    root = '~/Desktop/zackdata/';
end

% and different data folders, depending on the animal the data came from
if strcmp(monkey, 'Muzzy')
    folder = 'Muzzy/';
end

% check and see if the user wants to automatically rewrite files (if necessary), or if they want to be prompted each time
autoRewrite = 0;
if isa(varargin{end}, 'double')
    if varargin{end}==1
        autoRewrite = 1;
    end
    % remove autoRewrite tag from varargin
    varargin = varargin(1:end-1);
end

% step through each file that the user wants to process
for n=1:length(varargin);
    filename = varargin{n};
    
    %% Check if file exists yet
    processedFolder = [root 'processed' filesep folder];
    newfilename = [processedFolder filename '_Sp2.mat'];
    rexname = [filename '_Sp2'];
    
    % if necessary, prompt the user to ask if they want to rewrite files
    if exist(newfilename) & ~autoRewrite
        prompt = [rexname ' already exists in ' processedFolder '. Re-write data file (y/n)?\n'];
        response = input(prompt, 's');
        if strcmpi(response, 'n')
            fprintf('%s has been skipped.\n', filename)
            continue
        end
    end
            
    %% Load behavioral and neuro data, and fix trigger timing if needed!
    [data, ~, clusters] = LoadRex_MergeSpk2_fixTrigs([root folder filename]);
    
    %% Convert raw behavioral data (data.Events) to matrices of trial events (allcodes) and timing (alltimes)
    allcodes = [];
    alltimes = [];
    allh = [];
    allv = [];
    spikes = cell(length(clusters),1);
    
    % step through each trial
    for t=1:length(data)
        rawCodes = double([data(t).Events.Code]); % Extract codes
        theseCodes = rawCodes(rawCodes>=1000); % Only take real codes (no 1, 2, 800, 801, etc.)
        rawTimes = double([data(t).Events.Time]);
        theseTimes = rawTimes(rawCodes>=1000);
        
        % make times relative to start trigger (1502) which comes 1 state (~= 1 ms) before start of trial (1001)
        theseTimes = theseTimes - theseTimes(theseCodes==1001) + 1;
        
        if t==1 & theseCodes(1)==1502 % get rid of first trigger for uniformity
            theseCodes = theseCodes(2:end);
            theseTimes = theseTimes(2:end);
        end
        
        % save data into matrices
        allcodes = cat_variable_size_row(allcodes, theseCodes);
        alltimes = cat_variable_size_row(alltimes, theseTimes);
        allh = cat_variable_size_row(allh, data(t).Signals(1).Signal);
        allv = cat_variable_size_row(allv, data(t).Signals(2).Signal);
        
        % extract neuro data from multiple neurons
        for c = 1:length(clusters);
            oldSpikes = spikes{c};
            newSpikes = data(t).SpikeTimes(data(t).Units==clusters(c))';
            if isempty(newSpikes); % we never want to add [] because that won't add a new row
                newSpikes = NaN;
            end
            allSpikes = cat_variable_size_row(oldSpikes, newSpikes, 1);
            spikes{c} = allSpikes;            
        end
    end    
    allspk_clus = spikes;
    allspk = spikes{1};
    rexnumtrials = length(data);
    
    %% Save variables
    save(newfilename, 'allcodes', 'alltimes', 'allspk', 'allspk_clus', 'allh', 'allv', 'rexnumtrials', 'rexname');
    
    clearvars -EXCEPT root folder varargin n autoRewrite
    
    forDebug = 0;
end
