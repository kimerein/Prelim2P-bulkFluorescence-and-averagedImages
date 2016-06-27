function [mode_firstShuttered,mode_secondShuttered,baselineFrames,postOptoFrames,preShutterFrames,postShutterFrames]=findShutteredTrials(trials,thresh,nFramesAfterOpto,nFramesForBase)

% trials: cell array of cell arrays (frames) -- each frame in each trial
%       contains tif stack of images collected in that trial
% Suggestion: thresh=4*10^7; 
% thresh: Threshold below which frame is considered to be blank
% nFramesAfterOpto: Number of frames to average after opto stim or "fake" opto stim (i.e., just shutter)
% nFramesForBase: Number of frames preceding opto stim to average for baseline image

% Find stim (real opto) and stim (fake with just shutter) frames
trialByTrial_shuttered=cell(1,length(trials));
nShutteredTrials=1;
firstShuttered=nan(1,length(trials));
secondShuttered=nan(1,length(trials));
for i=1:length(trials)
    trial_frames=trials{i};
    sum_of_trial=nan(length(trial_frames),1);
    for j=1:length(trial_frames)
        sum_of_trial(j)=sum(sum(trial_frames{j},1),2);
    end
    shuttered_trials=find(sum_of_trial<thresh);
    if length(shuttered_trials)>2
        if i==1
            disp('num shuttered trials');
            disp(floor(length(shuttered_trials)/2));
        end
        nShutteredTrials=ceil(length(shuttered_trials)/2);
        shuttered_trials=shuttered_trials([1 end]);
    end
    trialByTrial_shuttered{i}=shuttered_trials;
    firstShuttered(i)=shuttered_trials(1);
    secondShuttered(i)=shuttered_trials(2);
end
mode_firstShuttered=mode(firstShuttered);
mode_secondShuttered=mode(secondShuttered);

% Assume that first shuttered trial is opto stim and second shuttered trial
% is just shutter (without opto stim)
% Find frame indices
baselineFrames=cell(1,length(trials));
postOptoFrames=cell(1,length(trials));
preShutterFrames=cell(1,length(trials));
postShutterFrames=cell(1,length(trials));
for i=1:length(trials)
    trial_frames=trials{i};
    shuttered_trials=trialByTrial_shuttered{i};
    if nShutteredTrials>2
        shuttered_trials=[mode_firstShuttered mode_secondShuttered];
    end
    % These are the frames before opto stim
    baselineFrames{i}=shuttered_trials(1)-1-nFramesForBase:shuttered_trials(1)-2;
    % These are the frames after opto stim
    postOptoFrames{i}=shuttered_trials(1)+nShutteredTrials-1+1:shuttered_trials(1)+nShutteredTrials-1+nFramesAfterOpto;
    % These are the frames before shutter (no opto stim)
    preShutterFrames{i}=shuttered_trials(2)-(nShutteredTrials-1)-1-nFramesForBase:shuttered_trials(2)-(nShutteredTrials-1)-2;
    % These are the frames after shutter (no opto stim)
    postShutterFrames{i}=shuttered_trials(2)-(nShutteredTrials-1)-1-nFramesForBase:shuttered_trials(2)-(nShutteredTrials-1)-2;
end