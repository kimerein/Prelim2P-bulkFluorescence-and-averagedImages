function [output,p_runningAfterOpto]=plotBehaviorEffects(fileDir)

% Set variables
encoderName='Wheel_Encoder'; % string of file name for analog encoder data
optoName='Opto_Stim'; % string of file name for voltage output to laser data
nFramesBeforeStim=4;
alignTraces=1;

% Read in running data
wheel_trials=readPhysFiles(fileDir,encoderName);

% Read in laser command data
laser_trials=readPhysFiles(fileDir,optoName);

% Plot trial-by-trial effects
trials_isRunning=findRunningFrames(wheel_trials,122);

% Sort trials by whether any running prior to opto stim
times=linspace(0,31.232,122);
isRunningAfterOpto=nan(1,length(trials_isRunning));
timeWindow=[1 31]; % in s
preStim_run_trials=zeros(1,length(trials_isRunning));
for i=1:length(trials_isRunning)
    currTrial=trials_isRunning{i};
    % If any trials are running prior to opto stim, then consider mouse to
    % be running before opto stim
    if any(currTrial(1:nFramesBeforeStim)>0)
        preStim_run_trials(i)=1;
    else
        preStim_run_trials(i)=0;
    end
    if any(currTrial(times>=timeWindow(1) & times<=timeWindow(2))>0)
        isRunningAfterOpto(i)=1;
    else
        isRunningAfterOpto(i)=0;
    end
end
p_runningAfterOpto=nansum(isRunningAfterOpto)./length(trials_isRunning);

% Plot position across trials according to whether any running prior to opto stim
% Running before opto
preRun_data=nan(length(wheel_trials),length(wheel_trials{1}.data));
preNonrun_data=nan(length(wheel_trials),length(wheel_trials{1}.data));
preRun_laser=nan(length(laser_trials),length(laser_trials{1}.data));
preNonrun_laser=nan(length(laser_trials),length(laser_trials{1}.data));
for i=1:length(wheel_trials)
    if preStim_run_trials(i)==1
        preRun_data(i,:)=wheel_trials{i}.data;
        preRun_laser(i,:)=laser_trials{i}.data;
    else
        preNonrun_data(i,:)=wheel_trials{i}.data;
        preNonrun_laser(i,:)=laser_trials{i}.data;
    end
end
% Align position data before opto stim
for i=1:size(preRun_data,1)
    preRun_data(i,:)=preRun_data(i,:)-nanmean(preRun_data(i,1:nFramesBeforeStim));
end
for i=1:size(preNonrun_data,1)
    preNonrun_data(i,:)=preNonrun_data(i,:)-nanmean(preNonrun_data(i,1:nFramesBeforeStim));
end
% Plot aligned data
times=linspace(0,31.232,length(wheel_trials{1}.data));
figure();
subplot(2,1,1);
plot(times,preRun_data','r');
subplot(2,1,2);
plot(times,preRun_laser','r');
title('Running Before Opto Stim');
figure();
subplot(2,1,1);
plot(times,preNonrun_data','b');
subplot(2,1,2);
plot(times,preNonrun_laser','b');
title('Non-Running Before Opto Stim');
output.preRun_data=preRun_data;
output.preRun_laser=preRun_laser;
output.preNonrun_data=preNonrun_data;
output.preNonrun_laser=preNonrun_laser;

% Calculate probability of running after opto stim
