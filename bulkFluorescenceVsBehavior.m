function bulkFluorescenceVsBehavior(fileDir,saveDir,saveName)

% Set variables
trialDuration=31.232;

% Get bulk fluorescence
fluorByTrial=getBulkFluorescenceOverTime(fileDir);

% Get behavioral (running) data
encoderName='Wheel_Encoder'; % string of file name for analog encoder data
% Read in running data
wheel_trials=readPhysFiles(fileDir,encoderName);

% Get opto stim data
optoName='Opto_Stim'; % string of file name for voltage output to laser data
% Read in laser command data
laser_trials=readPhysFiles(fileDir,optoName);

% Get shutter data
shutterName='Opto_Coming'; % string of file name for voltage output to laser data
% Read in laser command data
shutter_trials=readPhysFiles(fileDir,shutterName);

% Exclude shuttered trials from bulk fluorescence
% First, find shuttered trials
isShuttered=nan(length(shutter_trials),length(shutter_trials{1}.data));
for i=1:length(shutter_trials)
    curr=shutter_trials{i};
    data=curr.data;
    isShuttered(i,:)=data>1;
end
consensusShutter=mode(isShuttered,1);
% Second, nan these shuttered trials in bulk fluorescence
c=consensusShutter(1:floor(length(consensusShutter)/size(fluorByTrial,2)):end);
% Throw out trials when shutter may be mechanically opening or closing
c_firstHalf=c(1:floor(length(c)/3));
c_secondHalf=c(floor(length(c)/3)+1:end);
f=find(c_firstHalf>0.5,1,'first');
c(f-1)=1;
f=find(c_firstHalf>0.5,1,'last');
c(f+1)=1;
f=find(c_secondHalf>0.5,1,'first');
c(floor(length(c)/3)+f-1)=1;
f=find(c_secondHalf>0.5,1,'last');
% c(floor(length(c)/3)+f+1:floor(length(c)/3)+f+2)=1;
c(floor(length(c)/3)+f+1)=1;
c(1)=1;
fluorByTrial(:,c>0.5)=nan;

% Find running frames
isRunning=findRunningFrames(wheel_trials,length(c));
% Find trials in which mouse was running prior to opto stim
f=find(c_firstHalf>0.5,1,'first');
preOptoRunning=any(isRunning(:,1:f-1),2);

% Compare bulk fluorescence to running data
% Plot example trials
plotNtrials=15;
exampleTrials=randi(size(fluorByTrial,1),1,plotNtrials);
for i=1:length(exampleTrials)
    figure();
    subplot(3,1,1);
    times=linspace(0,trialDuration,size(fluorByTrial,2));
    plot(times,fluorByTrial(exampleTrials(i),:));
    axis tight;
    xlim([0 max(times)]);
    title('Bulk Fluorescence - Example Trial');
    subplot(3,1,2);
    times=linspace(0,trialDuration,length(wheel_trials{exampleTrials(i)}.data));
    plot(times,wheel_trials{exampleTrials(i)}.data);
    axis tight;
    xlim([0 max(times)]);
    title('Position - Example Trial');
    subplot(3,1,3);
    times=linspace(0,trialDuration,length(laser_trials{exampleTrials(i)}.data));
    plot(times,laser_trials{exampleTrials(i)}.data./100,'Color','c');
    hold on;
    plot(times,shutter_trials{exampleTrials(i)}.data,'Color','k');
    axis tight;
    xlim([0 max(times)]);
    title('Laser: Blue, Shutter: Black - Example Trial');
end

wheel_data=nan(length(wheel_trials),length(wheel_trials{1}.data));
for i=1:length(wheel_trials)
    wheel_data(i,:)=wheel_trials{i}.data;
end
laser_data=nan(length(laser_trials),length(laser_trials{1}.data));
for i=1:length(laser_trials)
    laser_data(i,:)=laser_trials{i}.data;
end
shutter_data=nan(length(shutter_trials),length(shutter_trials{1}.data));
for i=1:length(shutter_trials)
    shutter_data(i,:)=shutter_trials{i}.data;
end

% Plot average across all trials
figure();
subplot(3,1,1);
times=linspace(0,trialDuration,size(fluorByTrial,2));
fill([times fliplr(times)],[nanmean(fluorByTrial,1)+nanstd(fluorByTrial,[],1)./sqrt(size(fluorByTrial,1)) fliplr(nanmean(fluorByTrial,1)-nanstd(fluorByTrial,[],1)./sqrt(size(fluorByTrial,1)))],'c');
hold on;
plot(times,nanmean(fluorByTrial,1));
axis tight;
xlim([0 max(times)]);
title('Bulk Fluorescence - Average');
subplot(3,1,2);
times=linspace(0,trialDuration,length(wheel_trials{exampleTrials(1)}.data));
plot(times,nanmean(wheel_data,1));
axis tight;
xlim([0 max(times)]);
title('Position - Average');
subplot(3,1,3);
times=linspace(0,trialDuration,length(laser_trials{exampleTrials(1)}.data));
plot(times,nanmean(laser_data,1)./100,'Color','c');
hold on;
plot(times,nanmean(shutter_data,1),'Color','k');
axis tight;
xlim([0 max(times)]);
title('Laser: Blue, Shutter: Black - Average');

% Plot average across all trials separated by running or non-running before
% opto stim
tit='Running Before Opto';
plotData(trialDuration,fluorByTrial,1,wheel_data,laser_data,shutter_data,tit,preOptoRunning);
tit='Stationary Before Opto';
plotData(trialDuration,fluorByTrial,0,wheel_data,laser_data,shutter_data,tit,preOptoRunning);

% Save output
runBefore=1;
output.fluorByTrial_preRun=fluorByTrial(preOptoRunning==runBefore,:);
output.wheel_data_preRun=wheel_data(preOptoRunning==runBefore,:);
output.laser_data_preRun=laser_data(preOptoRunning==runBefore,:);
output.shutter_data_preRun=shutter_data(preOptoRunning==runBefore,:);
runBefore=0;
output.fluorByTrial_preNonrun=fluorByTrial(preOptoRunning==runBefore,:);
output.wheel_data_preNonrun=wheel_data(preOptoRunning==runBefore,:);
output.laser_data_preNonrun=laser_data(preOptoRunning==runBefore,:);
output.shutter_data_preNonrun=shutter_data(preOptoRunning==runBefore,:);
save([saveDir '\' saveName '.mat'],'output');

% Correlations

end

function plotData(trialDuration,fluorByTrial,runBefore,wheel_data,laser_data,shutter_data,tit,preOptoRunning)

figure();
subplot(3,1,1);
times=linspace(0,trialDuration,size(fluorByTrial,2));
fill([times fliplr(times)],[nanmean(fluorByTrial(preOptoRunning==runBefore,:),1)+nanstd(fluorByTrial(preOptoRunning==runBefore,:),[],1)./sqrt(size(fluorByTrial(preOptoRunning==runBefore,:),1)) fliplr(nanmean(fluorByTrial(preOptoRunning==runBefore,:),1)-nanstd(fluorByTrial(preOptoRunning==runBefore,:),[],1)./sqrt(size(fluorByTrial(preOptoRunning==runBefore,:),1)))],'c');
hold on;
plot(times,nanmean(fluorByTrial(preOptoRunning==runBefore,:),1));
axis tight;
xlim([0 max(times)]);
title(['Bulk Fluorescence - Average ' tit]);
subplot(3,1,2);
times=linspace(0,trialDuration,size(wheel_data,2));
plot(times,nanmean(wheel_data(preOptoRunning==runBefore,:),1));
axis tight;
xlim([0 max(times)]);
title(['Position - Average ' tit]);
subplot(3,1,3);
times=linspace(0,trialDuration,size(laser_data,2));
plot(times,nanmean(laser_data(preOptoRunning==runBefore,:),1)./100,'Color','c');
hold on;
plot(times,nanmean(shutter_data(preOptoRunning==runBefore,:),1),'Color','k');
axis tight;
xlim([0 max(times)]);
title(['Laser: Blue, Shutter: Black - Average ' tit]);

end

function trials_isRunning=findRunningFrames(wheel_trials,nFramesPerTrial)

ds=2000; % rate at which to downsample wheel encoder data before taking derivate for velocity
% velocity_thresh=0.01; % >= this velocity, consider mouse to be running
velocity_thresh=0.03; % >= this velocity, consider mouse to be running

trials_isRunning=nan(length(wheel_trials),nFramesPerTrial);
for i=1:length(wheel_trials)
    currRunningData=wheel_trials{i}.data;
    % Take derivate of downsampled encoder position data
    velocity=abs(diff(decimate(currRunningData,ds)));
    % Convert velocity into number of frames matching nFramesPerTrial
    resample_velocity=resample(velocity,nFramesPerTrial,length(velocity));
    isRunningFrame=resample_velocity>=velocity_thresh;
    trials_isRunning(i,:)=isRunningFrame;
end
end

