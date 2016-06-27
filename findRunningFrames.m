function trials_isRunning=findRunningFrames(wheel_trials,nFramesPerTrial)

ds=2000; % rate at which to downsample wheel encoder data before taking derivate for velocity
velocity_thresh=0.01; % >= this velocity, consider mouse to be running
% velocity_thresh=0.03; % >= this velocity, consider mouse to be running

trials_isRunning=cell(1,length(wheel_trials));
for i=1:length(wheel_trials)
    currRunningData=wheel_trials{i}.data;
    % Take derivate of downsampled encoder position data
    velocity=abs(diff(decimate(currRunningData,ds)));
    % Convert velocity into number of frames matching nFramesPerTrial
    resample_velocity=resample(velocity,nFramesPerTrial,length(velocity));
    isRunningFrame=resample_velocity>=velocity_thresh;
    trials_isRunning{i}=isRunningFrame;
end
