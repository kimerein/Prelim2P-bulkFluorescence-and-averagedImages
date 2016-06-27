function combineBulkFluorescenceVsBehavior(fileDir)

% Set variables
trialDuration=15.616;
baseSub=1; % whether to baseline subtract pre-opto stimulus baseline
norm=0; % whether to normalize to peak fluorescence on given trial

% Read in files
listing=dir([fileDir '\*.mat']);
fileTimes=cell(1,length(listing));
for i=1:length(listing)
    fileTimes{i}=listing(i).date;
end
[orderedTimes,order]=sort(fileTimes);
listing=listing(order);

% Read in data for each file
disp('num files');
disp(length(listing));
concat_fluorByTrial_preRun=[];
concat_wheel_data_preRun=[];
concat_laser_data_preRun=[];
concat_shutter_data_preRun=[];
concat_fluorByTrial_preNonrun=[];
concat_wheel_data_preNonrun=[];
concat_laser_data_preNonrun=[];
concat_shutter_data_preNonrun=[];
for i=1:length(listing)
    disp(i);
    fname=listing(i).name;
    currFile=load([fileDir '\' fname]);
    currFile=currFile.output;
    
    if i~=1 && size(currFile.fluorByTrial_preNonrun,2)~=size(concat_fluorByTrial_preNonrun,2)
        % Fit sizes
        currFile.fluorByTrial_preRun=resizeMatrix(currFile.fluorByTrial_preRun,concat_fluorByTrial_preRun);
        currFile.wheel_data_preRun=resizeMatrix(currFile.wheel_data_preRun,concat_wheel_data_preRun);
        currFile.laser_data_preRun=resizeMatrix(currFile.laser_data_preRun,concat_laser_data_preRun);
        currFile.shutter_data_preRun=resizeMatrix(currFile.shutter_data_preRun,concat_shutter_data_preRun);
        
        currFile.fluorByTrial_preNonrun=resizeMatrix(currFile.fluorByTrial_preNonrun,concat_fluorByTrial_preNonrun);
        currFile.wheel_data_preNonrun=resizeMatrix(currFile.wheel_data_preNonrun,concat_wheel_data_preNonrun);
        currFile.laser_data_preNonrun=resizeMatrix(currFile.laser_data_preNonrun,concat_laser_data_preNonrun);
        currFile.shutter_data_preNonrun=resizeMatrix(currFile.shutter_data_preNonrun,concat_shutter_data_preNonrun);    
    end    
    
    concat_fluorByTrial_preRun=[concat_fluorByTrial_preRun; currFile.fluorByTrial_preRun];
    concat_wheel_data_preRun=[concat_wheel_data_preRun; currFile.wheel_data_preRun];
    concat_laser_data_preRun=[concat_laser_data_preRun; currFile.laser_data_preRun];
    concat_shutter_data_preRun=[concat_shutter_data_preRun; currFile.shutter_data_preRun];
    
    concat_fluorByTrial_preNonrun=[concat_fluorByTrial_preNonrun; currFile.fluorByTrial_preNonrun];
    concat_wheel_data_preNonrun=[concat_wheel_data_preNonrun; currFile.wheel_data_preNonrun];
    concat_laser_data_preNonrun=[concat_laser_data_preNonrun; currFile.laser_data_preNonrun];
    concat_shutter_data_preNonrun=[concat_shutter_data_preNonrun; currFile.shutter_data_preNonrun];
end

% Baseline subtract
if baseSub==1
    % Find first opto stim
    consensusLaser=nanmean(concat_laser_data_preNonrun,1);
    c=consensusLaser(1:floor(length(consensusLaser)/size(concat_fluorByTrial_preRun,2)):end);
    f=find(c>50,1,'first');
    for i=1:size(concat_fluorByTrial_preRun,1)
        concat_fluorByTrial_preRun(i,:)=concat_fluorByTrial_preRun(i,:)-nanmean(concat_fluorByTrial_preRun(i,1:f-1));
    end
    for i=1:size(concat_fluorByTrial_preNonrun,1)
        concat_fluorByTrial_preNonrun(i,:)=concat_fluorByTrial_preNonrun(i,:)-nanmean(concat_fluorByTrial_preNonrun(i,1:f-1));
    end
end

% Peak-norm
if norm==1
    % Find first opto stim
    consensusLaser=nanmean(concat_laser_data_preNonrun,1);
    c=consensusLaser(1:floor(length(consensusLaser)/size(concat_fluorByTrial_preRun,2)):end);
    f=find(c>50,1,'first');
    for i=1:size(concat_fluorByTrial_preRun,1)
%         concat_fluorByTrial_preRun(i,:)=concat_fluorByTrial_preRun(i,:)./nanmax(concat_fluorByTrial_preRun(i,:));
        concat_fluorByTrial_preRun(i,:)=concat_fluorByTrial_preRun(i,:)./nanmean(concat_fluorByTrial_preRun(i,1:f-1));
    end
    for i=1:size(concat_fluorByTrial_preNonrun,1)
%         concat_fluorByTrial_preNonrun(i,:)=concat_fluorByTrial_preNonrun(i,:)./nanmax(concat_fluorByTrial_preNonrun(i,:));
        concat_fluorByTrial_preNonrun(i,:)=concat_fluorByTrial_preNonrun(i,:)./nanmean(concat_fluorByTrial_preNonrun(i,1:f-1));
    end
end

% Plot average across all trials separated by running or non-running before
% opto stim
tit='Running Before Opto';
plotData(trialDuration,concat_fluorByTrial_preRun,concat_wheel_data_preRun,concat_laser_data_preRun,concat_shutter_data_preRun,tit);
tit='Stationary Before Opto';
plotData(trialDuration,concat_fluorByTrial_preNonrun,concat_wheel_data_preNonrun,concat_laser_data_preNonrun,concat_shutter_data_preNonrun,tit);


end


function curr=resizeMatrix(curr,concated)

downSampFactor=ceil(size(curr,2)./size(concated,2));
temp=zeros(size(curr,1),size(concated,2));
for i=1:size(curr,1)
    ds=decimate(curr(i,:),downSampFactor);
    temp(i,1:length(ds))=ds;
end
curr=temp;
        
end

function plotData(trialDuration,fluorByTrial,wheel_data,laser_data,shutter_data,tit)

figure();
subplot(3,1,1);
times=linspace(0,trialDuration,size(fluorByTrial,2));
fill([times fliplr(times)],[nanmean(fluorByTrial,1)+nanstd(fluorByTrial,[],1)./sqrt(size(fluorByTrial,1)) fliplr(nanmean(fluorByTrial,1)-nanstd(fluorByTrial,[],1)./sqrt(size(fluorByTrial,1)))],'c');
hold on;
plot(times,nanmean(fluorByTrial,1));
axis tight;
xlim([0 max(times)]);
title(['Bulk Fluorescence - Average ' tit]);
subplot(3,1,2);
times=linspace(0,trialDuration,size(wheel_data,2));
plot(times,nanmean(wheel_data,1));
axis tight;
xlim([0 max(times)]);
title(['Position - Average ' tit]);
subplot(3,1,3);
times=linspace(0,trialDuration,size(laser_data,2));
plot(times,nanmean(laser_data,1)./100,'Color','c');
hold on;
plot(times,nanmean(shutter_data,1),'Color','k');
axis tight;
xlim([0 max(times)]);
title(['Laser: Blue, Shutter: Black - Average ' tit]);

end
