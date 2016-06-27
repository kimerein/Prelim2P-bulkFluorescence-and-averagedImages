function imageBeforeAndAfterStim(fileDir)

example=1;
thresh=4*10^7; % Threshold below which frame is considered to be blank
nFramesAfterOpto=1; % Number of frames to average after opto stim or "fake" opto stim (i.e., just shutter)
nFramesForBase=1; % Number of frames preceding opto stim to average for baseline image

% Read in files
listing=dir([fileDir '\*.tif']);
fileTimes=cell(1,length(listing));
for i=1:length(listing)
    fileTimes{i}=listing(i).date;
end
[orderedTimes,order]=sort(fileTimes);
listing=listing(order);

% Read in images for each trial
disp('num files');
disp(length(listing));
for i=1:length(listing)
    disp(i);
    fname=listing(i).name;
    info=imfinfo([fileDir '\' fname]);
    num_images=numel(info);
    for k=1:num_images
        frames{k}=double(imread([fileDir '\' fname],k,'Info',info));
    end
    trials{i}=frames;
end

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

% Assume that first shuttered trial is opto stim
baseImages=cell(1,length(trials));
afterRealStimImages=cell(1,length(trials));
afterOnlyShutterImages=cell(1,length(trials));
beforeShutterImages=cell(1,length(trials));
for i=1:length(trials)
    trial_frames=trials{i};
    shuttered_trials=trialByTrial_shuttered{i};
    if nShutteredTrials>2
        shuttered_trials=[mode_firstShuttered mode_secondShuttered];
    end
    % Calculate average image before opto stim
    summedIm=zeros(size(trial_frames{1}));
    n=0;
    for j=shuttered_trials(1)-1-nFramesForBase:shuttered_trials(1)-2
        summedIm=summedIm+trial_frames{j};
        n=n+1;
    end
    baseImages{i}=summedIm./n;
    % Calculate average image after opto stim
    summedIm=zeros(size(trial_frames{1}));
    n=0;
    for j=shuttered_trials(1)+nShutteredTrials-1+1:shuttered_trials(1)+nShutteredTrials-1+nFramesAfterOpto
        summedIm=summedIm+trial_frames{j};
        n=n+1;
    end
    afterRealStimImages{i}=summedIm./n;
    % Calculate average image after shutter only (no opto stim)
    summedIm=zeros(size(trial_frames{1}));
    n=0;
    for j=shuttered_trials(2)+1:shuttered_trials(2)+nFramesAfterOpto
        summedIm=summedIm+trial_frames{j};
        n=n+1;
    end
    afterOnlyShutterImages{i}=summedIm./n;
    % Calculate average image before shutter
    summedIm=zeros(size(trial_frames{1}));
    n=0;
    for j=shuttered_trials(2)-(nShutteredTrials-1)-1-nFramesForBase:shuttered_trials(2)-(nShutteredTrials-1)-2
        summedIm=summedIm+trial_frames{j};
        n=n+1;
    end
    beforeShutterImages{i}=summedIm./n;
end

% Plot example trial's images
figure();
imagesc(baseImages{example});
title('Baseline before opto stim');
colormap(gray);
figure();
imagesc(afterRealStimImages{example});
title('After opto stim');
colormap(gray);
figure();
imagesc(afterOnlyShutterImages{example});
title('After shutter only');
colormap(gray);
figure();
imagesc([baseImages{example} afterRealStimImages{example} afterOnlyShutterImages{example}]);
title('Altogether');
colormap(gray);
figure();
% imagesc((afterRealStimImages{example}-baseImages{example})./baseImages{example});
imagesc((afterRealStimImages{example}-baseImages{example}));
title('Difference after and before opto stim');    
colormap(gray);
figure();
% imagesc((afterOnlyShutterImages{example}-baseImages{example})./baseImages{example});
imagesc((afterOnlyShutterImages{example}-baseImages{example}));
title('Difference after and before shutter');
colormap(gray);

a=afterRealStimImages{example}-baseImages{example};
b=baseImages{example};
disp('corrcoef base and change w stim');
disp(corrcoef(a(1:end),b(1:end)));

a=afterOnlyShutterImages{example}-baseImages{example};
b=baseImages{example};
disp('corrcoef base and change w shutter');
disp(corrcoef(a(1:end),b(1:end)));

a=afterOnlyShutterImages{example}-beforeShutterImages{example};
b=beforeShutterImages{example};
disp('corrcoef before shutter and shutter minus before shutter');
disp(corrcoef(a(1:end),b(1:end)));

% Plot average images across trials
avIm=zeros(size(baseImages{1}));
for i=1:length(baseImages)
    avIm=avIm+baseImages{i};
end
av_baseImage=avIm./length(baseImages);

avIm=zeros(size(afterRealStimImages{1}));
for i=1:length(afterRealStimImages)
    avIm=avIm+afterRealStimImages{i};
end
av_afterStim=avIm./length(afterRealStimImages);

avIm=zeros(size(afterOnlyShutterImages{1}));
for i=1:length(afterOnlyShutterImages)
    avIm=avIm+afterOnlyShutterImages{i};
end
av_afterShutter=avIm./length(afterOnlyShutterImages);

avIm=zeros(size(beforeShutterImages{1}));
for i=1:length(beforeShutterImages)
    avIm=avIm+beforeShutterImages{i};
end
av_beforeShutter=avIm./length(beforeShutterImages);

figure();
imagesc(av_baseImage);
title('Av baseline before opto stim');
colormap(gray);
figure();
imagesc(av_afterStim);
title('Av after opto stim');
colormap(gray);
figure();
imagesc(av_afterShutter);
title('Av after shutter only');    
colormap(gray);
figure();
imagesc([av_baseImage av_afterStim av_afterShutter]);
title('Av altogether');    
colormap(gray);
figure();
% imagesc((av_afterStim-av_baseImage)./av_baseImage);
imagesc((av_afterStim-av_baseImage));
title('Av difference after and before opto stim');    
colormap(gray);
figure();
% imagesc((av_afterShutter-av_baseImage)./av_baseImage);
imagesc((av_afterShutter-av_baseImage));
title('Av difference before and after shutter');    
colormap(gray);
figure();
% imagesc([(av_afterStim-av_baseImage)./av_baseImage (av_afterShutter-av_baseImage)./av_baseImage]);
imagesc([(av_afterStim-av_baseImage) (av_afterShutter-av_beforeShutter)]);
title('Both av diffs');    
colormap(gray);
figure();
imagesc((av_afterStim-av_baseImage)-(av_afterShutter-av_beforeShutter));
title('Diff of av diffs');
colormap(gray);

a=av_afterStim-av_baseImage;
b=av_baseImage;
disp('corrcoef av base and change w stim');
disp(corrcoef(a(1:end),b(1:end)));

a=av_afterShutter-av_baseImage;
b=av_baseImage;
disp('corrcoef av base and change w shutter');
disp(corrcoef(a(1:end),b(1:end)));
        
    
    
    
    
    