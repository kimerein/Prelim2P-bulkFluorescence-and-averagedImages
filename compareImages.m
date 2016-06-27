function [output,outputrunning]=compareImages(fileDir)

% Set variables
encoderName='Wheel_Encoder'; % string of file name for analog encoder data

% Read in image files
image_trials=readTifFiles(fileDir);

% Read in running data
wheel_trials=readPhysFiles(fileDir,encoderName);

% Find opto stim and shutter trials
% [firstShutteredInd,secondShutteredInd,baselineFrames,postOptoFrames,preShutterFrames,postShutterFrames]=findShutteredFrames(image_trials,4*10^7,1,1);
% [firstShutteredInd,secondShutteredInd,baselineFrames,postOptoFrames,preShutterFrames,postShutterFrames]=findShutteredFrames(image_trials,1.7*10^7,1,1);
[firstShutteredInd,secondShutteredInd,baselineFrames,postOptoFrames,preShutterFrames,postShutterFrames]=findShutteredFrames(image_trials,0.6*10^7,1,1);

% Find frames where mouse is running
trials_isRunning=findRunningFrames(wheel_trials,length(image_trials{1}));

% Sort running vs. non-running image frames
% Also get running vs. non-running image average
running_image_trials=cell(1,length(image_trials));
nonRunning_image_trials=cell(1,length(image_trials));
exampleImagetrials=image_trials{1};
imageSize=size(exampleImagetrials{1});
running_avImage=zeros(imageSize);
running_countFrames=0;
nonRunning_avImage=zeros(imageSize);
nonRunning_countFrames=0;
allrunningimages=[];
allnonrunningimages=[];
for i=1:length(image_trials)
    trialImages=image_trials{i};
    runningFrames_thisTrial=trials_isRunning{i};
    currRunning_image_trials=cell(1,length(trialImages));
    currnonRunning_image_trials=cell(1,length(trialImages));
    for j=1:length(trialImages)
        if runningFrames_thisTrial(j)==1
            currRunning_image_trials{j}=trialImages{j};
            currnonRunning_image_trials{j}=nan(size(trialImages{j}));
            running_avImage=running_avImage+trialImages{j};
            running_countFrames=running_countFrames+1;
            allrunningimages(:,:,running_countFrames)=trialImages{j};
        else
            currnonRunning_image_trials{j}=trialImages{j};
            currRunning_image_trials{j}=nan(size(trialImages{j}));
            nonRunning_avImage=nonRunning_avImage+trialImages{j};
            nonRunning_countFrames=nonRunning_countFrames+1;
            allnonrunningimages(:,:,nonRunning_countFrames)=trialImages{j};
        end
    end
    running_image_trials{i}=currRunning_image_trials;
    nonRunning_image_trials{i}=currnonRunning_image_trials;
end
disp('Number of Running Frames');
disp(running_countFrames);
disp('Number of Non-Running Frames');
disp(nonRunning_countFrames);
running_avImage=running_avImage./running_countFrames;
nonRunning_avImage=nonRunning_avImage./nonRunning_countFrames;
outputrunning.allrunningimages=allrunningimages;
outputrunning.allnonrunningimages=allnonrunningimages;

% Display running vs. non-running image average
figure();
imagesc([running_avImage nonRunning_avImage]);
colormap('gray');
title('Left Image: Running, Right Image: Stationary');
            
figure();
imagesc(running_avImage-nonRunning_avImage);
colormap('gray');
title('Difference Between Running and Non-Running');

% Plot before vs. after opto stim and shutter for running vs. non-running
% Take average across trials of averaged frames within each trial
preOpto_running=zeros(imageSize);
postOpto_running=zeros(imageSize);
preOpto_nonrunning=zeros(imageSize);
postOpto_nonrunning=zeros(imageSize);
preShutter_running=zeros(imageSize);
postShutter_running=zeros(imageSize);
preShutter_nonrunning=zeros(imageSize);
postShutter_nonrunning=zeros(imageSize);

preMinusPostOpto_running=zeros(imageSize);
preMinusPostOpto_nonrunning=zeros(imageSize);
preMinusPostShutter_running=zeros(imageSize);
preMinusPostShutter_nonrunning=zeros(imageSize);

preOpto_running_n=0;
postOpto_running_n=0;
preOpto_nonrunning_n=0;
postOpto_nonrunning_n=0;
preShutter_running_n=0;
postShutter_running_n=0;
preShutter_nonrunning_n=0;
postShutter_nonrunning_n=0;

preMinusPostOpto_running_n=0;
preMinusPostOpto_nonrunning_n=0;
preMinusPostShutter_running_n=0;
preMinusPostShutter_nonrunning_n=0;

preOpto_run=nan(imageSize(1),imageSize(2),length(image_trials));
postOpto_run=nan(imageSize(1),imageSize(2),length(image_trials));
preOpto_nonrun=nan(imageSize(1),imageSize(2),length(image_trials));
postOpto_nonrun=nan(imageSize(1),imageSize(2),length(image_trials));
preShutter_run=nan(imageSize(1),imageSize(2),length(image_trials));
postShutter_run=nan(imageSize(1),imageSize(2),length(image_trials));
preShutter_nonrun=nan(imageSize(1),imageSize(2),length(image_trials));
postShutter_nonrun=nan(imageSize(1),imageSize(2),length(image_trials));
for i=1:length(image_trials)
    run_trialImages=running_image_trials{i};
    nonRun_trialImages=nonRunning_image_trials{i};
    
    % Running opto stim
    curr=averageImages(run_trialImages,baselineFrames{i});
    preOpto_run(:,:,i)=curr;
    if ~any(isnan(curr))
        preOpto_running=preOpto_running+curr;
        preOpto_running_n=preOpto_running_n+1;
    end
    
    curr=averageImages(run_trialImages,postOptoFrames{i});
    postOpto_run(:,:,i)=curr;
    if ~any(isnan(curr))
        postOpto_running=postOpto_running+curr;
        postOpto_running_n=postOpto_running_n+1;
    end
    
    curr=differenceOfaverageImages(run_trialImages,postOptoFrames{i},baselineFrames{i});
    if ~any(isnan(curr))
        preMinusPostOpto_running=preMinusPostOpto_running+curr;
        preMinusPostOpto_running_n=preMinusPostOpto_running_n+1;
    end
    
    % Nonrunning opto stim
    curr=averageImages(nonRun_trialImages,baselineFrames{i});
    preOpto_nonrun(:,:,i)=curr;
    if ~any(isnan(curr))
        preOpto_nonrunning=preOpto_nonrunning+curr;
        preOpto_nonrunning_n=preOpto_nonrunning_n+1;
    end
    
    curr=averageImages(nonRun_trialImages,postOptoFrames{i});
    postOpto_nonrun(:,:,i)=curr;
    if ~any(isnan(curr))
        postOpto_nonrunning=postOpto_nonrunning+curr;
        postOpto_nonrunning_n=postOpto_nonrunning_n+1;
    end
    
    curr=differenceOfaverageImages(nonRun_trialImages,postOptoFrames{i},baselineFrames{i});
    if ~any(isnan(curr))
        preMinusPostOpto_nonrunning=preMinusPostOpto_nonrunning+curr;
        preMinusPostOpto_nonrunning_n=preMinusPostOpto_nonrunning_n+1;
    end
    
    % Running shutter
    curr=averageImages(run_trialImages,preShutterFrames{i});
    preShutter_run(:,:,i)=curr;
    if ~any(isnan(curr))
        preShutter_running=preShutter_running+curr;
        preShutter_running_n=preShutter_running_n+1;
    end
    
    curr=averageImages(run_trialImages,postShutterFrames{i});
    postShutter_run(:,:,i)=curr;
    if ~any(isnan(curr))
        postShutter_running=postShutter_running+curr;
        postShutter_running_n=postShutter_running_n+1;
    end
    
    curr=differenceOfaverageImages(run_trialImages,postShutterFrames{i},preShutterFrames{i});
    if ~any(isnan(curr))
        preMinusPostShutter_running=preMinusPostShutter_running+curr;
        preMinusPostShutter_running_n=preMinusPostShutter_running_n+1;
    end
    
    % Nonrunning shutter
    curr=averageImages(nonRun_trialImages,preShutterFrames{i});
    preShutter_nonrun(:,:,i)=curr;
    if ~any(isnan(curr))
        preShutter_nonrunning=preShutter_nonrunning+curr;
        preShutter_nonrunning_n=preShutter_nonrunning_n+1;
    end
    
    curr=averageImages(nonRun_trialImages,postShutterFrames{i});
    postShutter_nonrun(:,:,i)=curr;
    if ~any(isnan(curr))
        postShutter_nonrunning=postShutter_nonrunning+curr;
        postShutter_nonrunning_n=postShutter_nonrunning_n+1;
    end
    
    curr=differenceOfaverageImages(nonRun_trialImages,postShutterFrames{i},preShutterFrames{i});
    if ~any(isnan(curr))
        preMinusPostShutter_nonrunning=preMinusPostShutter_nonrunning+curr;
        preMinusPostShutter_nonrunning_n=preMinusPostShutter_nonrunning_n+1;
    end
end
preOpto_running=preOpto_running./preOpto_running_n;
postOpto_running=postOpto_running./postOpto_running_n;
preOpto_nonrunning=preOpto_nonrunning./preOpto_nonrunning_n;
postOpto_nonrunning=postOpto_nonrunning./postOpto_nonrunning_n;
preShutter_running=preShutter_running./preShutter_running_n;
postShutter_running=postShutter_running./postShutter_running_n;
preShutter_nonrunning=preShutter_nonrunning./preShutter_nonrunning_n;
postShutter_nonrunning=postShutter_nonrunning./postShutter_nonrunning_n;
preMinusPostOpto_running=preMinusPostOpto_running./preMinusPostOpto_running_n;
preMinusPostOpto_nonrunning=preMinusPostOpto_nonrunning./preMinusPostOpto_nonrunning_n;
preMinusPostShutter_running=preMinusPostShutter_running./preMinusPostShutter_running_n;
preMinusPostShutter_nonrunning=preMinusPostShutter_nonrunning./preMinusPostShutter_nonrunning_n;

output.preOpto_run=preOpto_run;
output.postOpto_run=postOpto_run;
output.preOpto_nonrun=preOpto_nonrun;
output.postOpto_nonrun=postOpto_nonrun;

output.preShutter_run=preShutter_run;
output.postShutter_run=postShutter_run;
output.preShutter_nonrun=preShutter_nonrun;
output.postShutter_nonrun=postShutter_nonrun;

% Get p-value of paired comparison between images
pval_image=getImagePval(preOpto_run,postOpto_run);
displayPvalImage(pval_image,'Pval pre v post opto run',postOpto_running);
displayPvalImage(pval_image,'Pval pre v post opto run',preMinusPostOpto_running);
% displayPvalImage(pval_image,'Pval pre v post opto run - Delta F/F',preMinusPostOpto_running./preOpto_running);
pval_image=getImagePval(preOpto_nonrun,postOpto_nonrun);
displayPvalImage(pval_image,'Pval pre v post opto nonrun',postOpto_nonrunning);
displayPvalImage(pval_image,'Pval pre v post opto nonrun',preMinusPostOpto_nonrunning);
% displayPvalImage(pval_image,'Pval pre v post opto nonrun - Delta F/F',preMinusPostOpto_nonrunning./preOpto_nonrunning);
pval_image=getImagePval(preShutter_run,postShutter_run);
displayPvalImage(pval_image,'Pval pre v post shutter run',postShutter_running);
displayPvalImage(pval_image,'Pval pre v post shutter run',preMinusPostShutter_running);
% displayPvalImage(pval_image,'Pval pre v post shutter run - Delta F/F',preMinusPostShutter_running./preShutter_running);
pval_image=getImagePval(preShutter_nonrun,postShutter_nonrun);
displayPvalImage(pval_image,'Pval pre v post shutter nonrun',postShutter_nonrunning);
displayPvalImage(pval_image,'Pval pre v post shutter nonrun',preMinusPostShutter_nonrunning);
% displayPvalImage(pval_image,'Pval pre v post opto run - Delta F/F',preMinusPostShutter_nonrunning./preShutter_nonrunning);

% Display comparison images
figure();
imagesc([preOpto_running postOpto_running]);
colormap('gray');
title('Left: Pre Opto, Right: Post Opt (Running)');

figure();
imagesc([preOpto_nonrunning postOpto_nonrunning]);
colormap('gray');
title('Left: Pre Opto, Right: Post Opt (Non-Running)');

figure();
imagesc([preShutter_running postShutter_running]);
colormap('gray');
title('Left: Pre Shutter, Right: Post Shutter (Running)');

figure();
imagesc([preShutter_nonrunning postShutter_nonrunning]);
colormap('gray');
title('Left: Pre Shutter, Right: Post Shutter (Non-Running)');

figure();
imagesc([running_avImage-nonRunning_avImage preMinusPostOpto_running preMinusPostOpto_nonrunning preMinusPostShutter_running preMinusPostShutter_nonrunning]);
colormap('gray');
title('Run vs Nonrun, Post vs Pre Opto Run, Post vs Pre Opto Nonrun, Post vs Pre Shutter Run, Post vs Pre Shutter Nonrun');

figure();
imagesc([preMinusPostOpto_running preMinusPostShutter_running]);
colormap('gray');
title('Post vs Pre Opto Run, Post vs Pre Shutter Run');

figure();
imagesc([preMinusPostOpto_nonrunning preMinusPostShutter_nonrunning]);
colormap('gray');
title('Post vs Pre Opto Nonrun, Post vs Pre Shutter Nonrun');

end

function [avOfImages,allImages]=averageImages(imageStack,inds)

sumOfImages=zeros(size(imageStack{1}));
allImages=nan([size(imageStack{1}) length(inds)]);
j=0;
for i=1:length(imageStack)
    if ismember(i,inds)
        useInd=find(i,inds);
        if ~any(isnan(imageStack{i}))
            sumOfImages=sumOfImages+imageStack{i};
            allImages(:,:,useInd)=imageStack{i};
            j=j+1;
        end
    end
end
if j==0
    avOfImages=nan(size(imageStack{1}));
else
    avOfImages=sumOfImages./j;
end

end

function pval_image=getImagePval(im1,im2)

pval_image=nan(size(im1,1),size(im1,2));

for i=1:size(im1,1)
    for j=1:size(im1,2)
        if ~all(isnan(im1(i,j,:))) && ~all(isnan(im2(i,j,:))) && ~all(isnan(im1(i,j,:)) | isnan(im2(i,j,:)))
            pval_image(i,j)=signrank(reshape(im1(i,j,:),1,length(im1(i,j,:))),reshape(im2(i,j,:),1,length(im2(i,j,:))));
        else
            pval_image(i,j)=nan;
        end
    end
end

end

function displayPvalImage(im,tit,overlay)

if isempty(overlay)
    figure();
    imagesc(1-im);
    caxis([0.95,1]);
    colormap('gray');
    colorbar;
    title(tit);
else
    figure();
    imagesc(overlay);
    colormap(gray);
    colorbar;
    
%     decreaseIm=overlay;
%     decreaseIm(decreaseIm>0)=0;
%     decreaseIm=abs(decreaseIm);
    
    % Make a truecolor all-red image
%     red=cat(3,ones(size(overlay)),zeros(size(overlay)),zeros(size(overlay)));
%     hold on;
%     h1=imagesc(red);
%     hold off;
%     set(h1,'AlphaData',decreaseIm);
    
    % Make a truecolor all-blue image
    blue=cat(3,zeros(size(overlay)),zeros(size(overlay)),ones(size(overlay)));
    hold on;
    h=imagesc(blue);
    hold off;
    tryim=im;
    tryim(tryim>=0.05)=0.05;
    tryim=tryim-min(min(tryim));
    tryim=(tryim./max(max(tryim))).*0.3;
    tryim=0.3-tryim;
    set(h,'AlphaData',tryim);
    title(tit);
end

end

function diffOfAv=differenceOfaverageImages(imageStack,inds1,inds2)

sumOfImages1=zeros(size(imageStack{1}));
j1=0;
sumOfImages2=zeros(size(imageStack{1}));
j2=0;
for i=1:length(imageStack)
    if ismember(i,inds1)
        if ~any(isnan(imageStack{i}))
            sumOfImages1=sumOfImages1+imageStack{i};
            j1=j1+1;
        end
    end
    if ismember(i,inds2)
        if ~any(isnan(imageStack{i}))
            sumOfImages2=sumOfImages2+imageStack{i};
            j2=j2+1;
        end
    end
end
if j1==0 || j2==0
    diffOfAv=nan(size(imageStack{1}));
else
    diffOfAv=(sumOfImages1./j1)-(sumOfImages2./j2);
end

end
        
    

