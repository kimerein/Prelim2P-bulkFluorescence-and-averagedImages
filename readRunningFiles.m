function trials=readPhysFiles(fileDir,fileName)

% Read in files
listing=dir([fileDir '\' fileName '.mat']);
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