function trials=readPhysFiles(fileDir,fileName)

% Read in files
listing=dir([fileDir '\' fileName '*.mat']);
fileTimes=cell(1,length(listing));
for i=1:length(listing)
    fileTimes{i}=listing(i).date;
end
[orderedTimes,order]=sort(fileTimes);
listing=listing(order);

% Read in phys wave for each trial
disp('num files');
disp(length(listing));
for i=1:length(listing)
    disp(i);
    fname=listing(i).name;
    physFile=load([fileDir '\' fname]);
    fieldname=fieldnames(physFile);
    physWave=physFile.(fieldname{1});
    trials{i}=physWave;
end