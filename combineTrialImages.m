function combineTrialImages(fileDir,outputDir,outputName)

listing=dir([fileDir '\*.tif']);
fileTimes=cell(1,length(listing));
for i=1:length(listing)
    fileTimes{i}=listing(i).date;
end
[orderedTimes,order]=sort(fileTimes);
listing=listing(order);

j=1;
disp('num files');
disp(length(listing));
for i=1:length(listing)
    disp(i);
    fname=listing(i).name;
    info=imfinfo([fileDir '\' fname]);
    num_images=numel(info);
    for k=1:num_images
        frames{j}=imread([fileDir '\' fname],k,'Info',info);
        j=j+1;
    end
end

imwrite(frames{1},[outputDir '\' outputName '.tif']);
disp('num frames');
disp(length(frames));
for i=2:length(frames)
    disp(i);
    imwrite(frames{i},[outputDir '\' outputName '.tif'],'WriteMode','append');
end