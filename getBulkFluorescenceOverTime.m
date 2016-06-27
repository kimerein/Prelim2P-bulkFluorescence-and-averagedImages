function fluorByTrial=getBulkFluorescenceOverTime(fileDir)

% Read in image files
image_trials=readTifFiles(fileDir);

% Get trial-by-trial bulk fluorescence
fluorByTrial=nan(length(image_trials),length(image_trials{1}));
for i=1:length(image_trials)
    curr=image_trials{i};
    for j=1:length(curr)
        bulk=nanmean(nanmean(curr{j}));
        fluorByTrial(i,j)=bulk;
    end
end

