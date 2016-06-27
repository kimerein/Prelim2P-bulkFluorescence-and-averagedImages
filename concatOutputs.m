function [concatOutputFields,fnames]=concatOutputs(fileDir,fileName,variableName,concatDimension)

% fileDir is directory with files to concat
listing=dir([fileDir '\' fileName '.mat']);
for i=1:length(listing)
    a=load([fileDir '\' listing(i).name]);
    curr=a.(variableName);
    if i==1
        fnames=fieldnames(curr);
        concatOutputFields=cell(1,length(fnames));
    end
    if isempty(curr.(fnames{1}))
        curr.(fnames{1})=nan(size(curr.(fnames{2}),1),size(curr.(fnames{2}),2),0);
    end
    if isempty(curr.(fnames{2}))
        curr.(fnames{2})=nan(size(curr.(fnames{1}),1),size(curr.(fnames{1}),2),0);
    end
    for j=1:length(fnames)
        concatOutputFields{j}=cat(concatDimension,concatOutputFields{j},curr.(fnames{j}));
    end
end