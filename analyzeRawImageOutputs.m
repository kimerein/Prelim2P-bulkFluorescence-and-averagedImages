function analyzeRawImageOutputs(fileDir)

% % Compare running and non-running data
% fileName='*RUNNING';
% variableName='runningOutput';
% concatDimension=3;
% 
% % Concat outputs
% [runningOutput,fnames]=concatOutputs(fileDir,fileName,variableName,concatDimension);
% 
% % Analyze for p-values
% % Average images
% avIms=cell(1,length(fnames));
% for i=1:length(fnames)
%     av=nanmean(runningOutput{i},concatDimension);
%     avIms{i}=reshape(av,size(av,1),size(av,2));
% end
% figure(); 
% imagesc([avIms{1} avIms{2}]);
% colormap(gray);
% title('L: Running, R: Nonrunning');
% pval_im=getImagePval(runningOutput{1},runningOutput{2},0);
% displayPvalImage(pval_im,'Running v Nonrunning',avIms{1}-avIms{2});
% 
% 
% return

% Compare effects of opto stim and shutter
fileName='*V';
variableName='output';
concatDimension=3;

% Concat outputs
[output,fnames]=concatOutputs(fileDir,fileName,variableName,concatDimension);
disp('Read in outputs');

% Analyze for p-values
% Average images
avIms=cell(1,length(fnames));
for i=1:length(fnames)
    av=nanmean(output{i},concatDimension);
    avIms{i}=reshape(av,size(av,1),size(av,2));
end
disp('Averaged images');

figure();
imagesc([avIms{1} avIms{2}]);
colormap(gray);
title('L: Pre, R: Post Opto (Run)');
pval_im=getImagePval(output{1},output{2},1);
displayPvalImage(pval_im,'Pre v Post Opto (Run)',avIms{2}-avIms{1});
disp('Fig done');

figure();
imagesc([avIms{3} avIms{4}]);
colormap(gray);
title('L: Pre, R: Post Opto (Nonrun)');
pval_im=getImagePval(output{3},output{4},1);
displayPvalImage(pval_im,'Pre v Post Opto (Nonrun)',avIms{4}-avIms{3});
disp('Fig done');

figure();
imagesc([avIms{5} avIms{6}]);
colormap(gray);
title('L: Pre, R: Post Shutter (Run)');
pval_im=getImagePval(output{5},output{6},1);
displayPvalImage(pval_im,'Pre v Post Shutter (Run)',avIms{6}-avIms{5});
disp('Fig done');

figure();
imagesc([avIms{7} avIms{8}]);
colormap(gray);
title('L: Pre, R: Post Shutter (Nonrun)');
pval_im=getImagePval(output{7},output{8},1);
displayPvalImage(pval_im,'Pre v Post Shutter (Nonrun)',avIms{8}-avIms{7});
disp('Fig done');

end

function pval_image=getImagePval(im1,im2,paired)

pval_image=nan(size(im1,1),size(im1,2));

for i=1:size(im1,1)
    for j=1:size(im1,2)
        if ~all(isnan(im1(i,j,:))) && ~all(isnan(im2(i,j,:)))
            if paired==1
                pval_image(i,j)=signrank(reshape(im1(i,j,:),1,length(im1(i,j,:))),reshape(im2(i,j,:),1,length(im2(i,j,:))));
%                 pval_image(i,j)=ranksum(reshape(im1(i,j,:),1,length(im1(i,j,:))),reshape(im2(i,j,:),1,length(im2(i,j,:))));
            else
                pval_image(i,j)=ranksum(reshape(im1(i,j,:),1,length(im1(i,j,:))),reshape(im2(i,j,:),1,length(im2(i,j,:))));
            end
        else
            pval_image(i,j)=nan;
        end
    end
end

end

function displayPvalImage(im,tit,overlay)

im(isnan(im))=1;

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