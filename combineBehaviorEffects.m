function combineBehaviorEffects(fileDir)

fileName='a2a*output';
variableName='output';

[concatOutputFields,fnames]=concatOutputs(fileDir,fileName,variableName,1);

times=linspace(0,31.232,size(concatOutputFields{1},2));
% figure();
% subplot(2,1,1);
% curr=concatOutputFields{1}';
% curr=curr(:,~isnan(curr(1,:)));
% plot(times,curr,'r');
% title('Running Before Opto Stim');
% subplot(2,1,2);
% curr=concatOutputFields{2}';
% curr=curr(:,~isnan(curr(1,:)));
% plot(times,curr,'r');

figure();
subplot(2,1,1);
curr=concatOutputFields{3}';
curr=curr(:,~isnan(curr(1,:)));
plot(times,curr,'b');
title('Non-Running Before Opto Stim');
subplot(2,1,2);
curr=concatOutputFields{4}';
curr=curr(:,~isnan(curr(1,:)));
plot(times,curr,'b');

% figure(); 
% subplot(2,1,1);
% curr=concatOutputFields{1}';
% curr=curr(:,~isnan(curr(1,:)));
% plot(times,nanmean(curr,2),'r');
% title('Running Before Opto Stim');
% subplot(2,1,2);
% curr=concatOutputFields{2}';
% curr=curr(:,~isnan(curr(1,:)));
% plot(times,curr,'r');

figure(); 
subplot(2,1,1);
curr=concatOutputFields{3}';
curr=curr(:,~isnan(curr(1,:)));
plot(times,nanmean(curr,2),'b');
title('Non-Running Before Opto Stim');
subplot(2,1,2);
curr=concatOutputFields{4}';
curr=curr(:,~isnan(curr(1,:)));
plot(times,curr,'b');