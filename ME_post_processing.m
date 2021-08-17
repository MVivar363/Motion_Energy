function [totalME_norm_alignON, totalME_norm_alignRT] = ME_post_processing(data,totalME,residuals)

if nargin<3 || isempty(residuals)
    residuals = 0;
end

dur = data.duration;
max_dur = round(max(dur)*data.refresh); % in frames

% OR, only work with cohs in coh_set
% max_dur = round(max(dur(ismember(data.coh,coh_set)))*data.refresh); % in frames

% similarly:
% totalME = totalME(ismember(data.coh,coh_set));

% sanity check #1
%MVL: Basically checking that both the Original Length of the Trials is as
%long as the length of the Motion Energy vectors for trials (using Frames
%not time)
maxMElength = 0;
for n = 1:length(totalME)
    maxMElength = max([length(totalME{n}) maxMElength]);
end
if max_dur ~= maxMElength
    disp('houston we have a problem');
    keyboard
end

TotalME_alignON = nan(size(totalME,1),max_dur);
TotalME_alignRT = nan(size(totalME,1),max_dur);
for q = 1:size(totalME,1)
    TotalME_alignON(q,1:length(totalME{q})) = totalME{q};
    TotalME_alignRT(q,1:length(totalME{q})) = flipud(totalME{q});
end

%MVL: I can definitely clean this up
coh = data.coherence;
coh(data.direction==180) = -coh(data.direction==180); % signed coh???

if residuals
    % subtract out mean for each coh to get residuals
    % ** note this subtraction is done at each time step!
    cohs = unique(coh);
    meanMEbyCoh = nan(length(cohs),max_dur);
    for c = 1:length(cohs)
        meanMEbyCoh(c,:) = nanmean(TotalME_alignON(coh==cohs(c),:));
    end
    totalME_resid_alignON = nan(size(TotalME_alignON));
    totalME_resid_alignRT = nan(size(TotalME_alignRT));
    for n = 1:length(totalME)
        totalME_resid_alignON(n,:) = TotalME_alignON(n,:) - meanMEbyCoh(cohs==coh(n),:);
        totalME_resid_alignRT(n,:) = TotalME_alignRT(n,:) - meanMEbyCoh(cohs==coh(n),:);
    end

    totalME_norm_alignON = totalME_resid_alignON / max(abs(nanmean(totalME_resid_alignON,2)));
    totalME_norm_alignRT = totalME_resid_alignRT / max(abs(nanmean(totalME_resid_alignRT,2)));
else % in either case, normalize before pooling ^ _
    totalME_norm_alignON = TotalME_alignON / max(abs(nanmean(TotalME_alignON,2))); %MVL: scale all trials ME with the highest trial average ME.
    totalME_norm_alignRT = TotalME_alignRT / max(abs(nanmean(TotalME_alignRT,2)));
end

% sanity check #2: linear relationship between coh and ME?
ME = nanmean(totalME_norm_alignON,2); % avg 'momentary' ME, as opposed to cumulative; MVL: For each trial
% remove some nan trials (bad seeds?)
coh(isnan(ME)) = []; 
ME(isnan(ME)) = [];
[ME_m, ME_se] = calcGroupMean(ME, coh, unique(coh)); %MVL: Finding the true mean for each Coh (overall in task, not specific to trial) and its Standard Error
%keyboard
%b = regress(ME,[coh ones(size(coh))]);
figure(11); 
errorbar(unique(coh), ME_m, ME_se,'o-');
hold on; 
%plot(unique(coh),b(1)*unique(coh)+b(2),'g--');
xlabel('signed coh'); ylabel('avg motion energy (a.u.)');
title(['diam=' num2str(6) '  spd=' num2str(data.speed(1)/10)]);
%MVL: num2str(data.ap_diam(1)/10) -> for 'diam=' text above

