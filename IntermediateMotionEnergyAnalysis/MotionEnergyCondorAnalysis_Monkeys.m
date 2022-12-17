% Miguel Vivar-Lazo
% 02/24/2020
% Fetsch Lab
% Checking MotionEnergy of Partition Data by Dates

%%
% Load the files you want from the Condor Stuff
top = 1;  
bestFiles = 1; %best2worst(1:top); %contains the best to worst files in the condor (in motionenergycondoranalysis)

for ii = 1:length(bestFiles)
%filename = sprintf('CondorMEGenji%d.mat',bestFiles(ii)); %iterate through file, WHEN IT IS SOURCE FROM CONDOR
%filename = sprintf('Hanzo_Training_ALL_withME.mat');
%load(filename); %load file
intCoh = 0; %.512;

%zscore all the ME BY COHERENCE
[matME, matDVfromME, endDVfromME] = matrixConverterME(totalME);
residualME = normalizeME(matME,'center',coherence.*cosd(direction)); %Find residuals; cumsum this to get the average of matDVfromME
normMatME = residualME; %normalizeME(residualME,'zscore',coherence.*cosd(direction)); %Find residuals
  
clusters = 5;
avgMERight = cell(1,clusters);
avgMELeft = cell(1,clusters);
%fix cluster so its from best to worst
[~, indexMu] = sort(mu, 'descend');
organizeClusters = indexMu;
for i = organizeClusters %the cluster number does not reflect whether its the best or worst one, you have to look at 'mu' for that
    daysForCluster = uniqueDays(assignment == i); %assignments is from "Monkey_TrainingInDotsModels"
    MotionEnergyforDaysRight = normMatME(ismember(dates, daysForCluster) & choice == 2 & coherence == intCoh, :); %pick normalize ME for certain days (that belong to the cluster of interest)
    MotionEnergyforDaysLeft = normMatME(ismember(dates, daysForCluster) & choice == 1 & coherence == intCoh, :);
    %select between choices that are left vs right
    avgMERight{i} = nanmean(MotionEnergyforDaysRight); 
    avgMELeft{i} = nanmean(MotionEnergyforDaysLeft); 
end


%include time, and plot all cluster between right and left
time = 0:.0083:(size(normMatME,2)-1)*.0083;
stepColor = 1/(clusters);
color = 0;
figure(41+ii);
counter = 0;
hold on;
for i = organizeClusters(2:5)
    counter = counter+1;
    subplot(2,2,counter)
    plot(time, avgMERight{i}, 'r-'); hold on; % 'Color', [1-color, 0 ,0]);
    plot(time, avgMELeft{i}, 'b-'); hold off; %'Color', [0, 0 , 1 - color]); 
    color = color + stepColor;
end
hold off
xlim([0,1])
title(['Motion Energy: Right vs Left for ' num2str(ii)])
end

%% Graph Right vs Left
% but use the norm of the ME timeseries of each coherence, not jsut
% coherence (Warning: if not enough data for the later times, this causes
% the Normalization to be incredibly wrong because it makes the STE much
% bigger for later trials causes the normalization to be much bigger than
% earlier times that have more data, try using STD although that might be
% wrong in a statistical method standpoint) This problem generalizes for
% all usage of STE in any of the next analysis

% Load the files you want from the Condor Stuff
top = 1;  
bestFiles = 1; % best2worst(1:top); %contains the best to worst files in the condor (in motionenergycondoranalysis)

for ii = 1:length(bestFiles)
% filename = sprintf('CondorMEGenji%d.mat',bestFiles(ii)); %iterate through file
% load(filename); %load file
intCoh = 0;

%zscore all the ME BY COHERENCE
[matME, matDVfromME, endDVfromME] = matrixConverterME(totalME);
residualME = normalizeME(matDVfromME,'center',coherence.*cosd(direction)); %Find residuals; cumsum this to get the average of matDVfromME
DVmatME = residualME; %normalizeME(residualME,'zscore',coherence.*cosd(direction)); %Find residuals
   
clusters = 4;
avgMERight = cell(1,clusters);
avgMELeft = cell(1,clusters);
%fix cluster so its from best to worst
[~, indexMu] = sort(mu, 'descend');
organizeClusters = indexMu;
for i = organizeClusters %the cluster number does not reflect whether its the best or worst one, you have to look at 'mu' for that
    daysForCluster = uniqueDays(assignment == i); %assignments is from "Monkey_TrainingInDotsModels"
    MotionEnergyforDaysRight = DVmatME(ismember(dates, daysForCluster) & choice == 2 & coherence == intCoh, :); %pick normalize ME for certain days (that belong to the cluster of interest)
    MotionEnergyforDaysLeft = DVmatME(ismember(dates, daysForCluster) & choice == 1 & coherence == intCoh, :);
    %select between choices that are left vs right
    avgMERight{i} = nanmean(MotionEnergyforDaysRight); 
    avgMELeft{i} = nanmean(MotionEnergyforDaysLeft); 
end


%include time, and plot all cluster between right and left
time = 0:.0083:(size(normMatME,2)-1)*.0083;
stepColor = 1/(clusters);
color = 0;
figure(91+ii);
counter = 0;
hold on;
for i = organizeClusters
    counter = counter+1;
    subplot(2,2,counter)
    plot(time, avgMERight{i}, 'Color', [1-color, 0 ,0]); hold on;
    plot(time, avgMELeft{i}, 'Color', [0, 0 , 1 - color]); hold off;
    color = color + stepColor;
end
hold off
xlim([0,1])
title(['Motion Energy: Right vs Left for ' num2str(ii)])
end

%% Graph Correct vs Incorrect
% We will grab the residuals for each Coherence ME, and map the residual
% based on Correct for Incorrect for that Sign COherence

% Load the files you want from the Condor Stuff
top = 1;  
bestFiles = best2worst(1:top); %contains the best to worst files in the condor (in motionenergycondoranalysis)

for ii = 1:length(bestFiles)
    filename = sprintf('CondorMEGenji%d.mat',bestFiles(ii)); %iterate through file
    load(filename); %load file
    
    %we will normalize the residuals
    [matME, matDVfromME, endDVfromME] = matrixConverterME(totalME);
    residualMatME = ones(size(matME)); %Find residuals for all trials based on Coherence, Residuals are calcualted by Mean of ME through time. 
    for i = 1:length(signUniqueCoh) %Find residuals for each coherence (maybe not coherence just all coherence together?) BUT BY EACH TIME. So ME for coherence .512 at time sample 1 would all be normalize with all .512, time sample 2 would be normalize with all sample 2 in coherence .512
        tempME = matME(coherence.*cos(deg2rad(direction)) == signUniqueCoh(i),:); %pick ME that matter based on coherence
        avgME = nanmean(tempME); %average ME
        residualME = tempME - avgME; %residual
        avgResME = nanmean(residualME); %avg residuals
        stdResME = nanstd(residualME); %residuals std
        lengthResME = sum(~isnan(tempME)); %Residuals length
        steResME = stdResME./sqrt(lengthResME); %Residuals STE
        residualMatME(coherence.*cos(deg2rad(direction)) == signUniqueCoh(i),:) = (residualME - avgResME)./stdResME; %residuals, placing new values into original position
    end

    %fix cluster so its from best to worst
    [~, indexMu] = sort(mu, 'descend');
    organizeClusters = indexMu;
    cluster = length(organizeClusters);
    avgMERightCorr = cell(1,clusters);
    avgMERightInc = cell(1,clusters);
    avgMELeftCorr = cell(1,clusters);
    avgMELeftInc = cell(1,clusters);
    
    %Normalize all the residuals in order to combine all the Residuals
    %between coherences
    temp = unique(coherence);
    uniqueCohWout0 = temp(temp ~= 0);  
    for i = organizeClusters %the cluster number does not reflect whether its the best or worst one, you have to look at 'mu' for that
        MotionEnergyforDaysRightCorr = [];
        MotionEnergyforDaysRightInc = [];
        MotionEnergyforDaysLeftCorr = [];
        MotionEnergyforDaysLeftInc = [];
        for intCoh = uniqueCohWout0
            daysForCluster = uniqueDays(assignment == i); %assignments is from "Monkey_TrainingInDotsModels"
            %select between choice correct vs incorrect for sign coh
            MotionEnergyforDaysRightCorr = [MotionEnergyforDaysRightCorr ; residualMatME(ismember(dates, daysForCluster) & coherence.*cos(deg2rad(direction)) == intCoh & correct == 1, :)]; %pick normalize ME for certain days (that belong to the cluster of interest)
            MotionEnergyforDaysRightInc = [MotionEnergyforDaysRightInc ; residualMatME(ismember(dates, daysForCluster)& coherence.*cos(deg2rad(direction)) == intCoh & correct == 0, :)];
            %select betwee choice that are left and right and correct vs incorrect
            MotionEnergyforDaysLeftCorr = [MotionEnergyforDaysLeftCorr ; residualMatME(ismember(dates, daysForCluster) & coherence.*cos(deg2rad(direction)) == -intCoh & correct == 1, :)]; %pick normalize ME for certain days (that belong to the cluster of interest)
            MotionEnergyforDaysLeftInc = [MotionEnergyforDaysLeftInc ; residualMatME(ismember(dates, daysForCluster) & coherence.*cos(deg2rad(direction)) == -intCoh & correct == 0, :)];
        end
        %select between choices that are left vs right and correct vs incorrect
        avgMERightCorr{i} = nanmean(MotionEnergyforDaysRightCorr); 
        avgMERightInc{i} = nanmean(MotionEnergyforDaysRightInc);
        %select between choices that are left vs right and correct vs incorrect
        avgMELeftCorr{i} = nanmean(MotionEnergyforDaysLeftCorr); 
        avgMELeftInc{i} = nanmean(MotionEnergyforDaysLeftInc);
    end
    
    %Plot Incorrect vs Correct ME trials 
    color = 0;
    figure(51+ii);
    time = 0:.0083:(size(normMatME,2)-1)*.0083;
    hold on;
    for i = organizeClusters
        plot(time, avgMERightCorr{4}, '-', 'Color', [1-color, 0 ,0]); 
        plot(time, avgMERightInc{4}, '--', 'Color', [1-color, 0 ,0]); 
        color = color + 0;%stepColor;
    end
    hold off
    xlim([0,1])
    title(['Motion Energy: Right with Correct or Incorrect for ' num2str(ii)])
    
    %Plot Incorrect vs Correct ME trials 
    color = 0;
    figure(71+ii);
    hold on;
    for i = organizeClusters
        plot(time, avgMELeftCorr{4}, '-', 'Color', [0, 0 , 1 - color]);
        plot(time, avgMELeftInc{4}, '--', 'Color', [0, 0 , 1 - color]);
        color = color + 0;%stepColor;
    end
    hold off
    xlim([0,1])
    title(['Motion Energy: Right with Correct or Incorrect for ' num2str(ii)])
end
%% Understand the Correlation matrix of individual samples through time
% Load the files you want from the Condor Stuff
top = 1;  
bestFiles = best2worst(1:top); %contains the best to worst files in the condor (in motionenergycondoranalysis)

for ii = 1:length(bestFiles)
 % Find autocorrelation of each ME belonging to Coh
    uniqCoh = unique(coherence.*cos(deg2rad(direction)));
    autoCorrME = cell(size(uniqCoh));
    PCorrME = cell(size(uniqCoh));
    lengthCorrME = cell(size(uniqCoh));
    for i = 1:length(uniqCoh)
        matForCorr = matME(coherence.*cos(deg2rad(direction)) == uniqCoh(i), :);
        [corrMat, PcorrMat] = corrcoef(matForCorr, 'Rows', 'pairwise'); %remove nan
        %meanME(randperm(length(meanME)));%shuffle data to see if we can get it to not be autocorrelated
        autoCorrME{i} = corrMat;
        PCorrME{i} = PcorrMat;
        
        countCorrLength = 0; %What is the length of the correlation?
        for c = 1:size(PCorrME{i},2) %Go through all the sample times (RV) and check for significant correlation between the analyse sample
            tempCorrVec = PCorrME{i}(c,:); %Correlation to ith's RV
            tempSigCorr = find(tempCorrVec <= .05); %find significant correlations 
            tempCorr2Index = length(find(abs(c - tempSigCorr) <=10)) ; % find the ones around the ith RV of interest, how many are there?
            countCorrLength = countCorrLength + tempCorr2Index; %add one because the second 'diff' will cause you to always lose one corr sample            
        end
        
        lengthCorrME{i} = countCorrLength/size(PCorrME{i},2); %find the average length
        %seems to suggest that after a window a 4 there is not that much
        %correlation
    end
end  

%Find the average correlation matrix and construct heat map
avgCorr = NaN(size(autoCorrME{1},1), size(autoCorrME{1},2), length(autoCorrME));
for i = 1:length(autoCorrME)
    autoCorrME{i}(PCorrME{i} > .05 & PCorrME{i} ~= 1) = NaN;
    avgCorr(:,:,i) = autoCorrME{i};
end
avgCorr_real = nanmean(avgCorr,3);
%avgCorr_real(isnan(avgCorr_real)) = 0;
figure(38595839);
surf((1:length(avgCorr_real)).*.0084, (1:length(avgCorr_real)).*.0084, avgCorr_real);
colorbar
title('Correlation Map')
xlabel('Time Sample')
ylabel('Time Sample')
%% Pass a sliding window through the ME profiles in order to remove correlations
window = 6; %Three indeces to the left and 3 indeces to the right
NoCorrMatME = NaN(size(matME,1), ceil(size(matME,2)/(window*2+1)));
% will it be a sliding window that keeps the matrix the same size as
% original ME matrix, or one that downsamples? Depends on the outcome
% correlation matrix
indexAvg = 1: window*2+1: size(matME,2); %index that will be average with its surroundings
for ii = 1:length(bestFiles)
    for i = 1:length(indexAvg)
        if indexAvg(i) - window < 1
            NoCorrMatME(:, i) = matME(:, 1); %nanmean([NaN(size(matME, 1), abs(indexAvg(i)-window) + 1) matME(:, 1 : indexAvg(i) + window)], 2);
        elseif indexAvg(i) + window > size(matME, 2)
            NoCorrMatME(:, i) = matME(:, end); %nanmean([matME(:, indexAvg(i) - window : end) NaN(size(matME, 1), abs(size(matME, 2) - (indexAvg(i)+window))) ], 2);
        else
            NoCorrMatME(:, i) = matME(:, indexAvg(i)); %nanmean(matME(:, indexAvg(i) - window : indexAvg(i) + window), 2);
        end
    end
end  
%check that most correlations were eliminated
uniqCoh = unique(coherence.*cos(deg2rad(direction)));
CorrME = cell(size(uniqCoh));
PCorrME = cell(size(uniqCoh));
for i = 1:length(uniqCoh)
        matForCorr = NoCorrMatME(coherence.*cos(deg2rad(direction)) == uniqCoh(i), :);
        [corrMat, PcorrMat] = corrcoef(matForCorr, 'Rows', 'pairwise'); %remove nan
        %meanME(randperm(length(meanME)));%shuffle data to see if we can get it to not be autocorrelated
        CorrME{i} = corrMat;
        PCorrME{i} = PcorrMat;
       
end

%Instead of averaging maybe just pick out samples
%% Use a simple logistic regression to display which time points effect Choices most.
% Use a regression to show that over the days there is structure on how the
% ME timeseries gives rise to Monkey's answer

% Load the files you want from the Condor Stuff
top = 1;  
bestFiles = best2worst(1:top); %contains the best to worst files in the condor (in motionenergycondoranalysis)

for ii = 1:length(bestFiles)
    %zscore all the ME BY COHERENCE
    NoCorrMatME; %[matME, matDVfromME, endDVfromME] = matrixConverterME(totalME);
    normMatME = ones(size(NoCorrMatME));
    for i = 1:length(signUniqueCoh)
        tempME = NoCorrMatME(coherence.*cos(deg2rad(direction)) == signUniqueCoh(i),:); %pick ME that matter based on coherence
        %normMatME(coherence.*cos(deg2rad(direction)) == signUniqueCoh(i),:) = normalize(tempME, 'scale', 'std'); %normalize columns
        avgME = nanmean(tempME); %average ME
        stdME = nanstd(tempME); %standard dev
        lengthME = sum(~isnan(tempME)); %length
        steME = stdME./sqrt(lengthME); %standard error
        %steME(steME < .2) = 1; %if any ste are really small fix that
        normMatME(coherence.*cos(deg2rad(direction)) == signUniqueCoh(i),:) = (tempME - 0)./stdME; %normalize, placing new values into original position
    end

    %fix cluster so its from best to worst
    [~, indexMu] = sort(mu, 'descend');
    organizeClusters = indexMu;
    cluster = length(organizeClusters);
    coeff_Choices = cell(1,cluster);
    for i = organizeClusters %the cluster number does not reflect whether its the best or worst one, you have to look at 'mu' for that
        daysForCluster = uniqueDays(assignment == i); %assignments is from "Monkey_TrainingInDotsModels"
        %find weights for choices using ME using a logistic regression
        timeIWant = 7; %when window is 3, 13; ...%80;
        trialsWithTime = ~isnan(normMatME(: ,timeIWant)); %trials that have enough information
        normMatME_Partition = normMatME(trialsWithTime' & ismember(dates, daysForCluster) & coherence~= .75, :); %pick ME that has enough information in the time domain (before it turns into NaN)
        normMatME_Partition = normMatME_Partition(:,1:timeIWant); %squeeze out the NaN to the number you want
        choiceForLogistic = choice(trialsWithTime' & ismember(dates, daysForCluster) & coherence~= .75); %choices for the ME trials you want
        choiceForLogistic = choiceForLogistic - 1; %For some reason changing the outcomes (right = 2 is not right =1 and Left = 1 is now Left =1) works better
        choiceForLogistic(choiceForLogistic == 0) = 2;
        %shuffleTime = randperm(size(normMatME_Partition,2)); %i should shuffle the data like this
        %shuffleNormMatME = normMatME_Partition(:, shuffleTime);%shuffle data to remove autocorrelation (correlation might still exist that is not time series dependent
        %shuffleChoiceForLogistic = choiceForLogistic(shuffleTime);
        [tempCoeff, dev, stats] = mnrfit(normMatME_Partition, choiceForLogistic'); %The coefficients that go along with this Logistic Regression
        tempCoeff = tempCoeff(2:end); %forget bias term
        coeff_Choices{i} = tempCoeff;%(shuffleTime); %unscramble coeff
    end
end

%Plot Incorrect vs Correct ME trials 
color = 0;
time = 0: .0083.*(window*2 + 1): (timeIWant-1)*.0083.*(window*2 + 1);
plotIndex = 0;
figure(101+plotIndex);
hold on;
for i = organizeClusters
    plot(time, coeff_Choices{i}, 'o--', 'Color', [i == 1, i ==2 , i ==3]);
    axis tight
    color = color + stepColor;
    plotIndex = plotIndex +1;
end
hold off;
title(['Motion Energy: Right with Correct or Incorrect for ' num2str(ii)])
xlabel('time')
ylabel('Coefficient Values')
legend('Best Cluster', '2nd', '3rd', '4th')
%% Model to see if Outlier or Accumulation explains data best
[matME, matDVfromME, endDVfromME] = matrixConverterME(totalME);
normMatME = ones(size(matME));
for i = 1:length(signUniqueCoh)
        tempME = matME(coherence.*cos(deg2rad(direction)) == signUniqueCoh(i),:); %pick ME that matter based on coherence
        normMatME(coherence.*cos(deg2rad(direction)) == signUniqueCoh(i),:) = normalize(tempME, 'scale', 'std'); %normalize columns
end
normMatME(normMatME == inf) = 0;
[~, indexMu] = sort(mu, 'descend');
organizeClusters = indexMu;

%a lot of the variables here are defined on the above box
index = 1;
boundAcc = zeros(size(organizeClusters));
boundHeight = zeros(size(organizeClusters));

signCoherence = coherence.*cos(deg2rad(direction));
uniqSignCoh = unique(signCoherence);
for i = organizeClusters %the cluster number does not reflect whether its the best or worst one, you have to look at 'mu' for that
    daysForCluster = uniqueDays(assignment == i); %assignments is from "Monkey_TrainingInDotsModels"
    condition = ismember(dates, daysForCluster) & ismember(coherence, [0]);
    data{1} = (choice(condition) - 1).*2 - 1; %have to subtract -1 since the function file works with 0s and 1s
    data{2} = matDVfromME(condition, :); %DV for each trial, in normalize form
    type = 'non-RT';
    params = 2;
    tau05 = 0; %not being used but needed for function to work
    options = optimset('MaxFunEvals', 1000, 'MaxIter', 1000, 'TolFun', 1e-4, 'TolX', 1e-4);
    %
    for z = 1000:1000:10000
        params = z; %try many different param starting points
        [bestBound, err_] = fminsearch(@(params) SpecificBound(data,type,tau05,params), params, options);
        if z == 1000 || err_ < errorModel
            errorModel = err_;
            boundModelChoice = bestBound;
        end
    end
    [~, fitChoices] = SpecificBound(data,type,tau05,boundModelChoice); %Find the fit
    boundAcc(index) = sum(fitChoices == data{1})./length(fitChoices); %Percentage the model got right
    boundHeight(index) = boundModelChoice;
    
    signCoherenceInt = signCoherence(condition);
    uniqSignCoh = unique(signCoherenceInt);
    avgModelChoice = zeros(size(uniqSignCoh));
    avgRealChoice = zeros(size(uniqSignCoh));
    for c = 1:length(uniqSignCoh)
        avgModelChoice(c) = nanmean(fitChoices(signCoherenceInt == uniqSignCoh(c)));
        avgRealChoice(c) = nanmean(data{1}(signCoherenceInt == uniqSignCoh(c)));
    end
 
    figure(475894 + index);
    plot(uniqSignCoh, avgModelChoice, 'bo--', uniqSignCoh, avgRealChoice, 'ro--')
    xlabel('Coherence')
    ylabel('Prob Right')
    title('Probability of Right with Fit: Bound')
    legend('Model', 'Empirical')
    
    index = index + 1;
end
%% Outlier Fitting 
[matME, matDVfromME, endDVfromME] = matrixConverterME(totalME);
normMatME = ones(size(matME));
for i = 1:length(signUniqueCoh)
        tempME = matME(coherence.*cos(deg2rad(direction)) == signUniqueCoh(i),:); %pick ME that matter based on coherence
        normMatME(coherence.*cos(deg2rad(direction)) == signUniqueCoh(i),:) = normalize(tempME, 'scale', 'std'); %normalize columns
end
normMatME(normMatME == inf) = 0;
[~, indexMu] = sort(mu, 'descend');
organizeClusters = indexMu;

%a lot of the variables here are defined on the above box
index = 1;
OutlierAcc = zeros(size(organizeClusters));
OutlierHeight = zeros(size(organizeClusters));

signCoherence = coherence.*cos(deg2rad(direction));
uniqSignCoh = unique(signCoherence);
for i = organizeClusters %the cluster number does not reflect whether its the best or worst one, you have to look at 'mu' for that
    daysForCluster = uniqueDays(assignment == i); %assignments is from "Monkey_TrainingInDotsModels"
    condition = ismember(dates, daysForCluster) & ismember(coherence, [0]);
    data{1} = (choice(condition) - 1).*2 - 1; %have to subtract -1 since the function file works with 0s and 1s
    data{2} = matME(condition, :); %DV for each trial, in normalize form
    type = 'non-RT';
    params = 2;
    tau05 = 0; %not being used but needed for function to work
    options = optimset('MaxFunEvals', 1000, 'MaxIter', 1000, 'TolFun', 1e-4, 'TolX', 1e-4);
    %
    for z = 1000:1000:10000
        params = z; %try many different param starting points
        [bestOutlier, err_] = fminsearch(@(params) SpecificOutlier(data,type,tau05,params), params, options);
        if z == 1000 || err_ < errorOutlier
            errorOutlier = err_;
            outlierModelChoice = bestOutlier;
        end
    end
    [~, fitChoices] = SpecificOutlier(data,type,tau05,outlierModelChoice); %Find the fit
    OutlierAcc(index) = sum(fitChoices == data{1})./length(fitChoices); %Percentage the model got right
    OutlierHeight(index) = outlierModelChoice;
    
    signCoherenceInt = signCoherence(condition);
    uniqSignCoh = unique(signCoherenceInt);
    avgOutlierChoice = zeros(size(uniqSignCoh));
    avgRealChoice = zeros(size(uniqSignCoh));
    for c = 1:length(uniqSignCoh)
        avgOutlierChoice(c) = nanmean(fitChoices(signCoherenceInt == uniqSignCoh(c)));
        avgRealChoice(c) = nanmean(data{1}(signCoherenceInt == uniqSignCoh(c)));
    end
 
    figure(409 + index);
    plot(uniqSignCoh, avgOutlierChoice, 'bo--', uniqSignCoh, avgRealChoice, 'ro--')
    xlabel('Coherence')
    ylabel('Prob Right')
    title('Probability of Right with Fit: Outlier')
    legend('Model', 'Empirical')
    
    index = index + 1;
end
%% ------------------------------------------------------------------- %%
% --------------------------------------------------------------------- %
% ----------- Simulation Analysis ------------------------------------- %
% It seems that the monkeys simply adjust their boundary, but now we need simulation proof. Or a simple model
%for this to work you have to go into "MotionEnergyCondorAnalysis.m" and
%load the "GenerateDotsForPCA_2.mat" file and run it through the m file.

% bestFiles = best2worst(1); %contains the best to worst files in the condor (in motionenergycondoranalysis)
% filename = sprintf('motionEnergySim%d.mat',bestFiles); %iterate through file
% load(filename); 

[matME, matDVfromME, endDVfromME] = matrixConverterME(totalME);

%Ideal observer
numIncreases = 5;
avgMESimRight = cell(1,numIncreases);
avgMESimLeft = cell(1,numIncreases);
avgMESimRightOutlier = cell(1,numIncreases);
avgMESimLeftOutlier = cell(1,numIncreases);

%Bound and Outlier Updates
updateBound = 0;
updateOutlier = 0;
incBound = 1000;
incOutlier = 1000;
for ii = 1:numIncreases
    updateBound = updateBound + incBound;
    updateOutlier = updateOutlier + incOutlier;
    for i = 1:size(matDVfromME, 1)
        trialDV = matDVfromME(i,:);

        crossUB = find(trialDV >= updateBound, 1);
        crossLB = find(trialDV <= -updateBound, 1);
        if any(crossUB) || any(crossLB)
            indexCross = min([crossUB crossLB]);
            DVchoice(i) = trialDV(indexCross);
            rtDVchoice(i) = indexCross;
        else
            temp = trialDV(~isnan(trialDV));
            DVchoice(i) = temp(end);
            rtDVchoice(i) = sum(temp);
        end
        
        %Outlier
        trialME = matME(i,:); %momentary evidence zscored 
        crossUB = find(trialME >= updateOutlier, 1);
        crossLB = find(trialME <= -updateOutlier, 1);
        if any(crossUB) || any(crossLB)
            indexCross = min([crossUB crossLB]);
            outlierChoice(i) = trialME(indexCross);
            rtOutlier(i) = indexCross;
        else
            temp = trialME(~isnan(trialME));
            outlierChoice(i) = temp(end);
            rtOutlier(i) = sum(temp);
        end    
    end
    simChoice = (sign(DVchoice) + 1)./2; %turn to 1 and 0
    simChoiceOutlier = (sign(outlierChoice)+ 1)./2;
    
    %Motion Energy for Right vs Left in Accumulation Model
    MotionEnergyforSimRight = normMatME(simChoice == 1 & coherence == 0, :); %pick normalize ME for certain days (that belong to the cluster of interest)
    MotionEnergyforSimLeft = normMatME(simChoice == 0 & coherence == 0, :);
    %select between choices that are left vs right
    avgMESimRight{ii} = nanmean(MotionEnergyforSimRight);
    avgMESimLeft{ii} = nanmean(MotionEnergyforSimLeft); 
    
    %Motion Energy for Right vs Left in Outlier Model
    MotionEnergyforSimRightOutlier = normMatME(simChoiceOutlier == 1 & coherence == 0, :); %pick normalize ME for certain days (that belong to the cluster of interest)
    MotionEnergyforSimLeftOutlier = normMatME(simChoiceOutlier == 0 & coherence == 0, :);
    %select between choices that are left vs right
    avgMESimRightOutlier{ii} = nanmean(MotionEnergyforSimRightOutlier);
    avgMESimLeftOutlier{ii} = nanmean(MotionEnergyforSimLeftOutlier); 
end

time = 0:.0083:(size(normMatME,2)-1)*.0083;
%Plot Accumulation Results
stepColor = 1/(numIncreases);
color = 0;
figure(61);
hold on;
for i = 1:numIncreases
    plot(time, avgMESimRight{i}, 'Color', [1-color, 0 ,0]); 
    plot(time, avgMESimLeft{i}, 'Color', [0, 0 ,1-color]); 
    color = color + stepColor;
end
hold off
title('Motion Energy: Right vs Left for Sim In Accumulation')

%Plot Outlier Results
stepColor = 1/(numIncreases);
color = 0;
figure(62);
hold on;
for i = 1:numIncreases
    plot(time, avgMESimRightOutlier{i}, 'Color', [color, 0 ,0]); 
    plot(time, avgMESimLeftOutlier{i}, 'Color', [0, 0 ,color]);
    color = color + stepColor;
end
hold off
title('Motion Energy: Right vs Left for Sim In Outlier')
%% Also model what the outcome would look like if the sensitivity of animal change over time. WHile keeping bound fix
% One way to do this is assume the ME is the highest sensitivity (k) that
% can be accepted and if you multiply the samples by ratios between 0 to 1,
% with increasing values it would show monkey getting more sensitive over
% time. And graph those psychometric answers.

bestFiles = best2worst(1); %contains the best to worst files in the condor (in motionenergycondoranalysis)
filename = sprintf('motionEnergySim%d.mat',bestFiles); %iterate through file
load(filename); 

[matME, matDVfromME, endDVfromME] = matrixConverterME(totalME);
normMatME = ones(size(matME));
perceivedMatME = ones(size(matME));
for i = 1:length(signUniqueCoh) %scaling not normalizing (with K we are altering mean while keeping variance the same)
        tempME = matME(coherence.*cos(deg2rad(direction)) == signUniqueCoh(i),:); %pick ME that matter based on coherence
        avgME = nanmean(tempME); %average ME
        stdME = nanstd(tempME); %standard dev
        lengthME = sum(~isnan(tempME)); %length
        steME = stdME./sqrt(lengthME); %standard error
        normMatME(coherence.*cos(deg2rad(direction)) == signUniqueCoh(i),:) = (tempME - 0)./steME; %normalize, placing new values into original position
end
    
%Ideal observer
k = .1:.1:1;
avgMESimRight = cell(1,length(k));
avgMESimLeft = cell(1,length(k));
%Bound and Outlier Updates
bound = 1000;
for ii = 1:length(k)
    for i = 1:length(signUniqueCoh) %scaling not normalizing (with K we are altering mean while keeping variance the same)
        tempME = matME(coherence.*cos(deg2rad(direction)) == signUniqueCoh(i),:); %pick ME that matter based on coherence
        avgME = nanmean(tempME); %average ME
        stdME = nanstd(tempME); %standard dev
        lengthME = sum(~isnan(tempME)); %length
        steME = stdME./sqrt(lengthME); %standard error
        perceivedMatME(coherence.*cos(deg2rad(direction)) == signUniqueCoh(i),:) = (tempME - avgME*(1-k(ii)))./steME; %normalize, placing new values into original position
    end
    
    DV = cumsum(perceivedMatME , 2); %find pseudo DV that updates based on sensitivty (k) changes
    for i = 1:size(DV, 1)
        trialDV = DV(i,:);

        crossUB = find(trialDV >= bound, 1);
        crossLB = find(trialDV <= -bound, 1);
        if any(crossUB) || any(crossLB)
            indexCross = min([crossUB crossLB]);
            DVchoice(i) = trialDV(indexCross);
            rtDVchoice(i) = indexCross;
        end   
    end
    simChoice = (sign(DVchoice) + 1)./2; %turn to 1 and 0
    
    %Motion Energy for Right vs Left in Accumulation Model
    MotionEnergyforSimRight = normMatME(simChoice == 1 & coherence == .512, :); %pick normalize ME for certain days (that belong to the cluster of interest)
    MotionEnergyforSimLeft = normMatME(simChoice == 0 & coherence == .512, :);
    %select between choices that are left vs right
    avgMESimRight{ii} = nanmean(MotionEnergyforSimRight);
    avgMESimLeft{ii} = nanmean(MotionEnergyforSimLeft);  
end

time = 0:.0083:(size(normMatME,2)-1)*.0083;
%Plot Accumulation Results
stepColor = 1/(length(k));
color = 0;
figure(61);
hold on;
for i = 1:length(k)
    plot(time, avgMESimRight{i}, 'Color', [color, 0 ,0]); 
    plot(time, avgMESimLeft{i}, 'Color', [0, 0 ,color]);
    color = color + stepColor;
end
hold off
title('Motion Energy: Right vs Left for Sim In Accumulation')