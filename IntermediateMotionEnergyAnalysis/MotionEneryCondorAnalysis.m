% Miguel VIvar-Lazo
% Fetsch Lab
% 02/10/ 2020
% Analyzing DV End behavior of ME

load GenerateDotsForPCA_2.mat
%load GenjiTrainData_Condor_ForME.mat
organize = 0;
if ~organize
    trialCoh = NaN(1,length(c)); %vector for coherence
    direction = NaN(1,length(c)); %vector for directions
    for i = 1:length(c)
        trialCoh(i) = c{i}.stimulus.coherence; %log the trial
        direction(i) = c{i}.stimulus.direction; %log the direction
    end
    coherence = abs(trialCoh);
    clear dir %dir is used as a variable in previous file but we want to use it as its function 
end

%% After loading file go through condor outputs to find the best K and Speed value 
a=dir(['C:\Users\mviva\Documents\MATLAB\Fetsch Lab\Analyzing Motion Energy\motionEnergyCondor\motionEnergySimCondor', '\*.mat']);
%a=dir(['C:\Users\mviva\Documents\MATLAB\Fetsch Lab\Analyzing Motion Energy\motionEnergyCondor\motionEnergyGenjiCondor', '\*.mat']);
out=length(a); %length of the MAT files in the directory 

signCoh = coherence.*cos(deg2rad(direction));
signUniqueCoh = unique(signCoh); %unique coherences
params = 100;
data_ = NaN(4,length(signUniqueCoh), params); %contruct data file, rows represent the variables of interest (K, Speed, Avg Correct, Avg Right) column coherences they belong to, and slices the  

for i = 1:params
    if  isfile(sprintf('motionEnergySim%d.mat',i)) %sprintf('CondorMEGenji%d.mat',i)
        filename = sprintf('motionEnergySim%d.mat',i); %sprintf('CondorMEGenji%d.mat',i); %iterate through file 
        load(filename); %load file
        %Log speed, k, and error 
        [~, ~, endDVfromME] = matrixConverterME(totalME); %Convert ME information in matrix, totalME is the name of the variable from 'filename' we are interested in
        percCorrect = zeros(size(signUniqueCoh)); %Construct vector for percent correct
        probRight = zeros(size(signUniqueCoh)); %construct vector for percent to right
        for ii = 1:length(signUniqueCoh)
            right_ = (sign(endDVfromME(signCoh == signUniqueCoh(ii))) == 1); %choices for the right in the coherence we care about
            probRight(ii) = mean(right_); %Average to the right
        end
        alike = cos(deg2rad(direction)) == sign(endDVfromME); %ask if the DV end matches the direction of the stimulus
        percCorrect = mean(alike); %percentage that DV end matches coherence
        data_(1,:,i) = k;
        data_(2,:,i) = speed;
        data_(3,:,i) = percCorrect;
        data_(4,:,i) = probRight;    
    end
end

pCorrect = zeros(size(params));
figure(24);
hold on;
for i = 1:size(data_,3)
    plot(signUniqueCoh, data_(4,:,i)); %plot all 'prob to the right' stats
    pCorrect(i) = data_(3,1,i); %log percent correct for each parameter 
end
hold off;
xlabel('Coherence')
ylabel('Probability Right')
title('All Motion Energy Parameteres')

pCorrect(isnan(pCorrect)) = 0;
[~, best2worst] = sort(pCorrect,'descend'); %which is the most correct
%% Analyze best results (grab best 5 results)
chk = 5; %How many outcomes should i check?
avgMEbyCoh = zeros(chk,length(signUniqueCoh));
varMEbyCoh = zeros(chk,length(signUniqueCoh));
for i = 1:chk
    best = best2worst(i);
    bestRightAvg = data_(4,:,best); %find the best outcome
    figure(9);
    plot(signUniqueCoh, bestRightAvg); hold on; %plot it

    % Lets see how the best ME maps to coherence with Variance and Mean Relationship 
    filename = sprintf('motionEnergySim%d.mat',best);  %sprintf('CondorMEGenji%d.mat',best); %iterate through file
    load(filename); %load file
    [matME, matDVfromME, endDVfromME] = matrixConverterME(totalME); %Convert ME information in matrix
    for ii = 1:length(signUniqueCoh)
        temp = matME(signCoh == signUniqueCoh(ii),:);
        avgME = nanmean(temp(:));
        avgMEbyCoh(i,ii) = avgME; %find ME average

        varME = nanvar(temp(:));
        varMEbyCoh(i,ii) = varME; %find ME variance
    end
end
hold off;
xlabel('coherence')
ylabel('Prob')
title('probability of right')
legend('K=110, S=3','K=120 S=3','K=90, S=3','K=80, S=3','K=70, S=3')

bestSpeed = squeeze(data_(2,1,best2worst(1:chk))); %keep the speed and K valuies
bestK = squeeze(data_(1,1,best2worst(1:chk)));

figure(25);
plot(signUniqueCoh, avgMEbyCoh) %plot average ME based on coherence
xlabel('Coherence')
ylabel('Motion Energy (a.u)')
title('Average ME Relationship with Coherence')
axis tight
legend('K=110, S=3','K=120 S=3','K=90, S=3','K=80, S=3','K=70, S=3')
figure(27)
plot(signUniqueCoh, varMEbyCoh) %plot variance ME based on coherence
xlabel('Coherence')
ylabel('Motion Energy (a.u)')
title('Variance ME Relationship with Coherence')
legend('K=110, S=3','K=120 S=3','K=90, S=3','K=80, S=3','K=70, S=3')
