%Miguel Vivar-Lazo
%05/19/2020
%Fetsch Lab
%Analysis of Motion Filters
twoDSpatial = 0;

% Sptial Filter
speed = 6;
sigmult = speed*0.125+0.204;
sigma_c = 0.35 * sigmult; % controls period ...? ; mvl: or controls the width of the curve
sigma_g = 1 * sigmult;%0.05; % controls width of Gaussian envelope along orthogonal dir, 
%Even and Odd Symmetric Fourth-Order Cauchy 
if twoDSpatial == 0
    f1 = @(x) cos(atan(x/sigma_c)).^4 .* cos(4*atan(x/sigma_c)); % even
    f2 = @(x) cos(atan(x/sigma_c)).^4 .* sin(4*atan(x/sigma_c)); % odd
    % Temporal Filter
    n1=3; n2=5;
    k = 60; % Normally modulated by Speed of the dots (Degrees/second)
    % Low Pass Filter
    g1 = @(t) (k*t).^n1 .* exp(-k*t) .* (1/factorial(n1) - (k*t).^2/factorial(n1+2)); % fast
    g2 = @(t) (k*t).^n2 .* exp(-k*t) .* (1/factorial(n2) - (k*t).^2/factorial(n2+2)); % slow
else %This doesnt really work, stay indicator at 0
    f1 = @(x,y) cos(atan(x/sigma_c)).^4 .* cos(4*atan(x/sigma_c)) .* exp(-y.^2/(2*sigma_g^2)); % even
    f2 = @(x,y) cos(atan(x/sigma_c)).^4 .* sin(4*atan(x/sigma_c)) .* exp(-y.^2/(2*sigma_g^2)); % odd
end
%% Plot Filters
if twoDSpatial == 0
    xCord = linspace(-floor(6/2)/2, floor(6/2)/2,235); %235: Is the number of pixels in diameter %linspace(-1,1,101); %Space Axis
    tCord = 0:.0083:1.4940; %linspace(0,1,101); %Time axis
    evenf = f1(xCord); %Even
    oddf = f2(xCord); %Odd
    eveng = g1(tCord); %Even
    oddg = g2(tCord); %Odd
else %If you forget look at the spatial filter and notice that x direction will descrive orientation and y-dir only descrives length of the filter in the y direction 
    xCord = linspace(-1,1,101); %Space Axis
    yCord = linspace(-1,1,101); %Space Axis
    [x_, y_] = meshgrid(xCord, yCord);
    evenf = f1(x_, y_); %Even
    oddf = f2(x_, y_); %Odd
    figure; imagesc(evenf) %Even 2D spatial filter
    figure; imagesc(oddf) %odd 2D spatial filter
    figure;plot(exp(-yCord.^2/(2*1^2)))
    figure;plot(exp(-yCord.^2/(2*.05^2)))
end
figure(1);
subplot(2,2,1);
plot(xCord, evenf)
title('Even Spatial Filter')
subplot(2,2,2);
plot(xCord, oddf)
title('Odd Spatial Filter')
subplot(2,2,3);
plot(tCord, eveng)
title('Even Temporal Filter')
subplot(2,2,4);
plot(tCord, oddg)
title('Odd Spatial Filter')
%% Combine filters using outer dot product -> Plot
stfEE = evenf' * eveng; %Phase = 0
stfEO = evenf' * oddg; %Phase = 0
stfOE = oddf' * eveng; %Phase = 180
stfOO = oddf' * oddg; %Phase = 180

figure(2);
subplot(2,2,1);
imagesc(stfEE');
%colormap('gray')
colorbar
xlabel('Position')
ylabel('Time')
title('SpatialTemporal Filter: Even Fast')
subplot(2,2,2);
imagesc(stfEO');
%colormap('gray')
colorbar
xlabel('Position')
ylabel('Time')
title('SpatialTemporal Filter: Even Slow')
subplot(2,2,3);
imagesc(stfOE');
%colormap('gray')
colorbar
xlabel('Position')
ylabel('Time')
title('SpatialTemporal Filter: Odd Fast')
subplot(2,2,4);
imagesc(stfOO');
%colormap('gray')
colorbar
xlabel('Position')
ylabel('Time')
title('SpatialTemporal Filter: Odd Slow')
%% Frequency Domain of Spatial Filters
freqDomainEE = fftshift(abs(fft2(stfEE')));
freqDomainEO = fftshift(abs(fft2(stfEO')));
freqDomainOE = fftshift(abs(fft2(stfOE')));
freqDomainOO = fftshift(abs(fft2(stfOO')));

figure(11);
subplot(2,2,1);
imagesc(freqDomainEE);
colormap('gray')
xlabel('Position')
ylabel('Time')
title('Frequency Even Fast')
subplot(2,2,2);
imagesc(freqDomainEO);
colormap('gray')
xlabel('Position')
ylabel('Time')
title('Frequency Even Slow')
subplot(2,2,3);
imagesc(freqDomainOE);
colormap('gray')
xlabel('Position')
ylabel('Time')
title('Frequency Odd Fast')
subplot(2,2,4);
imagesc(freqDomainOO);
colormap('gray')
xlabel('Position')
ylabel('Time')
title('Frequency Odd Slow')
%% Before Adding Filters, see what happens to simple sine (even) cos(odd) waves when you add them.
% You should notice a similar property ass adding the odd and even filters.
% In that amplitudes double for the favoring direction.
angles = 0:359;
figure(10); 
plot(angles, cosd(angles)); 
hold on; 
plot(angles, sind(angles)); 
plot(angles, cosd(angles)+sind(angles)); %The difference in phase between cos+sin & cos-sin is still 90 degrees, why?
plot(angles, cosd(angles)-sind(angles)); 
plot(angles, cosd(angles).^2+sind(angles).^2, '--');
plot(angles, (cosd(angles)+sind(angles)).^2 + (cosd(angles)-sind(angles)).^2, '--');
hold off;

legend('Cos','Sine','Cos + Sin','Cos - Sin', 'Square', 'Square After + & -')

%% Add the Spatiotemporal Filters to get Orientated Filters
RightEven = stfEE + stfOO; 
RightOdd = stfEO - stfOE; %Odd is Slow
LeftEven = stfEE - stfOO; 
LeftOdd = stfEO + stfOE;

%Every other sort of combination seems to give similar filters as the
%spatialtemporal ones, in order to create one that is different and unique
%for L or Right motion it needs both spatial and temporal phase shift, if
%only one phase shift exist the outcome appears as a shift in the one of
%the spatiotemporal filers.
% StationaryEven = stfEO + stfOO;
% StationaryOdd = stfEO - stfOO;

figure(3);
subplot(2,2,1);
imagesc(RightEven');
%colormap('gray')
colorbar
xlabel('Position')
ylabel('Time')
title('Orientated Right Filter: Even')
subplot(2,2,2);
imagesc(RightOdd');
%colormap('gray')
colorbar
xlabel('Position')
ylabel('Time')
title('Orientated Right Filter: Odd')
subplot(2,2,3);
imagesc(LeftEven');
%colormap('gray')
colorbar
xlabel('Position')
ylabel('Time')
title('Orientated Left Filter: Even')
subplot(2,2,4);
imagesc(LeftOdd');
%colormap('gray')
colorbar
xlabel('Position')
ylabel('Time')
title('Orientated Left Filter: Odd')
%% Orientated Filters in Frequency Domain
%fftshift(log(1+abs(F)))
freqDomainRE = fftshift(abs(fft2(RightEven')));
freqDomainRO = fftshift(abs(fft2(RightOdd')));
freqDomainLE = fftshift(abs(fft2(LeftEven')));
freqDomainLO = fftshift(abs(fft2(LeftOdd')));

time = linspace(-60,60,180); %180 is the number of frames in time, and 60 is .5 of 120 hz
space = 1:235;%linspace(-60,60,180);
figure(4);
subplot(2,2,1);
imagesc(space, time, freqDomainRE);
colormap('gray')
xlabel('Position')
ylabel('Time')
title('Frequency Orientated Right Filter: Even')
subplot(2,2,2);
imagesc(space, time, freqDomainRO);
colormap('gray')
xlabel('Position')
ylabel('Time')
title('Frequency Orientated Right Filter: Odd')
subplot(2,2,3);
imagesc(space, time, freqDomainLE);
colormap('gray')
xlabel('Position')
ylabel('Time')
title('Frequency Orientated Left Filter: Even')
subplot(2,2,4);
imagesc(space, time, freqDomainLO);
colormap('gray')
xlabel('Position')
ylabel('Time')
title('Frequency Orientated Left Filter: Odd')
%% Quadrature Pairing
% Right = RightEven.^2 + RightOdd.^2; %Out of phase, Even and Odd -> should convert to energy
% Left = LeftEven.^2 + LeftOdd.^2; %And be phase independent
% 
% RightFreqD = fftshift(abs(fft2(Right')));
% LeftFreqD = fftshift(abs(fft2(Left')));
% 
% figure(5);
% subplot(2,2,1);
% imagesc(Right');
% colormap('gray')
% xlabel('Position')
% ylabel('Time')
% title('Right Energy Filter')
% subplot(2,2,2);
% imagesc(Left');
% colormap('gray')
% xlabel('Position')
% ylabel('Time')
% title('Left Energy Filter')
% subplot(2,2,3);
% imagesc(RightFreqD);
% colormap('gray')
% xlabel('Position')
% ylabel('Time')
% title('Frequency Right Energy Filter')
% subplot(2,2,4);
% imagesc(LeftFreqD);
% colormap('gray')
% xlabel('Position')
% ylabel('Time')
% title('Frequency Left Energy Filter')
%% Now analyze the Random Dots in Frequency Domain
% One trick is to only draw a matrix for X direction and T time, that way
% simplifying the 3D space to the two components that really matter.
%load('GeneratedDotPositionData_withMEProcess_V2.mat')
% First convert dot positions to matrix
 nTrials = length(save_struct); %length(best_data.dotPos); 
 MotionMat = cell(1,nTrials);
for c = 1:nTrials
    diameterPPD = ceil(save_struct(c).dots_struct.aperture(3) * screen_struct.pix_per_deg); %  ceil(best_data.aperture(3) * best_data.ppd); %diameter for pixels per degrees
    slice = 3; %Number of slices in a frame, I remember this but it might be important to add to 'best_data'
    maxFrames = 0;
    frames =  size(save_struct(c).dots_struct.dot_pos, 3); %size(best_data.dotPos{c}, 3);
    MotionMat{c} = zeros(frames, diameterPPD);
    maxFrames = max([maxFrames, frames]);
    for i = 1:frames
        xCord = save_struct(c).dots_struct.dot_pos(:,1,i); % best_data.dotPos{c}(:,1,i); %Dots Position in X Direction
        xPos = xCord + ceil(diameterPPD/2); %Shift by the correct radius, since the values frmo dorPos are center around 0 and we can work with negative positions in matlab
        xPos(isnan(xPos)) = []; %Elimate nans
        MotionMat{c}(i, xPos) = 1; %add dots (1s) to matrix, for which frame I am in
    end
end
%% Now simply convert the images to frequency domain and see how they look
freqMotionMat_Danique = cell(1,length(MotionMat));
for c = 1:length(MotionMat)
    freqMotionMat_Danique{c} = fftshift(abs(fft2(MotionMat{c})));
end

% time = linspace(-60,60,180); %180 is the number of frames in time, and 60 is .5 of 120 hz
% space = 1:235;%linspace(-60,60,180);
figure;imagesc(freqMotionMat_Danique{1260}); colormap('gray'); %imagesc(space, time, freqMotionMat{21999}); colormap('gray')
%% Now average over all Frequency Domains, For Left Stim and Right Stim
leftFreqMotionMat = zeros(maxFrames, diameterPPD);
NleftFreqMotionMat = zeros(maxFrames, diameterPPD);
rightFreqMotionMat = zeros(maxFrames, diameterPPD);
NrightFreqMotionMat = zeros(maxFrames, diameterPPD);
%Only go for the the high Coherence to get a good estimate of frequency
%component
for c = 1:length(freqMotionMat)
    if direction(c) == 180
        leftFreqMotionMat(1:size(freqMotionMat{c},1), 1:size(freqMotionMat{c},2)) = leftFreqMotionMat(1:size(freqMotionMat{c},1), 1:size(freqMotionMat{c},2)) + freqMotionMat{c}; %Averaging is going to be tricky
        NleftFreqMotionMat(1:size(freqMotionMat{c},1), 1:size(freqMotionMat{c},2)) = NleftFreqMotionMat(1:size(freqMotionMat{c},1), 1:size(freqMotionMat{c},2)) + 1;
        %Construct matrix with that counts the total number of trials that
        %belong to a pixal
    elseif direction(c) == 0
        rightFreqMotionMat(1:size(freqMotionMat{c},1), 1:size(freqMotionMat{c},2)) = rightFreqMotionMat(1:size(freqMotionMat{c},1), 1:size(freqMotionMat{c},2)) + freqMotionMat{c};
        NrightFreqMotionMat(1:size(freqMotionMat{c},1), 1:size(freqMotionMat{c},2)) = NrightFreqMotionMat(1:size(freqMotionMat{c},1), 1:size(freqMotionMat{c},2)) + 1;
    end
end

figure(20);
imagesc(leftFreqMotionMat(1:25,:)./NleftFreqMotionMat(1:25,:));
colormap('gray')
colorbar
figure(21);
imagesc(rightFreqMotionMat(1:25,:)./NrightFreqMotionMat(1:25,:));
colormap('gray')
colorbar



