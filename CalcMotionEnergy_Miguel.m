%MVL: Motion Energy -> Its going to be better... at least nicer (hopefully)

% function totalME = ME_par_simple_forPldaps(data,plotflag,saveDotPositions)
% Takes in a data structure containing aperture[X, Y, diam], 
% rseed, coh, speed, dir, dur, and computes motion energy.
%   -CRF beginning 09/2013
clear
% example file:
load('Users\mviva\Documents\MATLAB\Fetsch Lab\BasicSaccadeData\Hanzo\hanzo20181221DotsBasic_setup1207.PDS','-mat');

plotflag = 0;
anti_aliasing = 0;

% constants
data.ppd            = PDS.initialParameters{4}.display.ppd; % pixels per deg
data.screenRect     = PDS.initialParameters{4}.display.winRect;
data.refresh        = PDS.initialParameters{4}.display.frate;

data.aperture       = PDS.initialParameters{3}.stimulus.aperture(1:3); % x,y,diam
data.speed          = PDS.initialParameters{3}.stimulus.speed; % deg/s
data.dotSize        = PDS.initialParameters{3}.stimulus.dotDiam; % pixels (n x n squares, usually)
data.prefDir        = min(PDS.initialParameters{3}.stimulus.direction); % defines 'positive' motion energy (arbitrary for behavior-only blocks)

screenCenter        = PDS.initialParameters{4}.display.ctr;
apertureDiameter    = data.aperture(3); % diam of aperture in degrees
diameterPPD         = apertureDiameter * data.ppd;	% diam of aperture in pixels
apertureCent        = [screenCenter(1) + data.aperture(1)*data.ppd, screenCenter(2) - data.aperture(2)*data.ppd];  % aperture center in screen coords
left                = apertureCent(1) - diameterPPD/2; %Aperture length from center to left
top                 = apertureCent(2) - diameterPPD/2; %Aperture length from center to top
right               = apertureCent(1) + diameterPPD/2; %Aperture length from center to right
bottom              = apertureCent(2) + diameterPPD/2; %Aperture length from center to bottom
apertureRect        = [left top right bottom]; % aperture in 'rect' format, in Screen (psychtoolbox) pixel values
dotSize             = data.dotSize; %Size of each dot drawn

data.coherence = nan(size(PDS.data))';
data.direction = nan(size(PDS.data))';
data.duration = nan(size(PDS.data))';
totalME = cell(size(PDS.data))';

% can try to use parfor eventually, to speed this up
for n = 1:length(PDS.data) % trials    
    data.coherence(n) = PDS.conditions{n}.stimulus.coherence;
    data.direction(n) = PDS.conditions{n}.stimulus.direction;
    data.duration(n) = PDS.data{n}.stimulus.dotDuration;

    dotPos = PDS.data{n}.stimulus.dotPos;
    nframes = size(dotPos,3); %Extract number of frames
    % old ME code assumes dotXY is organized as [nframes, ndots, xy], but
    % PLDAPS generateDots saves it as [ndots, xy, nframes], so:
    %dotPos = permute(dotPos,[3 1 2]);
        
    % make a 'pixel matrix' of 1s and 0s covering the aperture for each frame
    pixelMat = zeros(ceil(diameterPPD), ceil(diameterPPD), nframes);
       
    %Allocate space for X and Y positions of Dots 
    for f = 1:nframes
        %We round up because it looks like Psychtoolbox rounds up for
        %position of Circle (Square actually)
        xPix = ceil(dotPos(:,1,f) + apertureCent(1)); %Calculate each position of dots in aperture , iterate through frames
        xPix = xPix(~isnan(xPix)); %Remove nans
        %Remember that for Y, psychtoolbox is inverted, so 0 is top 
        %Also we round down because PsychtoolBox seems to do this
        yPix = floor(dotPos(:,2,f) + apertureCent(2)); %Calculate each position of dots in aperture , iterate through frames
        yPix = yPix(~isnan(yPix)); %Remove nans
        
        %Change Coordinate System of DotPos from Whole Screen to Pos
        %relative to Aperture
        tempXPix = ceil(xPix - apertureRect(1)); %Subtract Position of Dots by Left Side of Aperture
        tempYPix = ceil(yPix - apertureRect(2)); %Subtract Position of Dots by Right side of Aperture 
        tempYPix = ceil(diameterPPD) - tempYPix; %Flip the Y axis 
        
        %Calculate Radius and Circumference Values of Draw Dots and the
        %values in between Circumference and Center
        dotRadius = (dotSize-1)/2; %Find the radius of the DrawDot with respect to the center dot
        degrees = 1:1:360; %Number of degrees we care for 'draw dot' area
        xLim = (1:.5:dotRadius)'.*cos(degrees.*(pi/180)); %All values of a Circle, with center 0 and Radius belonging to DrawDots, in the circumference and between Circumference and Cetner
        xLim = reshape(xLim',1, size(xLim,1)*size(xLim,2)); %Have to transpose xLim' so reshape construct the correct vector, Vector with all Values between Center and Circumference
        yLim = (1:.5:dotRadius)'.*sin(degrees.*(pi/180)); %Same as X
        yLim = reshape(yLim',1, size(yLim,1)*size(yLim,2)); %Same as Y
        
        indexInDrawDotsX_ = tempXPix + xLim.*ones(length(tempXPix),1); %Values of Pixels within the DrawDots Circumference: X
        indexInDrawDotsX  = reshape(indexInDrawDotsX_',1, size(indexInDrawDotsX_,1)*size(indexInDrawDotsX_,2)); %Make into vector
        indexInDrawDotsY_ = tempYPix + yLim.*ones(length(tempYPix),1); %Values of Pixels within the DrawDots Circumference: Y
        indexInDrawDotsY  = reshape(indexInDrawDotsY_',1, size(indexInDrawDotsY_,1)*size(indexInDrawDotsY_,2)); %Make into vector
        
        %Eliminate any values outside aperuture -> Occures becuase of the
        %way I add Indexes Values onto TempXPix and TempYPix
        indexInDrawDotsXFinal = indexInDrawDotsX(indexInDrawDotsX > 0 & indexInDrawDotsX <= ceil(diameterPPD) & ...
            indexInDrawDotsY > 0 & indexInDrawDotsY <= ceil(diameterPPD)); 
        indexInDrawDotsYFinal = indexInDrawDotsY(indexInDrawDotsY > 0 & indexInDrawDotsY <= ceil(diameterPPD) & ...
            indexInDrawDotsX > 0 & indexInDrawDotsX <= ceil(diameterPPD));
        
        %Convert X-Y indexes into Single index
        tempIndex = sub2ind(size(pixelMat(:,:,f)), round(indexInDrawDotsXFinal), round(indexInDrawDotsYFinal)); %All Single value index
        tempIndex = unique(tempIndex); % Eliminate repeats constructed by tight Circle intervals
        tempMat = pixelMat(:, :, f); %For Holding Matrix purposes
        tempMat(tempIndex) = 1; %Insert Ones for dot Positions into pixelMat
        pixelMat(:, :, f) = tempMat; %Permanantly safe it
        
        %Here we will just concentrate in plotting correctly the
        %circumference intensity values, to anti-alias correctly
        %We are going to use the Anti-Aliasing circle technique from Xiaolin Wu
        % 1991, "An efficient Antialiasing Technique". best explain by wiki
        %Since technique is normally use for 1D but cirlce is 2D we need to
        %combine out x and y coordinates in x^2 + y^2 fashion to get a
        %radius value we can calcualte the intensities from
%         if anti_aliasing %This might be too taxing and consuming to do,
%         look into using a filte instead
%             xCircumference = indexInDrawDotsX_(:,end); %Circumference X coordinate
%             yCircumference = indexInDrawDotsY_(:,end); %Circumference Y coordinate
%             circumferenceValues     = xCircumference.^2 + yCircumference.^2; %MVL: Anti-aliasing technique (Keeping the residual of the circumference line)
%             circumferenceResidues   = floor(circumferenceResidual) - circumferenceResidual; %Residual from floor of its value
%             
%             %First two loops pick out the points
%             for i = xCircumference
%                 for ii = yCircumference
%                     if xCircumference
%                     end
%                 end
%             end
%             
%         end
             
        if plotflag
            figure(1); % xlim([1 data.screen_rect(3)]);ylim([1 data.screen_rect(4)]); 
            hold on;
            contour(squeeze(pixelMat(:,:,f))');
            set(gca);
            drawnow
        end
    end
    
    % Construct the filters for motion energy
    
    %MVL: Sigmult is being altered by the speed because the speed on which
    %a stimuli moves. This is because we would like to make a filter that
    %will create non-aliasing images. The way in which this is modulated
    %(line of code below) is still a mystery to me
    
    %MVL: Since Chris stated that the code seems to be off, the place to
    %probably start is at here and fix the Sigmult in order so Coherence
    %and motion energy have a linear relationship
    sigmult = data.speed*0.125+0.204; % from ad hoc linear regression 
    
    % first adjust filters for speed
    % k affects the temporal filters that we now know could be altered
    % because of the higher sampling frequency, but I forgot if we increase
    % it or decrease it. 
    if data.speed<4
        k = 45; %MVL: What is K?
    elseif data.speed>13
        k=75; 
    else
        k=60; 
    end

    % SPATIAL FILTERS f1,f2 (see Adelson & Bergen 1985, Kiani et al. 2008, Bollimunta et al. 2012)
    sigma_c = 0.35 * sigmult; % controls period ...?
    sigma_g = 0.05 * sigmult; % controls width of Gaussian envelope along orthogonal dir
    %MVL: It seems that the Spatial Filters are of 2D nature, which would
    %make sense. Also, for both f1 and f2 if you exclude the gaussian
    %(exp(-y..)) you get the gabora filters from the paper. Does the
    %gaussian addition make it a 3D Gabor? Why wouldnt you need another 1D
    %gabor?
    f1 = @(x,y) cos(atan(x/sigma_c)).^4 .* cos(4*atan(x/sigma_c)) .* exp(-y.^2/(2*sigma_g^2)); % even
    f2 = @(x,y) cos(atan(x/sigma_c)).^4 .* sin(4*atan(x/sigma_c)) .* exp(-y.^2/(2*sigma_g^2)); % odd
    
    %MVL: So interestingly enough, X is equal to xp and Y is equal to yp
    %(this is after the meshgrid modification)... Also why divided by 20,
    %arbitrary I am assuming. 
    %(dont divide by 20) only 2
    x = linspace(-apertureDiameter/20,apertureDiameter/20,ceil(d_ppd)); % need same size as aperture for Fourier method
    
    %MVL: x and y are just straight lines
    y = x;
    [x,y] = meshgrid(x,y);
    
    %MVL: Oh I see, so what we are doing here is simply rotating the axis to the
    %right or left, which would make sense on why much has not change to
    %the images above
    % rotate to pref/null axis
    %MVL: It rotates the values of X and Y so that its in the direction of
    %motion, this is not the filter but simply the values that get
    %implemented in the filter. 
    xp = x*cos(data.prefDir*pi/180) - y*sin(data.prefDir*pi/180);
    yp = x*sin(data.prefDir*pi/180) + y*cos(data.prefDir*pi/180);
                                 %^ this was a typo in Bollimunta
        
    % TEMPORAL FILTERS g1,g2 (see Kiani et al. 2008, Adelson & Bergen 1985)
    %k = 60; % 60 is standard, but instead set above, depending on speed
    
    %MVL: This is the factorial equation we saw in the Adelson paper
    %(parameters from paper as well)
    n1=3; n2=5;
    g1 = @(t) (k*t).^n1 .* exp(-k*t) .* (1/factorial(n1) - (k*t).^2/factorial(n1+2)); % fast
    g2 = @(t) (k*t).^n2 .* exp(-k*t) .* (1/factorial(n2) - (k*t).^2/factorial(n2+2)); % slow
    added_blankframe_no = 30;   % pad with zeros, for least the length of the temporal filter
    
    %MVL: Also a straight line (diagnal), on the right side it is the
    %amount of time that went by (using seconds per frame) in the trial's
    %dots (nframes) plus 'added_blankframe_no' minus 1. Why the extra
    %components are added, is currently a mystery. --> I know where the
    %minus one comes from, to calcuate the correct number of partitions to
    %the vector t. (If no minus one (- 1) you get a vector t of length
    %nframes+added)blankframe_no + 1)
    %MVL: I really think the addition of zeros is for the FFT or
    %interpolation purposes
    t = 0 : 1/data.refresh : (nframes+added_blankframe_no-1)*1/data.refresh;
    
    %MVL: Before we move on, interestingly spatial filters because a 2D
    %filter or were transformed in 2D filters while the temporal filters
    %was kept as a 1D filter (filtered? might not be the right word)
    
        % now convolve with dot movie (pixelMat):
    
    %MVL: Convolution time 
    % first expand dimensions to allow pointwise multiplication in x,y,t
    A = f1(xp,yp); %MVL: Adding the values we calculated (xp and xy) onto the filter functions to make
    B = f2(xp,yp); %2D filters of the spatial filters from 'motionTutorial.m'
    A(abs(A)<1e-6) = 0; 
    B(abs(B)<1e-6) = 0; % threshold to avoid floating point weirdness
    
    %MVL: Grabbing the new constructed filters and extending them through
    %time (t that we we constructed), we transpose because of the meshgrid
    %function earlier), want to use the same filter for each frame
    f1_t = repmat(A',[1,1,length(t)]); % f1(xp,yp) extended in time
    f2_t = repmat(B',[1,1,length(t)]); % f2(xp,yp) extended in time
  % (transpose here^ because functions evaluated on a meshgrid give x,y as col,row)
    %MVL: I think the transpose is for the later 'outer factor' that when
    %done in 1D filters is A' .* B in order to make Matrix.
    
    %MVL: Here they construct first two 1D filters (g1(t) and g2(t)) and
    %then turn that into a 1D filter in g1_xy and g2_xy that has a filter
    %for each time frame. 
    g1_xy = shiftdim(repmat(g1(t),[size(x,1) 1 size(x,2)]),2); % g1(t) extended in x and y
    g2_xy = shiftdim(repmat(g2(t),[size(x,1) 1 size(x,2)]),2); % g2(t) extended in x and y
    %MVL: the two above lines reform the time filter as a vector for the
    %Number of frames per Dot stimulus (3D component), each vector in each
    %dimension contains the same number for that Dimensions frame Time
    %Value.
    
    %MVL: Up to this moment they simply have been creating Filters in the
    %2D space (the temporal filters are a Vector which dimensions, so it
    %could be consider 1D). Now they will implement the layout set out in
    %the Delson papaer
    
    % note: switching to terminology from Eero Simoncelli's motionTutorial.m
    even_fast = f1_t.*g1_xy;
    odd_slow = f2_t.*g2_xy;
    odd_fast = f2_t.*g1_xy;
    even_slow = f1_t.*g2_xy;
    
    left1 =  odd_fast+even_slow; 
    left2 =  -odd_slow+even_fast;
    right1 = -odd_fast+even_slow;
    right2 = odd_slow+even_fast;
        
    %MVL: Things get confusing here again -> simple actually fft for
    %simplicity in convoluting and add zeros to keep it consistent when
    %adding zeros to 't'
    % fourier transform and multiply (courtesy NaYoung So)
    padZeros = zeros(size(pixelMat,1),size(pixelMat,2),added_blankframe_no);
    pixelMatN = cat(3,padZeros,pixelMat); %MVL: Adds the zero matrix onto PixelMat-> to the front of the 3D matrix, keep 'cat' function in mind for later
    fright1 = fftn(right1); 
    fright2 = fftn(right2);
    fleft1  = fftn(left1);  
    fleft2  = fftn(left2);
    fstim = fftn(pixelMatN);
    
    %MVL: This all makes sense. ACtually pretty intuitive and I'm dumb for
    %not thinking about it earlier. No need to abs if youre going to
    %square, thats what the Energy transformation is for
    right1Resp = abs(ifftn( fstim .* fright1 ));
    right2Resp = abs(ifftn( fstim .* fright2 ));
    left1Resp  = abs(ifftn( fstim .* fleft1  ));
    left2Resp  = abs(ifftn( fstim .* fleft2  ));
    response = (right1Resp.^2 + right2Resp.^2) - (left1Resp.^2 + left2Resp.^2);

    if plotflag
        figure(21);
        cmax = max([max(max(max(response))) abs(min(min(min(response))))]);
        for q = 1:size(response,3)
            clf; 
            contourf(response(:,:,q)); 
            colorbar;
            title(num2str(q)); 
            caxis([-cmax cmax]);
            pause(0.05)
        end    
    end
    
    totalMotion = squeeze(sum(sum(response(:,:,:))));
    
    %MVL: Confirms my belief that 30 blank zeros get added to the front and
    %here they are removed.
    totalME{n} = totalMotion(1+added_blankframe_no:length(totalMotion));
    
%     NEED TO CHECK FILTER PARAMS, THESE ME VS TIME PLOTS SHOULD MAKE SENSE
    if n == 100
        figure(5); 
        plot(t(1:length(totalME{n})),totalME{n}); 
        title(sprintf('dir = %d, coh=%1.2f, dur=%1.2f',data.direction(n),data.coherence(n),data.duration(n)));
        xlabel('time (s)'); 
        ylabel('motion energy (arb units)');
%       pause
    end
end

residuals = 0;
[totalME_norm_alignON, totalME_norm_alignRT] = ME_post_processing(data,totalME,residuals);


%% Was going to try using the Optic Flow Method (Bookmarked and Saved in Computer) to calculate motion energy 
