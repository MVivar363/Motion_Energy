function totalME = MotionEnergyCondor(condorJobId) %,saveDotPositions)

Job = condorJobId + 1; %Have to add +1 because the computer has a based of 0

%data will saved in a MAT file and loaded appropriately
load('Condor_MultiDirectionData_0824_0911.mat', 'best_data'); %load file that has all the information and variables for this set up
data = best_data;

AP_D    = data.aperture(3); % diam of aperture in degrees
d_ppd 	= AP_D * data.ppd;	% diam of aperture in pixels
left = data.screen_rect(3)/2 - AP_D/2*data.ppd;
top = data.screen_rect(4)/2 - AP_D/2*data.ppd;
right = data.screen_rect(3)/2 + AP_D/2*data.ppd;
bottom = data.screen_rect(4)/2 + AP_D/2*data.ppd;
AP_RECT = [left top right bottom]; % aperture in 'rect' format
ap_off  = [AP_RECT(1)+d_ppd/2 AP_RECT(2)+d_ppd/2];  % aperture center in screen coords
totalME = cell(size(data.coherence))';
totalLME = cell(size(data.coherence))';
totalRME = cell(size(data.coherence))';

%mvl: select parameters for the Job
%K = 30 to 120
%data.speed = 3 to 12
% k_ = 30:10:120;
% speed = 3:1:12;
% [kTable, speedTable] = meshgrid(k_, speed); 
direction = 0:15:359;

%k = kTable(Job); %Pick paramteres for this particular Job
%data.speed = speedTable(Job); %Pick paramteres for this particular Job
data.prefDir = direction(Job); %pick parameter from list

tic %start the time
% can try to use parfor eventually, to speed this up
for n = 1:length(data.coherence) %length(PDS.data) % trials
    dotxy = data.dotPos{n};
    % old ME code assumes dotXY is organized as [nframes, ndots, xy], but
    % PLDAPS generateDots saves it as [ndots, xy, nframes], so:
    dotxy = permute(dotxy,[3 1 2]);
    
    nframes = size(dotxy,1);   
    
    % make a 'pixel matrix' of 1s and 0s covering the aperture for each frame
    pixelMat = zeros(ceil(d_ppd),ceil(d_ppd),nframes);
    for f = 1:nframes
        xPix = dotxy(f,:,1)'+ceil(d_ppd/2); xPix(isnan(xPix))=[];
        yPix = dotxy(f,:,2)'+ceil(d_ppd/2); yPix(isnan(yPix))=[];
        
        switch data.dotSize % making some guesses about how pixels for each dot extend from 'center' coordinate
            % math problem for Miguel: simplify this :)
            case 2 % 4 pixels per dot
                xPix = [xPix ; xPix   ; xPix+1 ; xPix+1];
                yPix = [yPix ; yPix+1 ; yPix   ; yPix+1];
            case 3 % 9 pixels per dot
                xPix = [xPix-1 ; xPix-1 ; xPix-1 ; xPix   ; xPix ; xPix   ; xPix+1 ; xPix+1 ; xPix+1];
                yPix = [yPix-1 ; yPix   ; yPix+1 ; yPix-1 ; yPix ; yPix+1 ; yPix-1 ; yPix   ; yPix+1];
            case 4 % 16 pixels per dot
                xPix = [xPix-1 ; xPix-1 ; xPix-1 ; xPix-1 ; xPix   ; xPix ; xPix   ; xPix   ; xPix+1 ; xPix+1 ; xPix+1 ; xPix+1 ; xPix+2 ; xPix+2 ; xPix+2 ; xPix+2];
                yPix = [yPix-1 ; yPix   ; yPix+1 ; yPix+2 ; yPix-1 ; yPix ; yPix+1 ; yPix+2 ; yPix-1 ; yPix   ; yPix+1 ; yPix+2 ; yPix-1 ; yPix   ; yPix+1 ; yPix+2];
        end
                
        % remove pixels outside the ap (again, guessing that PTB does this)
        deleteMe = xPix>ceil(d_ppd) | yPix>ceil(d_ppd) | xPix<1 | yPix<1;
        xPix(deleteMe) = []; yPix(deleteMe) = [];
        for p = 1:length(xPix)
            pixelMat(xPix(p),yPix(p),f) = 1;
        end
    end
    
    % mvl: What we are going to have to do is construct a table that will %
    % select the Paramteres we want given the 'condorJobId' number %
    
    % construct the filters for motion energy
    sigmult = data.speed*0.125+0.204; % from ad hoc linear regression 
    % first adjust filters for speed

    % SPATIAL FILTERS f1,f2 (see Adelson & Bergen 1985, Kiani et al. 2008, Bollimunta et al. 2012)
    sigma_c = 0.35 * sigmult; % controls period ...? ; mvl: or controls the width of the curve
    sigma_g = 0.05 * sigmult; % controls width of Gaussian envelope along orthogonal dir, MVL: seems to contrll just a bit of width but more the description of both end of the tails (left and right, making it bigger in amplitude but also longer in tail or wider)
    
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
    %x = linspace(-AP_D/20,AP_D/20,ceil(d_ppd)); % need same size as aperture for Fourier method
    x = linspace(-AP_D/2,AP_D/2,ceil(d_ppd)); %mvl: Going to use Degrees and not pixels
    
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
    k = 60; % 60 is standard, but instead set above, depending on speed
    
    %MVL: This is the factorial equation we saw in the Adelson paper
    %(parameters from paper as well)
    n1=3; n2=5;
    g1 = @(t) (k*t).^n1 .* exp(-k*t) .* (1/factorial(n1) - (k*t).^2/factorial(n1+2)); % fast
    g2 = @(t) (k*t).^n2 .* exp(-k*t) .* (1/factorial(n2) - (k*t).^2/factorial(n2+2)); % slow
    added_blankframe_no = 30;   % pad with zeros, for least the length of the temporal filter
    %mvl: K makes the period or width of curves smaller. Frequency higher.
    
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
    %A(abs(A)<1e-6) = 0; B(abs(B)<1e-6) = 0; % threshold to avoid floating point weirdness
    
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
    %the Adelson papaer
    
    % note: switching to terminology from Eero Simoncelli's motionTutorial.m
    even_fast = f1_t.*g1_xy; %mvl: shouldnt we be doing a convolution here?
    odd_slow = f2_t.*g2_xy;  %mvl: discrete convolution is equivalent to taking a dot product between the filter weights and the values underneath the filter. So i THINK this works as a convolution
    odd_fast = f2_t.*g1_xy;  %mvl: I'm pretty sure this is not a dot product. Unfortunately, so should we be doing a convolution here?
    even_slow = f1_t.*g2_xy; 
    %mvl: by extending the 2D filter in time and Time Filter in some weird way and then multiplying each 2D Space index with its equivalent 1D time vector, 
    %we are performing the correct dot product (its confuding but look at page 6 in Adleson and Begen paper. they do a Nx1 * 1xN product giving you a 2D matrix. 
    %What we want here is to grab a 2D space vector and a 1D time vector, then do a proper dot
    %product that turns it into a 3D filter. Which is what they do here in
    %a smart and intereting way.
     
    
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
    fright1 = fftn(right1); fright2 = fftn(right2);
    fleft1  = fftn(left1);  fleft2  = fftn(left2);
    fstim = fftn(pixelMatN);
    
    %MVL: This all makes sense. 
    right1Resp = abs(ifftn( fstim .* fright1 ));
    right2Resp = abs(ifftn( fstim .* fright2 ));
    left1Resp  = abs(ifftn( fstim .* fleft1  ));
    left2Resp  = abs(ifftn( fstim .* fleft2  ));
    response = (right1Resp.^2 + right2Resp.^2) - (left1Resp.^2 + left2Resp.^2);
    leftResponse = (left1Resp.^2 + left2Resp.^2);
    rightResponse = (right1Resp.^2 + right2Resp.^2);

    if 0 %plotflag
        figure(21);
        cmax = max([max(max(max(response))) abs(min(min(min(response))))]);
        for q = 1:size(response,3)
            clf; contourf(response(:,:,q)); colorbar;
            title(num2str(q)); caxis([-cmax cmax]);
            pause(0.05)
        end    
    end
    
    totalMotion = squeeze(sum(sum(response(:,:,:))));
    leftMotion = squeeze(sum(sum(leftResponse(:,:,:))));
    rightMotion = squeeze(sum(sum(rightResponse(:,:,:))));
    keyboard
    %MVL: Confirms my belief that 30 blank zeros get added to the front and
    %here they are removed.
    totalME{n} = totalMotion(1+added_blankframe_no:length(totalMotion));
    totalLME{n} = leftMotion(1+added_blankframe_no:length(leftMotion));
    totalRME{n} = rightMotion(1+added_blankframe_no:length(rightMotion));
end

toc %stop the time
%system('hostname') %What does this do?
save(['motionEnergyTheta' num2str(Job)], 'totalME', 'totalLME', 'totalRME') %save data





