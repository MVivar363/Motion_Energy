% Miguel VIvar-Lazo
% Fetsch Lab
% 02/21/2020
% Getting Data Ready For Condor

%% For RDM Training their might be sets where we present only one target and not both, eliminate one targets
choice(twoTarget == 0)          = [];
coherence(twoTarget == 0)       = [];
correct(twoTarget == 0)         = [];
dates(twoTarget == 0)           = [];
direction(twoTarget == 0)       = [];
dotPos(twoTarget == 0)          = [];
rt(twoTarget == 0)              = [];
taskRT(twoTarget == 0)          = [];
twoTarget(twoTarget == 0)       = [];

%% Assuming you have everthing ready to be loaded into a best_data structure
best_data.dotPos = dotPos;
best_data.ppd = 39.0951; %Most of these numbers are standard in RigA
best_data.screen_rect = [0 0 1920 1080];
best_data.refresh = 120;
best_data.aperture = [0 0 6 6]; %[x,y,dim]
best_data.speed = 6; %deg/s , although in the constructing of dots we say speed is 6, in ME file it should be 3 because of how it changes spatial parameters
best_data.dotSize = 4; %Radius or Diameter?
best_data.prefDir = 0;
best_data.coherence = coherence;
best_data.direction = direction;
best_data.duration = RT; 

% After this data should be ready to put in the Condor machince