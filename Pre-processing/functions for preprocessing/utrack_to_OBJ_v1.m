
%% 
%
% This functions transfers data in tracksFinal to the OBJ structure. 
% Single objects in compound tracks are saves as independent objects. 
% The index to the orginal compound track
% is stored into OBJ.family. Events are stored into OBJ.events as follows:
% Column 1: frame at which compound track (family) starts
% Column 2: frame at which compound track (family)ends
% Column 3: frame at which single object starts
% Column 4: frame at which single object ends
% Column 5: parent object. (NaN in case of TRUE start. Absolute index of
% parent object in case of splitting).
% Column 6: child object. (NaN in case of TRUE end. Absolute index of child
% object in case of merging).


function utrack_to_OBJ_v1(fileName, time, path)
%% Input variables Variables
% filenames

movieName =[fileName '.tif'];
tracksFinal = [fileName '_tracking'];
outputFile = [fileName '_gui2_steps'];

%background variables
cellprofile= 0;                    % calculation of the mask of the cell (if the imaged area is bigger then the cell).
numBkgPx=200;                      % number of pixels to be summed up to calculate the regional background
maskradius=7;                      % parameter to define the mask around a center of intensity.
squaresize = 25;                   % area around the object that is used to find the proximal pixels in the background calculatio17
roiSize = 7;                       % region of interest for intensity calculation
bkgSize = 20;                      % region to draw a square around the object in the gui.

% filling the gap variables:
fillgp.wind=50;                    % max number of frames to be filled in between the trace
fillgp.frames=0;                  % number of frames added before/after in the trace intensity calculation in function of the option chosen.
fillgp.opt=0;                      % FillTheGap3 option (1=fill before and after, 2= fill to the end; =3 fill after of the number of frames indicated; =4 fill before of the indicated value).

% drawing
plotting = 0;                      % If setted to 1 are plotted the masks and the background point for each plane of the image.

% classification
minTrackLen=5;                     % minimum track lengh required to classify an object as class 1. Shorter lived objects are classified as class 2.

%% Load files:
display('loading files....');    tic ;
cd(path);                               %open directory
load(tracksFinal);                          %load tracksFinal file
% establishing frame numbers.
nframes = size(kalmanInfoLink,2);

% time
if ischar(time)
    % Loading time matrix
    timeMatrix =  dlmread(time,'\t',6,1);
    TME = timeMatrix(:,1)/1000;        %to have the time variable in sec
    for y=1:numel(TME)-1
        k(y) = TME(y+1)-TME(y);
    end
    Fr = mean(k);
    clear k y x timeMatrix tfname

elseif isnumeric(time) 
    TME = (0:nframes-1)*time;
    TME = TME';
    Fr = time;
end

%% Building the InFormatiOn array (IFO).
IFO.outfname = outputFile;
IFO.fpath = path; %movie path-- same as outputfile path
IFO.roiSize = roiSize;
IFO.bkgSize = bkgSize;       %essential for the gui

IFO.numFrames = nframes;
IFO.activeCh = [1 0 0];
IFO.numCh = 1;
IFO.ifnameRB = movieName;
[IFO.iRBmin, IFO.iRBmax] = minmaxpix(IFO.fpath, IFO.ifnameRB);

IFO.numBkgPx = numBkgPx;
IFO.maskradius = maskradius;    % parameter to define the mask around a center of intensity.
IFO.squaresize = squaresize;
clear outputFile pathname roiSize bkgSize movieName numBkgPx maskradius squaresize


%% Loading the values from tracksFinal to OBJ array variable (x,y,trj).
%To fill the variable is important to take into account the starting point
%(which frame), information stored in seqOfEvents first column.

% insertion of the timing part:
fprintf ('filling the variables..\n');tic;

nn=1;

OBJ.events=[];
lastTrackNum=0;

for ii =1:numel(tracksFinal);
    
    startFrame=tracksFinal(ii).seqOfEvents(1,1);
    
    for subobj=1:size(tracksFinal(ii).tracksCoordAmpCG,1)
        track2(nn).data = reshape(tracksFinal(ii).tracksCoordAmpCG(subobj,:),8,size(tracksFinal(ii).tracksCoordAmpCG,2)/8);
        
        OBJ.xR(nn,1:IFO.numFrames)=NaN;
        OBJ.yR(nn,1:IFO.numFrames)=NaN;
        OBJ.trjRB(nn,1:IFO.numFrames)=NaN;
        
        OBJ.xR(nn,startFrame:startFrame+size(track2(nn).data,2)-1) = track2(nn).data(1,:);
        OBJ.yR(nn,startFrame:startFrame+size(track2(nn).data,2)-1) = track2(nn).data(2,:);
        OBJ.trjRB(nn,startFrame:startFrame+size(track2(nn).data,2)-1) = track2(nn).data(4,:);
        
        nn=nn+1;
    end
    
    %getting family info
    %putting sequence of events data into OBJ.events
    
    
    localEvents=[];
    
    eventsNum=size(tracksFinal(ii).seqOfEvents,1);
    subObjNum=max(tracksFinal(ii).seqOfEvents(:,3));
    
    visStart=min(tracksFinal(ii).seqOfEvents(:,1));
    visEnd=max(tracksFinal(ii).seqOfEvents(:,1));
    
    localEvents(1:subObjNum,1:6)=NaN; %filling localEvents with NaN
    
    for event=1:eventsNum;
        
        eventFrame=(tracksFinal(ii).seqOfEvents(event,1));
        startOrEnd=(tracksFinal(ii).seqOfEvents(event,2));
        localTrackNum=(tracksFinal(ii).seqOfEvents(event,3));
        linkedTrack=(tracksFinal(ii).seqOfEvents(event,4))+lastTrackNum;
        
        if (isnan(linkedTrack)==1 && startOrEnd==1) %if a true start is present
            localEvents(localTrackNum,1)=visStart;
            localEvents(localTrackNum,3)=eventFrame;
            OBJ.family(localTrackNum+lastTrackNum,1)=ii; %putting family information into OBJ.family
        end
        
        if (isnan(linkedTrack)==1 && startOrEnd==2) %if a true end is present
            localEvents(localTrackNum,2)=visEnd;
            localEvents(localTrackNum,4)=eventFrame;
        end
        
        if (isnan(linkedTrack)==0 && startOrEnd==1) %if a merging is present
            localEvents(localTrackNum,1)=visStart;
            localEvents(localTrackNum,3)=eventFrame;
            localEvents(localTrackNum,5)=linkedTrack;
            OBJ.family(localTrackNum+lastTrackNum,1)=ii; %putting family information into OBJ.family
        end
        
        if (isnan(linkedTrack)==0 && startOrEnd==2) %if a splitting is present
            localEvents(localTrackNum,2)=visEnd;
            localEvents(localTrackNum,4)=eventFrame;
            localEvents(localTrackNum,6)=linkedTrack;
        end
        
    end
    
    OBJ.events = [OBJ.events; localEvents]; %appending local event data to OBJ.events
    lastTrackNum=lastTrackNum+size(localEvents,1);
    
end

IFO.numObj = numel(track2);
ntracks = length(track2);

clear ii nn track2 subobj


%% Preallocation of memory
%Preallocation of the memory to store the data.

OBJ.mxR = NaN(ntracks, nframes);
OBJ.myR = NaN(ntracks, nframes);
OBJ.btrjRB = zeros(ntracks, nframes);
OBJ.kn = zeros(ntracks,2);
OBJ.prp = zeros(ntracks,3);
OBJ.mstdRB = zeros(ntracks,nframes);
OBJ.bstdRB = zeros(ntracks,nframes);
OBJ.satR = zeros(ntracks);
OBJ.satL = zeros(ntracks);
OBJ.class = ones (1,size(OBJ.trjRB,1)); %introduced in version 1.1 to classify objects outside cell mask!
OBJ.classManual (1:numel(tracksFinal))=NaN;
clear ntracks nframes
toc

%% Part to fill the NaN values, fill the variables.

% 1) Filling NaN values in x,y.
% The guy has to address a position x,y for each object, each frame. It
% will addressed the first x,y position the object is detected to frames
% before it happearence and the last x,y position for the object after its
% disappearence.
for  u = 1:size (OBJ.xR,1);
    notnum  = find(isnan(OBJ.xR(u,:)));
    numbers = find(~isnan(OBJ.xR(u,:)));
    for j=1:numel(notnum)
        for m=1:numel(numbers)
            if abs(notnum(j)-numbers(m))  == min (abs(numbers-notnum(j)));
                OBJ.xR (u,notnum(j)) = OBJ.xR (u,numbers(m));
                OBJ.yR (u,notnum(j)) = OBJ.yR (u,numbers(m));
            end
        end
    end
end
clear u j m notnum numbers

% 2) Filling up the intensity values using FillTheGap function.
%[OBJ.trjRB] = FillTheGap3 (OBJ.trjRB,fillgp.wind,fillgp.frames,fillgp.opt);

% 3) round x,y
OBJ.mxR = round(OBJ.xR);
OBJ.myR = round(OBJ.yR);

% 4) Exclusion of objects at the edges
im = imread(IFO.ifnameRB,1);
[exclobj,OBJ] = edgeObjFind (OBJ, IFO, im, plotting);
OBJ.exclobj = exclobj;
clear exclobj im
toc
%% Calculation of the background and signal intensities

fprintf('background!\n');

for u=1:IFO.numFrames
  tic
 fprintf(['Starts ', num2str(u), ' frame\n'])
 % loading frame:
 im = imread(IFO.ifnameRB,u);
    % generation of the mask around center of intensities
[mask1] = maskgeneratorSingleplane2(im, OBJ.mxR(:,u), OBJ.myR(:,u),OBJ.trjRB(:,u), IFO.maskradius, plotting);
    %generation of the mask of the cell
if cellprofile==1
    [mask2] = mask_cell_profile3(im, plotting);
  else
       mask2 = ones(size(im,1),size(im,2));
  end
  bm_imRB = mask1.*mask2;         %multiplication of the 2 masks to get the final mask (area out of cell and areas of signal excluded from background calculations).

   toc
end

%% Giving classes to the objects
%assigning class=2 to fast objects

for  u = 1:size(OBJ.trjRB,1)
    startFr=OBJ.events(u,3);
    endFr=OBJ.events(u,4);
    
    if endFr-startFr<minTrackLen
        
        OBJ.class(u)=2;
    else
        OBJ.class(u)=1;       
    end
    
end

%setting class of objects outside mask = 0, objects inside mask = 1

if exist ([fileName '_msk.tif'],'file') == 2 %check if mask file exists
mask3=imread([fileName '_msk.tif']);
for u = 1:size (OBJ.trjRB,1)
    startFr=OBJ.events(u,3);
    startX=OBJ.mxR(u,startFr);
    startY=OBJ.myR(u,startFr);
    if mask3(startY,startX)==0
        OBJ.class(u)=0;
    end
end
end

fprintf ('Have a nice day!!!\n')

toc;
%%
save (IFO.outfname, 'OBJ', 'IFO', 'TME', 'im', 'bm_imRB', '-v7.3'); % to save files larger than 2 GB
close all;
