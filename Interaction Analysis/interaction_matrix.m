% Copyright (c) Lanoiselée and Calebiro 2023
% Please cite:
% Plasma membrane preassociation drives beta-arrestin coupling to receptors and activation
% Grimes, Koszegi, Lanoiselée et al.
% Cell 186, 1-18 (2023)
% https://doi.org/10.1016/j.cell.2023.04.018
% 
% This script is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY;
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.



function [int_matrix]=interaction_matrix(X,Y, frameStart,frameEnd,searchRadius)
% X and Y are row cells containing coordinated of the two channels for which colocalisations are computed 
% frameStart: first frame to visualize
% frameEnd: last frame to visualize

% searches for particles in two channels that fall closer than searchRadius (pixel).
% When two particles fall within searchRadius
% Generates a matrix 'interaction_matrix' with the following structure:
% column 1: frame at which interaction starts
% column 2: frame at which interaction ends
% column 3: channel number of first obj
% column 4: object number (in first channel)
% column 5: flag for true start (1, obj exist also at frame before
% interaction begins, 0 new appearance)
% column 6: flag for true end (1, obj exist also at frame after
% interaction ends, 0 object disappears)
% column 7: channel number of second obj
% column 8: object number (in second channel)
% column 9: flag for true start
% column 10: flag for true end

exclude_duplicate_int=0; %flag for excluding shorter duplicate interactions. It assumes that a particle in one channel at a given time can interac only with one particle in the second channel. If one particle interacts simulatenously with two particles (i.e. due to true + random colocalization), the shorter interaction is removed.
int_matrix=nan(max([size(X{1},1),size(X{2},1)]),10);

class=cell(1,numel(X));
for ch=1:numel(X)
    class{1,ch}=zeros(size(X{ch},1));
    duration=sum(~isnan(X{ch}),2);
    class{1,ch}(duration==0)=nan;
    class{1,ch}(duration<5)=2;
    class{1,ch}(duration>=5)=1;
end

%% interaction matrix
eventn=0;
for t=frameStart:frameEnd
    %% distance calculation
    Idx_nearest_neighbors = rangesearch([X{2}(:,t),Y{2}(:,t)],[X{1}(:,t),Y{1}(:,t)],searchRadius);
    for channel=1:size(X{1},1)
        nearpointIndexes=Idx_nearest_neighbors{channel}; %objs in ch2 that at frame t intract with obj(i) in ch 1
        %%
        if ~isempty(nearpointIndexes)
            for w=1:numel(nearpointIndexes)
                
                intobj_1=channel;
                intobj_2=nearpointIndexes(w);
                %% modify list of event to update interaction last time
                test_exist=int_matrix(:,4)==intobj_1 & int_matrix(:,8)==intobj_2 & int_matrix(:,2)==t-1;
                if any(test_exist)==1 % if at least one element of text_exist is nonzero
                    int_matrix(test_exist,2)=t;% interaction time is updated
                else
                    %% otherwise, i.e. new interaction:
                    if any([1 2 8]==class{1}(intobj_1)) && any([1 2 8]==class{2}(intobj_2)) %exclude objects falling outside mask (i.e. class=0)
                        eventn=eventn+1;
                        int_matrix(eventn,1) = t;
                        int_matrix(eventn,2) = t;
                        int_matrix(eventn,3) = 1; %temporary for ch1
                        int_matrix(eventn,4) = intobj_1;
                        int_matrix(eventn,5) = nan;%Ch(1).OBJ.events(intobj_1,5); %flag for true start, data taken direcly from Ch.OB.events (number of parent object or NaN for true start)
                        int_matrix(eventn,6) = nan;%Ch(1).OBJ.events(intobj_1,6); %flag for true end, data taken direcly from Ch.OB.events (number of child object or NaN for true end)
                        int_matrix(eventn,7) = 2; %temporary for ch2
                        int_matrix(eventn,8) = intobj_2;
                        int_matrix(eventn,9) = nan;%Ch(2).OBJ.events(intobj_2,5); %flag for true start, data taken direcly from Ch.OB.events (number of parent object or NaN for true start)
                        int_matrix(eventn,10) = nan;%Ch(2).OBJ.events(intobj_2,6); %flag for true end, data taken direcly from Ch.OB.events (number of child object or NaN for true end)
                        
                        class{1}(intobj_1)=8; %assigns class 8 to interacting objects
                        class{2}(intobj_2)=8;
                    end
                end
            end
        end
    end
end

%% exlude duplicate interactions
%in case an object interacts simultaneously with two objects, it removes
%the shorter interaction (i.e. presumably a random colocalization)
if exclude_duplicate_int==1
    %for objs in ch1
    for channel=1:eventn
        obj_ch1=int_matrix(channel,4);
        int_family=find(int_matrix(:,4)==obj_ch1); %index to all interactions that involve obj_ch1
        int_family=int_family';
        if numel(int_family)>1 %as long as there is more than one interaction involving the same obj_ch1
            for int1=int_family
                for int2=int_family(int_family~=int1)
                    %check if interaction is contained in any other interaction (i.e.
                    %it start after and ends before another one
                    if int_matrix(int1,1)>=int_matrix(int2,1) && int_matrix(int1,2)<=int_matrix(int2,2)
                        %then remove the interaction
                        int_matrix(int1,1:8)=NaN;
                    end
                end
            end
        end
    end
    %now for objs in ch2
    for channel=1:eventn
        obj_ch2=int_matrix(channel,8);
        int_family=find(int_matrix(:,4)==obj_ch2); %index to all interactions that involve obj_ch1
        int_family=int_family';
        if numel(int_family)>1 %as long as there is more than one interaction involving the same obj_ch1
            for int1=int_family
                for int2=int_family(int_family~=int1)
                    %check if interaction is contained in any other interaction (i.e.
                    %it start after and ends before another one
                    if int_matrix(int1,1)>=int_matrix(int2,1) && int_matrix(int1,2)<=int_matrix(int2,2)
                        %then remove the interaction
                        int_matrix(int1,1:8)=NaN;
                    end
                end
            end
        end
    end
end
% remove empty entries
int_matrix(all(isnan(int_matrix),2),:)=[];
end
