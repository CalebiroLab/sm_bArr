% OBJ_reshape:
% this new routine reshapes obj trajectories and intensities to eliminate
% the gaps due to:
% 1) obj merging and splitting
% 2) obj temporary dispappearance (due to blinking or mssed detection)

function OBJ_reshape (filename, global_folders)

load ([global_folders.localfolder, filesep, filename, '_gui2_steps.mat']);


OBJ2.events=OBJ.events; %create a copy of OBJ
OBJ2.xR=OBJ.xR; %create a copy of OBJ
OBJ2.yR=OBJ.yR; %create a copy of OBJ
OBJ2.mxR=OBJ.mxR; %create a copy of OBJ
OBJ2.myR=OBJ.myR; %create a copy of OBJ
OBJ2.trjRB=OBJ.trjRB; %create a copy of OBJ
OBJ2.btrjRB=OBJ.btrjRB;
OBJ2.family=OBJ.family; %create a copy of OBJ
OBJ2.class=OBJ.class; %create a copy of OBJ
OBJ2.classManual=OBJ.classManual; %create a copy of OBJ

nan_intensities=find(isnan(OBJ2.trjRB)); %set coordinates to NaN if intensity is NaN!!!!
OBJ2.xR(nan_intensities)=NaN;
OBJ2.yR(nan_intensities)=NaN;

objnum=size(OBJ2.xR,1);
framenum=size(OBJ2.xR,2);

%going through all objects and fuse spltting object with parent object,
%this will be reiterated

%first identify the first potential merging point
last_obj=objnum;

while any(~isnan(OBJ2.events(:,5))) || any(~isnan(OBJ2.events(:,6))) %as long as there are splitting or merging objects
ob=1;

while ob<=last_obj
 
if ~isnan(OBJ2.events(ob,1)) 
    
merge_fstart=[]; 
split_fstart=[];
merge_obj=[];
split_obj=[];

%find all potential merging objs and frames
ob_fstart=OBJ2.events(ob,3);%the starting frame of the object ob
merging_objs_temp=find(OBJ2.events(:,6)==ob);
merging_objs=[-1; merging_objs_temp]; %all objs that merge onto ob, plus the hypothetical -1 obj to consider the possibility that the obj was already merged from the beginiing  
merging_objs_fr=[ob_fstart-1; (OBJ2.events(merging_objs_temp,4))]; %the frames at which the merging onjects end, plus the begining of ob for the -1 case
merging_objs_fr=merging_objs_fr';

%find all potential splitting objs and frames
ob_fend=OBJ2.events(ob,4);
splitting_objs_temp=find(OBJ2.events(:,5)==ob);
splitting_objs=[splitting_objs_temp; -1]; %all objs that spit from ob after the merge_fstart, plus -1, hypothetical obj that splits after the end of the ob
%if ~isempty(splitting_objs)

splitting_objs_fr=[(OBJ2.events(splitting_objs_temp,3)); ob_fend+1] ; %the frames at which the  splitting objs start

%now identify an appropriate pair of merge and split
%start with the earliest splitting frame and relative obj
splitting_objs_fr=splitting_objs_fr';
[split_fstart split_obj_index]=min(splitting_objs_fr);
split_obj=splitting_objs(split_obj_index); %the identified splitting object

%now take the closer merging frame and obj
merge_fend=max(merging_objs_fr(merging_objs_fr<split_fstart)); %the selected merging frame, i.e. the one closest to the split frame

if ~isempty(merge_fend)
merging_obj_index=find(merging_objs_fr==merge_fend,1); %the index to the corresponding obj
merge_obj=merging_objs(merging_obj_index);
end

if merge_obj ~=-1 %in case there was a previous merge, identify the merge_fstart (otherwise it�s NaN)
merge_fstart=OBJ2.events(merge_obj,3);
end

if isempty(merge_obj)
    
    merge_obj=-1;
end

%now building the fused object
if split_obj>-1 || merge_obj>-1 %there is at least a detected split or merge  
    if split_obj==-1 %if there is no splitting object, create a new empty one
        split_obj=last_obj+1;
        last_obj=last_obj+1;
        OBJ2.xR(split_obj,1:framenum)=NaN; %create a new empty obj
        OBJ2.yR(split_obj,1:framenum)=NaN; %idem for yR
        OBJ2.mxR(split_obj,1:framenum)=NaN; %create a new empty obj
        OBJ2.myR(split_obj,1:framenum)=NaN; %idem for yR
        OBJ2.events(split_obj,1:6)=OBJ2.events(ob,1:6);
        OBJ2.trjRB(split_obj,1:framenum)=NaN;
        OBJ2.btrjRB(split_obj,1:framenum)=NaN;
        OBJ2.class(split_obj)=OBJ2.class(ob);
        %OBJ2.classManual(split_obj)=NaN;
        OBJ2.family(split_obj,1)=OBJ2.family(ob,1);
    end  

if merge_obj>-1 %i.e. a previous merging object was identified
OBJ2.xR(split_obj,merge_fstart:merge_fend)=OBJ2.xR(merge_obj,merge_fstart:merge_fend); %move from merging object to target obj
OBJ2.yR(split_obj,merge_fstart:merge_fend)=OBJ2.yR(merge_obj,merge_fstart:merge_fend); %move from merging object to target obj
OBJ2.xR(merge_obj,merge_fstart:merge_fend)=NaN; %delete the orginal merging object
OBJ2.yR(merge_obj,merge_fstart:merge_fend)=NaN; %delete the original merging object
OBJ2.family(merge_obj)=NaN; %delete the original merging object
OBJ2.class(merge_obj)=NaN; %delete the original merging object

OBJ2.mxR(split_obj,merge_fstart:merge_fend)=OBJ2.mxR(merge_obj,merge_fstart:merge_fend); %move from merging object to target obj
OBJ2.myR(split_obj,merge_fstart:merge_fend)=OBJ2.myR(merge_obj,merge_fstart:merge_fend); %move from merging object to target obj
OBJ2.mxR(merge_obj,merge_fstart:merge_fend)=NaN; %delete the orginal merging object
OBJ2.myR(merge_obj,merge_fstart:merge_fend)=NaN; %delete the original merging object

OBJ2.trjRB(split_obj,merge_fstart:merge_fend)=OBJ2.trjRB(merge_obj,merge_fstart:merge_fend); %move intensity data from merging object to target obj
OBJ2.btrjRB(split_obj,merge_fstart:merge_fend)=OBJ2.btrjRB(merge_obj,merge_fstart:merge_fend); 
OBJ2.trjRB(merge_obj,merge_fstart:merge_fend)=NaN; %delete the original intensity data of merging object
OBJ2.btrjRB(merge_obj,merge_fstart:merge_fend)=NaN;
OBJ2.events(split_obj,3)=OBJ2.events(merge_obj,3); %....updates the start frame
OBJ2.events(split_obj,5)=OBJ2.events(merge_obj,5);
OBJ2.events(merge_obj,1:6)=NaN; % update the OBJ.events

%update any other object that mrged/split with the deleted merge obj
old_merge_obj=find(OBJ2.events(:,6)==merge_obj);
OBJ2.events(old_merge_obj,6)=split_obj;
old_split_obj=find(OBJ2.events(:,5)==merge_obj);
OBJ2.events(old_split_obj,5)=split_obj;
end

OBJ2.xR(split_obj,merge_fend+1:split_fstart-1)=OBJ2.xR(ob,merge_fend+1:split_fstart-1);
OBJ2.yR(split_obj,merge_fend+1:split_fstart-1)=OBJ2.yR(ob,merge_fend+1:split_fstart-1);
OBJ2.mxR(split_obj,merge_fend+1:split_fstart-1)=OBJ2.mxR(ob,merge_fend+1:split_fstart-1);
OBJ2.myR(split_obj,merge_fend+1:split_fstart-1)=OBJ2.myR(ob,merge_fend+1:split_fstart-1);
OBJ2.trjRB(split_obj,merge_fend+1:split_fstart-1)=OBJ2.trjRB(ob,merge_fend+1:split_fstart-1);
OBJ2.btrjRB(split_obj,merge_fend+1:split_fstart-1)=OBJ2.btrjRB(ob,merge_fend+1:split_fstart-1);


if merge_obj==-1
OBJ2.events(split_obj,3)=OBJ2.events(ob,3); %update starting frame in case there was no merge;
OBJ2.events(split_obj,5)=NaN; %update starting frame in case there was no merge;
end
end
end


ob=ob+1;

end

end

%%%%%%%%%%%
%now filling the gaps in x and y coordinates
objnum=size(OBJ2.xR,1);

for ob=1:objnum  
    obj_fstart=find(~isnan(OBJ2.xR(ob,:)),1); %do not connsider NaNs before the beginning of an object
    f=obj_fstart;
    
    while f<framenum 
        if isnan(OBJ2.xR(ob,f))% a gap start is found
        gap_fstart=f;
        f2=f+1;
        gap_end_fnd=0;
            while gap_end_fnd==0 && f2<=framenum
                if ~isnan(OBJ2.xR(ob,f2)) % a gap end is found
                    gap_fend=f2-1;
                    OBJ2.xR(ob,gap_fstart-1:gap_fend+1)=linspace(OBJ2.xR(ob,gap_fstart-1),OBJ2.xR(ob,gap_fend+1),gap_fend-gap_fstart+3);
                    OBJ2.yR(ob,gap_fstart-1:gap_fend+1)=linspace(OBJ2.yR(ob,gap_fstart-1),OBJ2.yR(ob,gap_fend+1),gap_fend-gap_fstart+3);
                    gap_end_fnd=1;
                end
                f2=f2+1;
            end
        end
        f=f+1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%


clear cp f f2 framenum gap_end_fnd gap_fend gap_fstart last_obj merge_fend merge_fstart merg_obj merging_obj_index merging_objs merging_objs_fr merging_objs_temp nan_intensities
clear ob ob_fend ob_fstart obj_duration1 obj_duration2 obj_fstart objnum old_merge_obj old_split_obj split_fstart split_obj split_obj_index splitting_objs splitting_objs_fr
clear splitting_objs_temp
 IFO.numObj2=size(OBJ2.xR,1);
save ([global_folders.localfolder, filesep, filename, '_gui2_steps'], 'OBJ','OBJ2','IFO','TME', '-v7.3');

end