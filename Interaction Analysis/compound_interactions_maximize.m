% Copyright (c) Lanoiselée and Calebiro 2023
% Please cite:
% Plasma membrane preassociation drives beta-arrestin coupling to receptors and activation
% Grimes, Koszegi, Lanoiselée et al.
% Cell 186, 1-18 (2023)
% https://doi.org/10.1016/j.cell.2023.04.018
% 
% This script is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY;
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
function [interaction_data] = compound_interactions_maximize(int_matrix,gap_wind)

%description of int_families strcuture array that contain most relevant
%information
%int_families.
%       .ints (the indexes of interactions in a family)
%       .int_linking_matrix (the linking matrix, where columns are the reference interactions, rows are the possible forward linking interractions, ones indicate a possible interaction)
%
%       .int_linking_matrix_mask_init (the single biunivocal initialized interaction matrix. The starting point for the optimization)
%       .int_linking_matrix_mask_current (the cuurent and best so far)
%       .int_linking_matrix_mask_test (the new test matrix after the perturb
%
%       .linked_ints_init
%       .nonlinked_ints_init
%       .linked_int_durations_init - .linked_int_durations_init
%       .nonlinked_int_durations_init - .nonlinked_int_durations_init
%       .linked_int_durationsAvg_init - .linked_int_durationsAvg_init
%
%       .linked_ints_current
%       .nonlinked_ints_current
%       .linked_int_durations_current - .linked_int_durations_current
%       ..nonlinked_int_durations_current - .nonlinked_int_durations_current

%       .linked_int_durationsAvg_cuurent - .linked_int_durationsAvg_current
%
%       .linked_ints_test
%       .nonlinked_ints_test
%       .linked_int_durations_test - .linked_int_durations_test
%       .nonlinked_int_durations_test - .nonlinked_int_durations_test
%       .linked_int_durationsAvg_test - .linked_int_durationsAvg_test
%
%
%       .track(iter).linked_int_durations
%       .track(iter).nonlinked_int_durations
%       .linked_int_durationsAvg_track

iter_max=2000;


%% %%%%%%%%%%%%%%%%%%%%%%
%STEP1 - generate the int_linking_matrix starting from the int_matrix (the
%previous interaction matrix). Each column correspond to an interaction and
%each row to an interaction that could link with it. The diagonal must
%contain zeros

int_num=size(int_matrix,1); %the number of interactions

int_families=[]; %structure array int_families.int, contains all the interactions belonging to the same family of interaction, i.e. all the interactions that could be linked with each other, later used for maximization!!!!!!!


%% filling the int_linking_matrix with ones for all possible linking interactions

[int_linking_matrix,int_linking_matrix_durations,int_families] = interaction_linking(int_matrix,gap_wind);

%% sort and keep only unique int_families.ints
int_families_num=numel(int_families);

for i=1:int_families_num
    int_families(i).ints=sort(unique(int_families(i).ints));
end

in_families_save=int_families;


%% set the diagonals to zeros!!!! (one interaction cannot link to itself!!!)
% diag_fill=zeros(size(int_matrix,1),1);
% int_linking_matrix=spdiags(diag_fill,0,int_linking_matrix);

%int_linking_matrix=int_linking_matrix-diag(int_linking_matrix);
% int_linking_matrix_durations=spdiags(diag_fill,0,int_linking_matrix_durations);

% int_linking_matrix_durations=int_linking_matrix_durations-diag(int_linking_matrix_durations);

%% %%%%%%%%%%%%%
%generate a initialization duration matrix in which one level linking
%produces maximal linked durations
int_linking_matrix_durations_temp=int_linking_matrix_durations;
int_linking_matrix_durations_init=sparse(int_num,int_num);%(1:int_num,1:int_num)=0;
int_linking_matrix_mask_init=sparse(int_num,int_num);%(1:int_num,1:int_num)=0;

% int_linking_matrix_durations_init(1:int_num,1:int_num)=0;
% int_linking_matrix_mask_init(1:int_num,1:int_num)=0;


[row,col]=find(int_linking_matrix_durations_temp>0);
list_temp=int_linking_matrix_durations_temp(sub2ind(size(int_linking_matrix_durations_temp),row,col));
list_dur=[row,col,list_temp];
list_dur=sortrows(list_dur,3,'descend');

number_int=size(list_dur,1);
while number_int>0
    int_linking_matrix_mask_init(list_dur(1,1),list_dur(1,2))=1; %building the init mask matrix!
    int_linking_matrix_durations_init(list_dur(1,1),list_dur(1,2))=list_dur(1,3); %building the init duration
    
    A=list_dur(:,1)==list_dur(1,1) | list_dur(:,2)==list_dur(1,2);
    number_int=number_int-sum(A);
    list_dur(A,:)=[];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %STEP2 %pick up only a family of interactions and for eack of them find the best


for ifam=1:int_families_num
    iter=1;
    %generate an int_linking_matrix just for the ints in the family
    ints_in_family=int_families(ifam).ints;
    ints_in_family_num=numel(ints_in_family);
    int_families(ifam).int_linking_matrix=int_linking_matrix(ints_in_family,ints_in_family);
    
    int_families(ifam).int_linking_matrix_mask_init=int_linking_matrix_mask_init(ints_in_family,ints_in_family); %the single biunivocal linking matrix for which the durations of
    %the linked ints are maximal, this is the starting matrix
    %calculate the all linked durations for the starting matrix
    [int_families(ifam).linked_ints_init int_families(ifam).linked_int_durations_init int_families(ifam).nonlinked_ints_init int_families(ifam).nonlinked_int_durations_init]=linked_interactions_durations(int_families(ifam).int_linking_matrix_mask_init,ints_in_family,int_matrix);
    int_families(ifam).track(iter).linked_int_durations =int_families(ifam).linked_int_durations_init; %to keep track of the evolution of the optimizaation
    int_families(ifam).track(iter).nonlinked_int_durations=int_families(ifam).nonlinked_int_durations_init;
    
    %calculate the average lenght of the linked interactions %Willl be
    %used for maximization
    int_families(ifam).linked_int_durationsAvg_init=sum(int_families(ifam).linked_int_durations_init)/numel(int_families(ifam).linked_int_durations_init);
    int_families(ifam).compound_int_num_init=numel(int_families(ifam).nonlinked_ints_init)+numel(int_families(ifam).linked_ints_init);
    int_families(ifam).linked_int_durationsAvg_track=int_families(ifam).linked_int_durationsAvg_init;
    int_families(ifam).nonlinked_int_num_track=numel(int_families(ifam).nonlinked_ints_init);
    int_families(ifam).compound_int_num_track=int_families(ifam).compound_int_num_init;
    %%%%%%%%
    int_families(ifam).int_linking_matrix_mask_current=int_families(ifam).int_linking_matrix_mask_init; %the currently selected biunivocal interaction matrix
    int_families(ifam).linked_ints_current=int_families(ifam).linked_ints_init;
    int_families(ifam).linked_int_durations_current=int_families(ifam).linked_int_durations_init;
    int_families(ifam).nonlinked_ints_current=int_families(ifam).nonlinked_ints_init;
    int_families(ifam).nonlinked_int_durations_current=int_families(ifam).nonlinked_int_durations_init;
    int_families(ifam).linked_int_durationsAvg_current=int_families(ifam).linked_int_durationsAvg_init;
    int_families(ifam).compound_int_num_current=int_families(ifam).compound_int_num_init;
    if numel(int_families(ifam).ints)>1 %i.e. there are atl least two interactions in a family
        stop=0;
        while iter<iter_max && stop==0
            %randomly perturb the starting matrix (->test matrix)
            [int_families(ifam).int_linking_matrix_mask_test] = matrixPerm(int_families(ifam).int_linking_matrix_mask_current,int_families(ifam).int_linking_matrix,1);
            
            %recalculate the all linked durations for the perturbed (test) matrix
            
            [int_families(ifam).linked_ints_test,int_families(ifam).linked_int_durations_test,int_families(ifam).nonlinked_ints_test,int_families(ifam).nonlinked_int_durations_test]=linked_interactions_durations(int_families(ifam).int_linking_matrix_mask_test,ints_in_family,int_matrix);
            
            int_families(ifam).linked_int_durationsAvg_test=sum(int_families(ifam).linked_int_durations_test)/numel(int_families(ifam).linked_int_durations_test);
            
            int_families(ifam).compound_int_num_test=numel(int_families(ifam).nonlinked_ints_test)+numel(int_families(ifam).linked_ints_test);
            
            if int_families(ifam).compound_int_num_test<=int_families(ifam).compound_int_num_current %chaned in test!!!!!!!!!
                
                int_families(ifam).int_linking_matrix_mask_current=int_families(ifam).int_linking_matrix_mask_test; %the currently selected biunivocal interacionn matrix!
                int_families(ifam).linked_ints_current=int_families(ifam).linked_ints_test;
                int_families(ifam).linked_int_durations_current=int_families(ifam).linked_int_durations_test;
                int_families(ifam).nonlinked_ints_current=int_families(ifam).nonlinked_ints_test;
                int_families(ifam).nonlinked_int_durations_current=int_families(ifam).nonlinked_int_durations_test;
                int_families(ifam).linked_int_durationsAvg_current=int_families(ifam).linked_int_durationsAvg_test;
                int_families(ifam).compound_int_num_current=int_families(ifam).compound_int_num_test;
            end
            int_families(ifam).track(iter).linked_int_durations =int_families(ifam).linked_int_durations_test; %to keep track of the evolution of the optimizaation
            int_families(ifam).track(iter).nonlinked_int_durations=int_families(ifam).nonlinked_int_durations_test;
            
            int_families(ifam).linked_int_durationsAvg_track=[int_families(ifam).linked_int_durationsAvg_track int_families(ifam).linked_int_durationsAvg_test];
            
            int_families(ifam).nonlinked_int_num_track=[int_families(ifam).nonlinked_int_num_track numel(int_families(ifam).nonlinked_ints_test)];
            int_families(ifam).compound_int_num_track=[int_families(ifam).compound_int_num_track int_families(ifam).compound_int_num_test];
            
            iter=iter+1;
            %             disp(num2str(int_families(ifam).compound_int_num_current))
            if numel(int_families(ifam).compound_int_num_current)==1
                stop=1;
            end
        end
    end
    
end
%toc
%%%%%%%%%%%%%
%retrieving the optimized data and store them in compound_ints and
%compound_ints_durations - the real output of the routine
w=1;
compound_ints_all=[];

for ifam=1:int_families_num
    for l=1:numel(int_families(ifam).linked_ints_current)
        if ~isempty(int_families(ifam).linked_ints_current(l).ints)
            compound_ints(w).ints=int_families(ifam).linked_ints_current(l).ints;
            compound_ints_all=[compound_ints_all compound_ints(w).ints];
            w=w+1;
        end
    end
end

for ifam=1:int_families_num
    for l=1:numel(int_families(ifam).nonlinked_ints_current)
        if ~isempty(int_families(ifam).nonlinked_ints_current(l))
            compound_ints(w).ints=int_families(ifam).nonlinked_ints_current(l);
            compound_ints_all=[compound_ints_all compound_ints(w).ints];
            w=w+1;
        end
    end
end


%%%%%%%%%%%%%%

list_compound_start=Inf*ones(numel(compound_ints),1);
list_compound_end=ones(numel(compound_ints),1);

%filling the OBJ3 matrix for each channel

% list_compound_int=[];

list_compound_int=compound_ints;
for n=1:numel(compound_ints)
    for i=1:numel(compound_ints(n).ints)
        current_int=compound_ints(n).ints(i);
        current_int_start=int_matrix(current_int,1);
        current_int_end=int_matrix(current_int,2);
        current_obj=int_matrix(current_int,4)+int_matrix(current_int,8);
        if ~isnan(current_obj)
            if current_int_start<= list_compound_start(n)
                list_compound_start(n)=current_int_start;
            end
            if current_int_end>=list_compound_end(n)
                list_compound_end(n)=current_int_end;
            end
        end
    end
    if list_compound_start(n)> list_compound_end(n)
        list_compound_start(n)=NaN;
        list_compound_end(n)=NaN;
    end
end



comp_ints_starts=list_compound_start;
comp_ints_ends=list_compound_end;

compound_ints_durations=comp_ints_ends-comp_ints_starts+1;
interaction_data=struct;
interaction_data.list_compound_start= list_compound_start;
interaction_data.list_compound_end= list_compound_end;
interaction_data.list_compound_int= list_compound_int;
interaction_data.compound_ints_durations=compound_ints_durations;

end
