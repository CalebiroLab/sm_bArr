% Copyright (c) Lanoiselée and Calebiro 2023
% Please cite:
% Plasma membrane preassociation drives beta-arrestin coupling to receptors and activation
% Grimes, Koszegi, Lanoiselée et al.
% Cell 186, 1-18 (2023)
% https://doi.org/10.1016/j.cell.2023.04.018
% 
% This script is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY;
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

function [linked_ints, linked_int_durations, nonlinked_ints, nonlinked_int_durations] = linked_interactions_durations(int_linking_matrix_masked,ints_in_family,int_matrix)

%Gives a list containing the duration of all the linked interactions as
%well as a separate list for all the non linked interactions
%   Detailed explanation goes here

%%%%%this part calculates  the duration of all interactions
%%%%%with a given permutation, after limking them. IMPORTANT:

linked_ints.ints=[];
nonlinked_ints=[];
linked_int_durations=0;
link_number=0;
q=0;
r=1;

ints_num=numel(ints_in_family);
%now for each possible starting interaction
poss_start_ints=find(sum(int_linking_matrix_masked,2)==0); %find all possible starting interaction (i.e. all those that
%have no backward link)

if ~isempty(poss_start_ints) %in case there are no not linked interactions (circular linking, skip calculations)
for ci=1:numel(poss_start_ints)
current_int_index=poss_start_ints(ci);
current_int=ints_in_family(current_int_index); %the real interaction number
current_int_duration=int_matrix(current_int,2)-int_matrix(current_int,1)+1; %the duration of the current int
linked_int_durations_temp=current_int_duration;
int_start=1;

while current_int_index<=ints_num
joining_int_index=find(int_linking_matrix_masked(:,current_int_index)==1);
if ~isempty(joining_int_index)

joining_int=ints_in_family(joining_int_index); %the real interaction number
joining_int_duration=int_matrix(joining_int,2)-int_matrix(joining_int,1)+1; %the duration of the current int

if int_start==1%the start of an interaction
q=q+1;
linked_int_durations(q)=linked_int_durations_temp;
linked_ints(q).ints=current_int;
int_start=0;
end

linked_int_durations(q)=linked_int_durations(q)+joining_int_duration;
linked_ints(q).ints=[linked_ints(q).ints joining_int];

link_number=link_number+1;

current_int_index=joining_int_index;
else
current_int_index=ints_num+1;
end

end

end



%%%%%%%%%
%calculate the total duration of non linked interactions
%THE PARAMETER WILL BE MINIMIZED TO FIND THE BEST SOLUTION
%find interactions that are neither linked forwards nor backwards
%(these interactions should have a column with zeros and a row with
%zeros!!!!)

non_linking_ints_array=or(sum(int_linking_matrix_masked,1),sum(int_linking_matrix_masked,2)'); %if zero, the int is not linked!
non_linking_ints_indexes=find(non_linking_ints_array==0);

nonlinked_int_durations=[];

for nl_int=non_linking_ints_indexes
non_linking_int=ints_in_family(nl_int); %the real non-linking interaction number
nonlinked_int_durations(r)=int_matrix(non_linking_int,2)-int_matrix(non_linking_int,1)+1; %the duration of the non-linking interaction

nonlinked_ints(r)=non_linking_int;

r=r+1;
end

else
linked_int_durations=NaN; %in case of circularly linked objects!!!!
nonlinked_int_durations=NaN; %idem
linked_ints=NaN;
nonlinked_ints=NaN;
end
end
