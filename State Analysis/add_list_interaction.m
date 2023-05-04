% Copyright (c) Lanoiselée and Calebiro 2023
% Please cite:
% Plasma membrane preassociation drives beta-arrestin coupling to receptors and activation
% Grimes, Koszegi, Lanoiselée et al.
% Cell 186, 1-18 (2023)
% https://doi.org/10.1016/j.cell.2023.04.018
% 
% This script is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY;
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
function [list_state] = add_list_interaction(X,interaction_data,channel_number,couple_int,list_state)
% function that reads the interaction data and generate a
% binary list_state.interaction field where 1 where the molecule is colocalising 
% another list_state.interaction_partner field is generated that stores
% with which molecule in the other channel is colocalizing.
couple_int=sort(couple_int);
channel={'C1','C2','C3','C4'};
%%
list_interact_ch1=zeros(size(X{1},1),size(X{1},2));
list_interact_ch2=zeros(size(X{2},1),size(X{2},2));

list_interact_all_partner_ch1=cell(size(X{1},1),size(X{1},2));
list_interact_all_partner_ch2=cell(size(X{2},1),size(X{2},2));
list_compound_idx_ch1=cell(size(X{1},1),size(X{1},2));
list_compound_idx_ch2=cell(size(X{2},1),size(X{2},2));
%% run through trajectories
N=size(X{1},2);
Long_int_list=1:numel(interaction_data.list_compound_int);
for ch=couple_int
interaction_data.list_compound_durations=interaction_data.list_compound_end-interaction_data.list_compound_start+1;
end

%% search interactions
for n_comp=1:numel(Long_int_list)
compInt=Long_int_list(n_comp);
%%% define trajectory before and after interaction
%% pre_trj
[~,int_start_idx]=min(interaction_data.int_matrix(interaction_data.list_compound_int{compInt},1));
int_numb_start=interaction_data.list_compound_int{compInt}(int_start_idx); %% the first 'interaction number' in 'compound_ints'
compound_start_frame=interaction_data.int_matrix(int_numb_start,1);                           %% frame at which interaction starts
preObj_ch1=interaction_data.int_matrix(int_numb_start,4);                        %% Ch1 object at which interaction starts
preObj_ch2=interaction_data.int_matrix(int_numb_start,8);                        %% Ch2 object at which interaction starts

%% post_trj
%int_numb_end=max(interaction_data.list_compound_int(1,compInt).ints);   %% the last 'interaction number' in 'compound_ints'
[~,int_end_idx]=max(interaction_data.int_matrix(interaction_data.list_compound_int{compInt},2));
int_numb_end=interaction_data.list_compound_int{compInt}(int_end_idx); %% the first 'interaction number' in 'compound_ints'

compound_end_frame=interaction_data.int_matrix(int_numb_end,2);                          %% frame at which interaction ends
postObj_ch1=interaction_data.int_matrix(int_numb_end,4);                         %% Ch1 object at which interaction ends
postObj_ch2=interaction_data.int_matrix(int_numb_end,8);                         %% Ch2 object at which interaction ends

%% write 1 when the particle is interacting
list_interact_ch1(preObj_ch1,compound_start_frame:compound_end_frame)=1;
list_interact_ch2(preObj_ch2,compound_start_frame:compound_end_frame)=1;
%% write the idx of partners on the other channel when the particle is colocalising
list_interact_all_partner_ch1(preObj_ch1,compound_start_frame:compound_end_frame)=cellfun(@(x)[x,preObj_ch2], list_interact_all_partner_ch1(preObj_ch1,compound_start_frame:compound_end_frame),'UniformOutput',false);
list_interact_all_partner_ch2(preObj_ch2,compound_start_frame:compound_end_frame)=cellfun(@(x)[x,preObj_ch1], list_interact_all_partner_ch2(preObj_ch2,compound_start_frame:compound_end_frame),'UniformOutput',false);

list_compound_idx_ch1(preObj_ch1,compound_start_frame:compound_end_frame)=cellfun(@(x)[x,n_comp],list_compound_idx_ch1(preObj_ch1,compound_start_frame:compound_end_frame),'UniformOutput',false);
list_compound_idx_ch2(preObj_ch2,compound_start_frame:compound_end_frame)=cellfun(@(x)[x,n_comp],list_compound_idx_ch2(preObj_ch2,compound_start_frame:compound_end_frame),'UniformOutput',false);
end

%% reordering the interaction partners from longest to shortest interaction compound
list_interact_all_partner_ch1=cellfun(@(x,y) order_by_compound_duration(x,y,interaction_data.list_compound_durations), list_interact_all_partner_ch1,list_compound_idx_ch1,'UniformOutput',false);
list_interact_all_partner_ch2=cellfun(@(x,y) order_by_compound_duration(x,y,interaction_data.list_compound_durations), list_interact_all_partner_ch2,list_compound_idx_ch2,'UniformOutput',false);
list_interact_partner_ch1=cellfun(@(x)find_first_or_empty(x),list_interact_all_partner_ch1);%,'UniformOutput',false);
list_interact_partner_ch2=cellfun(@(x)find_first_or_empty(x),list_interact_all_partner_ch2);%,'UniformOutput',false);

if channel_number==couple_int(1)
list_state.interaction{couple_int(1),couple_int(2)}=list_interact_ch1;
list_state.interaction_partner{couple_int(1),couple_int(2)}=list_interact_partner_ch1;
elseif channel_number==couple_int(2)
list_state.interaction{couple_int(1),couple_int(2)}=list_interact_ch2;
list_state.interaction_partner{couple_int(1),couple_int(2)}=list_interact_partner_ch2;
end

end

function [interact_partner_ordered]=order_by_compound_duration(interact_partner,compound_list,compound_duration_list)
if numel(compound_list)>1
[~,idx_sort]=sort(compound_duration_list(compound_list),'descend');
interact_partner_ordered=interact_partner(idx_sort);
else
interact_partner_ordered=interact_partner;    
end
end

function result=find_first_or_empty(x)
if ~isempty(x)
result=x(1);
else
result=nan;
end
end
