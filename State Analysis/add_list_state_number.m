% Copyright (c) Lanoiselée and Calebiro 2023
% Please cite:
% Plasma membrane preassociation drives beta-arrestin coupling to receptors and activation
% Grimes, Koszegi, Lanoiselée et al.
% Cell 186, 1-18 (2023)
% https://doi.org/10.1016/j.cell.2023.04.018
% 
% This script is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY;
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
function [list_state] = add_list_state_number(list_state_ch,channel_number,couple_int,parameter)
% Generates a field 'list_state.state_number' where for each trajectory and
% each frame a state number is associated depending on the binary vectors
% from 'parameter.state_vec_fields' to be considered. 
channel={'C1','C2','C3','C4'};
n_state=size(parameter.state,1);
if ~isempty(couple_int)
    list_channel=couple_int;
    if channel_number==couple_int(1)
        other_channel=couple_int(2);
    elseif channel_number==couple_int(2)
        other_channel=couple_int(1);
    end
else
    list_channel=channel_number;
end
interaction_pos=find(ismember(parameter.state_vec_fields,'interaction'));

for ch=list_channel
    %% store binary_vec in list for both channels
    for n_state_dim=1:numel(parameter.state_vec_fields)
        inter=list_state_ch{ch}.binary_vec(:,:,n_state_dim);
        list{1,ch}(:,n_state_dim)=reshape(double(inter)',[],1);
    end
    if ~isempty(couple_int)
        if ch==channel_number
            M=size(list_state_ch{ch}.on_membrane,1);
            N=size(list_state_ch{ch}.on_membrane,2);
            mem_keep=list_state_ch{ch}.on_membrane;
        elseif ch==other_channel
            M_other_channel=size(list_state_ch{ch}.on_membrane,1);
            N_other_channel=size(list_state_ch{ch}.on_membrane,2);
        end
    else
        M=size(list_state_ch{ch}.on_membrane,1);
        N=size(list_state_ch{ch}.on_membrane,2);
        mem_keep=list_state_ch{ch}.on_membrane;
    end 
end
%% translating interaction partner from subscript to indices
if ~isempty(couple_int)
    test= list_state_ch{channel_number}.interaction_partner{couple_int(1),couple_int(2)};
    test2=N*(test-1)+repmat(1:N_other_channel,M,1);
    list_interaction_partner_ind=reshape(test2',[],1);
end
list_state_number=nan(M,N);
list_state_number_ind=reshape( list_state_number',[],1);

%Order of states follow the order from parameter.state_vec_fields
List_N=1:size(list{1,channel_number},1);
List_N=List_N(1,list{1,channel_number}(:,interaction_pos)==1);
%% get state vector of partner for interaction
list_not_exist_sub=isnan(mem_keep);
list_not_exist_ind2=reshape(list_not_exist_sub',[],1);
list_state_number_inter=nan(sum(~list_not_exist_ind2),1);
list{1,channel_number}(list_not_exist_ind2,:)=[];
if ~isempty(couple_int)
    state_partner=nan(N*M,numel(parameter.state_vec_fields));
    list_int_not_exist_ind=isnan(list_interaction_partner_ind);
    list_interaction_partner_ind(list_int_not_exist_ind)=[];
    state_partner(List_N',:)=list{1,other_channel}(list_interaction_partner_ind,:);
    state_partner(list_not_exist_ind2,:)=[];
end

%%  set state number based on combinations of binary vectors defined in 'parameter.state'
for state_num=1:n_state
    for  n_sub_state=1:numel(parameter.state{state_num,2})  
        if ~isempty(couple_int)
            if parameter.state{state_num,2}{n_sub_state}(interaction_pos,1)==1 %&& ~isempty(List_N)
                list_found_ch1=ismember(list{1,channel_number},parameter.state{state_num,2}{n_sub_state}','rows');
                list_found_ch2=ismember(state_partner,parameter.state{state_num,3}{n_sub_state}','rows');
                list_state_number_inter(list_found_ch1 & list_found_ch2,1)=state_num;
            else
                list_found=ismember(list{1,channel_number},parameter.state{state_num,2}{n_sub_state}','rows');
                list_state_number_inter(list_found'==1,1)=state_num;
            end
        else
            list_found=ismember(list{1,channel_number},parameter.state{state_num,2}{n_sub_state}','rows');
            list_state_number_inter(list_found'==1,1)=state_num;
        end
    end
end
list_state_number_ind(~list_not_exist_ind2)=list_state_number_inter;
list_state_number=reshape(list_state_number_ind,N,M)';
list_state=list_state_ch{channel_number};
list_state.state_number=list_state_number;
end
