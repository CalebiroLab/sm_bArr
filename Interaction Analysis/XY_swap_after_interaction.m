% Copyright (c) Lanoiselée and Calebiro 2023
% Please cite:
% Plasma membrane preassociation drives beta-arrestin coupling to receptors and activation
% Grimes, Koszegi, Lanoiselée et al.
% Cell 186, 1-18 (2023)
% https://doi.org/10.1016/j.cell.2023.04.018
% 
% This script is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY;
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
function [X,Y,int_matrix,interaction_data] = XY_swap_after_interaction(X,Y,int_matrix,interaction_data)
% script for creating new OBJ from interaction files
%% add interaction idx at 11th column
int_matrix(:,11)=ones(size(int_matrix,1),1);

%% add binary for being considered in swaping at 12th column
int_matrix(:,12)=ones(size(int_matrix,1),1);
%% use the longest interaction compound for the swap
list_mol_ch=cell(1,2);
for ch=1:2
    list_mol_ch{ch}=unique(int_matrix(:,4*ch));
    list_mol_ch{ch}(isnan(list_mol_ch{ch}))=[];
    
    % list all interactions with similar numbers in ch
    location_ch{ch}=int_matrix(:,4*ch)==list_mol_ch{ch}';%n-th row is n-th coloc m-th column is m-th molecule 
    
end
%
list_compound_int=struct2cell(interaction_data.list_compound_int(:));
list_compound_int_new=list_compound_int;

%% append int_matrix with compound number in 13th column
int_matrix(:,13)=nan(size(int_matrix,1),1);
for n_compound_copy=1:numel(list_compound_int)
    int_matrix(list_compound_int{1,n_compound_copy},13)=n_compound_copy;
end

%%
list_compound_dur=(interaction_data.list_compound_end-interaction_data.list_compound_start+1)';
n_comp=size(list_compound_int,2); %number of compound interactions
for ch=1:2
    %% for ch1
    for m=1:numel(list_mol_ch{ch})
        list_beg_end_int=[int_matrix(location_ch{ch}(:,m),1),int_matrix(location_ch{ch}(:,m),2)];
        
        % test if two interactions with the same molecule do overlap
%         test=(((list_beg_end_int(:,1)<=list_beg_end_int(:,1)') & (list_beg_end_int(:,2)>list_beg_end_int(:,1)') & (list_beg_end_int(:,2)>=list_beg_end_int(:,1)'))...
%             |((list_beg_end_int(:,1)<list_beg_end_int(:,2)') & (list_beg_end_int(:,2)>list_beg_end_int(:,2)')))...
%             |...
%             (((list_beg_end_int(:,1)<=list_beg_end_int(:,1)') & (list_beg_end_int(:,2)>list_beg_end_int(:,1)'))...
%             &((list_beg_end_int(:,1)<=list_beg_end_int(:,2)') & (list_beg_end_int(:,2)>list_beg_end_int(:,2)')));
%         
        test=(list_beg_end_int(:,1)<=list_beg_end_int(:,2)') & (list_beg_end_int(:,2)>=list_beg_end_int(:,1)');
        
        
        
        test=test & logical((1-tril(ones(numel(list_beg_end_int(:,1))))));% remove lower triangle including diagonal
        if any(test(:))
            ID_mol=unique(int_matrix(location_ch{ch}(:,m),4*ch));
            ID_ints=find(location_ch{ch}(:,m)==1)';
            
            [row,col]=find(test);
            ID_int_conflict=[ID_ints(1,row)',ID_ints(1,col)'];
            for n_conflict=1:size(ID_int_conflict,1)
                ID_compound_int1=int_matrix(ID_int_conflict(n_conflict,1),13);
                ID_compound_int2=int_matrix(ID_int_conflict(n_conflict,2),13);
                % if entries are not found in compounds then longest between the two interactions is chosen
                if isnan(ID_compound_int1)
                    dur1=int_matrix(ID_int_conflict(n_conflict,1),2)-int_matrix(ID_int_conflict(n_conflict,1),1)+1;
                else
                    dur1=list_compound_dur(1,ID_compound_int1);
                end
                if isnan(ID_compound_int2)
                    dur2=int_matrix(ID_int_conflict(n_conflict,2),2)-int_matrix(ID_int_conflict(n_conflict,2),1)+1;
                else
                    dur2=list_compound_dur(ID_compound_int2);
                end
                
                if dur1>dur2
                    int_matrix(ID_int_conflict(n_conflict,2),12)=0;
                elseif dur2>=dur1
                    int_matrix(ID_int_conflict(n_conflict,1),12)=0;
                end
                
            end
        end
    end
    
end
%% create and sort compound interaction list

list_last_val=zeros(n_comp,1);% list of last frame number of latest coloc
for n=1:n_comp
    list_last_val(n,1)=max(int_matrix(list_compound_int{1,n},2));% last frame number of latest colocalisation
end
[sorted_list_val,order]=sort(list_last_val,'descend');

list_compound_int_ord=list_compound_int(order);
list_compound_int_ord_new=list_compound_int_ord;

int_matrix_new=int_matrix;
for n=1:n_comp
    num_interaction_comp=numel(list_compound_int_ord{1,n}); %number of interactions in the compound
    if  num_interaction_comp>1
        for n_int_con=num_interaction_comp-1:-1:1
            
            int_matrix_mem=int_matrix_new;
            int_current=list_compound_int_ord{1,n}(1,n_int_con);
            int_next= list_compound_int_ord{1,n}(1,n_int_con+1);
            ID_int_part_current=int_matrix(int_current,[4,8]);
            ID_int_part_next=int_matrix(int_next,[4,8]);
            beg_next_int=int_matrix(int_next,1);
            end_next_int=int_matrix(int_next,2);
            
            
            if int_matrix_new(int_next,12)==1
                
                %% look at
                if (ID_int_part_current(1)-ID_int_part_next(1)==0) && (ID_int_part_current(2)-ID_int_part_next(2)~=0)
                    ch=2;
                elseif (ID_int_part_current(1)-ID_int_part_next(1)~=0) && (ID_int_part_current(2)-ID_int_part_next(2)==0)
                    ch=1;
                end
                % if particle ID in channel 1 is conserved but changed in channel 2
                %% Swap Coordinates
                % save traj coordinates of next int
                xch_portion_save=X{ch}(ID_int_part_next(ch),beg_next_int:end);
                ych_portion_save=Y{ch}(ID_int_part_next(ch),beg_next_int:end);
                % swap coordinates for the two sequences posterior to next int beginning
                X{ch}(ID_int_part_next(ch),beg_next_int:end)=X{ch}(ID_int_part_current(ch),beg_next_int:end);
                Y{ch}(ID_int_part_next(ch),beg_next_int:end)=Y{ch}(ID_int_part_current(ch),beg_next_int:end);
                X{ch}(ID_int_part_current(ch),beg_next_int:end)=xch_portion_save;
                Y{ch}(ID_int_part_current(ch),beg_next_int:end)=ych_portion_save;
                %                 %% Swap trjRB
                %                 % keep in memory traj coordinates of next int
                %                 trjRBch_portion_save=X{ch}(ID_int_part_next(ch),beg_next_int:end);
                %                 % swap coordinates for the two sequences posterior to next int beginning
                %                 trjRBch{ch}(ID_int_part_next(ch),beg_next_int:end)=trjRBch{ch}(ID_int_part_current(ch),beg_next_int:end);
                %                 trjRBch{ch}(ID_int_part_current(ch),beg_next_int:end)=trjRBch_portion_save;
                %% replace in int matrix ID of molecule
                int_matrix_new(int_matrix_mem(:,2)>=beg_next_int & int_matrix_mem(:,4*ch)==ID_int_part_next(ch),4*ch)=ID_int_part_current(ch);
                int_matrix_new(int_matrix_mem(:,2)>=beg_next_int & int_matrix_mem(:,4*ch)==ID_int_part_current(ch),4*ch)=ID_int_part_next(ch);
                int_matrix_mem=int_matrix_new;
            end
        end
    end
end

%% Creating new structure Ch

interaction_data.list_compound_int=list_compound_int_new;
interaction_data.list_compound_start=cellfun(@(x)min(int_matrix_new(x,1)),list_compound_int_new);
interaction_data.list_compound_end=cellfun(@(x)max(int_matrix_new(x,2)),list_compound_int_new);

int_matrix= int_matrix_new;


end

