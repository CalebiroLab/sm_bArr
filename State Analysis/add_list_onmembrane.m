% Copyright (c) Lanoiselée and Calebiro 2023
% Please cite:
% Plasma membrane preassociation drives beta-arrestin coupling to receptors and activation
% Grimes, Koszegi, Lanoiselée et al.
% Cell 186, 1-18 (2023)
% https://doi.org/10.1016/j.cell.2023.04.018
% 
% This script is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY;
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
function [list_state] = add_list_onmembrane(list_state,couple_int,parameter)
% Generates a field 'list_state.on_membrane' with 1 when molecule is
% visible, -1 at the frame just before and after being visible and 0
% otherwise

M=size(list_state.exist,1);
N=size(list_state.exist,2);
on_membrane=nan(size(list_state.exist));
if any( strcmp(parameter.state_vec_fields,'interaction'))
    on_membrane(list_state.exist==1 | list_state.interaction{couple_int(1),couple_int(2)}==1)=1;
else
    on_membrane(list_state.exist==1)=1;
end

for m=1:M
    if any(list_state.exist(m,:)==1)
        [ List_min_max ] = Make_list_min_max_index_equal(on_membrane(m,:)'==1,1 );
        if ~isempty(List_min_max)
            for nom=1:size(List_min_max,1)
                if List_min_max(nom,1)>1
                    
                    on_membrane(m,List_min_max(nom,1)-1)=0;
                    list_state.binary_vec(m,List_min_max(nom,1)-1,:)=-1;
                end
                if List_min_max(nom,2)<N
                    on_membrane(m,List_min_max(nom,2)+1)=0;
                    list_state.binary_vec(m,List_min_max(nom,2)+1,:)=-1;
                end
            end
        end
    end
end

list_state.binary_vec(:,:,numel(parameter.state_vec_fields))=on_membrane;
exist_in_cell=list_state.exist;
exist_in_cell(on_membrane==0)=1;
list_state.exist_in_cell=exist_in_cell;
list_state.on_membrane=on_membrane;

end
