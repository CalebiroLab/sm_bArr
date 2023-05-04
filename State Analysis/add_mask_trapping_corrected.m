% Copyright (c) Lanoiselée and Calebiro 2023
% Please cite:
% Plasma membrane preassociation drives beta-arrestin coupling to receptors and activation
% Grimes, Koszegi, Lanoiselée et al.
% Cell 186, 1-18 (2023)
% https://doi.org/10.1016/j.cell.2023.04.018
% 
% This script is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY;
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
function [list_state] =  add_mask_trapping_corrected(list_state,state_name)
%Correction of trapping at CCP
% turns fluctuation in/out cccp while confined that are dure to localisation error to only trapped in CCP
% generates a field 'list_state.([state_name,'_trapping_corrected'])' that
% is equal to one when the trapping in the given domains have been
% corrected
%%
for m=1:size(list_state.trapped,1)
    [ List_min_max_trapped ] = Make_list_min_max_index_equal( list_state.trapped(m,:)',1 );
    if ~isempty(List_min_max_trapped)
        for kj=1:size(List_min_max_trapped,1)
            [ List_min_max_mask ] = Make_list_min_max_index_equal( list_state.(state_name)(m,List_min_max_trapped(kj,1):List_min_max_trapped(kj,2))',1 );
            if size(List_min_max_mask,1)>1
                interval=List_min_max_trapped(kj,1)-1+(List_min_max_mask(1,1):List_min_max_mask(end,2));
                list_state.(state_name)(m,interval)=1;
            end
        end
    end
end
list_state.([state_name,'_trapping_corrected'])=1;
end

