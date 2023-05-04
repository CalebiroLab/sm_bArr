% Copyright (c) Lanoiselée and Calebiro 2023
% Please cite:
% Plasma membrane preassociation drives beta-arrestin coupling to receptors and activation
% Grimes, Koszegi, Lanoiselée et al.
% Cell 186, 1-18 (2023)
% https://doi.org/10.1016/j.cell.2023.04.018
% 
% This script is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY;
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
function [list_state] = add_list_binary_vec(list_state,couple_int,parameter)
% generates a field 'list_state.binary_vec' where the binary fields are
% combined according to the order specified in 'parameter.state_vec_fields'
list=zeros(size(list_state.exist,1),size(list_state.exist,2),numel(parameter.state_vec_fields)-1);
for n_dim=1:numel(parameter.state_vec_fields)-1
    if strcmp(parameter.state_vec_fields{n_dim},'interaction')
        list(:,:,n_dim)=double(list_state.(parameter.state_vec_fields{n_dim}){couple_int(1),couple_int(2)});
    else
        list(:,:,n_dim)=double(list_state.(parameter.state_vec_fields{n_dim}));
    end
end
list_state.binary_vec=list;
end
