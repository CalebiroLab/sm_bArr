% Copyright (c) Lanoiselée and Calebiro 2023
% Please cite:
% Plasma membrane preassociation drives beta-arrestin coupling to receptors and activation
% Grimes, Koszegi, Lanoiselée et al.
% Cell 186, 1-18 (2023)
% https://doi.org/10.1016/j.cell.2023.04.018
% 
% This script is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY;
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
function [list_state] = add_list_exist(X,Y,list_state,parameter)
% generates a binary field 'list_state.exist' with values 1 where trajectories coordinates are not nan 
% and trajectories shorter than 'parameter.min_traj_length' are considered
% not existing (value=0)
list_exist=~isnan(X) & ~isnan(Y);
list_exist(sum(double(list_exist),2)<parameter.min_traj_length,:)=0;
list_state.exist=list_exist;

end
