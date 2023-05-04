% Copyright (c) Lanoiselée and Calebiro 2023
% Please cite:
% Plasma membrane preassociation drives beta-arrestin coupling to receptors and activation
% Grimes, Koszegi, Lanoiselée et al.
% Cell 186, 1-18 (2023)
% https://doi.org/10.1016/j.cell.2023.04.018
% 
% This script is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY;
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
function [list_state] = add_list_istrapped(X,Y,list_state,parameter)
% Script that detects transient trapping of each trajectories using the
% code from
% Detecting Transient Trapping from a Single Trajectory: A Structural Approach
% Y Lanoiselée, J Grimes, Z Koszegi, D Calebiro
% Entropy 23 (8), 1044 (2021)

% Result is stored in the field 'list_state.trapped'
dataset.X=X;
dataset.Y=Y;
[list_trapped] = Detect_transient_trapping_multiscale(dataset,parameter);
list_state.trapped=list_trapped;
end
