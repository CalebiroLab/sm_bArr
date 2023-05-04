% Copyright (c) Lanoiselée and Calebiro 2023
% Please cite:
% Plasma membrane preassociation drives beta-arrestin coupling to receptors and activation
% Grimes, Koszegi, Lanoiselée et al.
% Cell 186, 1-18 (2023)
% https://doi.org/10.1016/j.cell.2023.04.018
% 
% This script is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY;
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
function [list_state] = add_markov_chain_forward(list_state,parameter)
% compute the transition in the forward markov chain defined by the states defined
% in the parameter file
% Result is stored in the field 'list_state.markov_chain_forward'
n_state=size(parameter.state,1);
list=list_state.state_number;

M=size(list,1);
Mat=zeros(n_state,n_state);
for m=1:M
    for n=2:parameter.markov_chain_Nmax
        if isnan(list_state.state_number(m,n,1))==0 && isnan(list_state.state_number(m,n-1,1))==0
            from=list_state.state_number(m,n-1);
            to=list_state.state_number(m,n);
            
            Mat(from,to)=Mat(from,to)+1;
        end
    end
end
list_state.markov_chain_forward=Mat;
end
