% Copyright (c) Lanoiselée and Calebiro 2023
% Please cite:
% Plasma membrane preassociation drives beta-arrestin coupling to receptors and activation
% Grimes, Koszegi, Lanoiselée et al.
% Cell 186, 1-18 (2023)
% https://doi.org/10.1016/j.cell.2023.04.018
% 
% This script is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY;
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
function [list_state] = add_markov_chain_forward_supdur(list_state,parameter)
% compute the transition in the forward markov chain defined by the states defined
% in the parameter file
% if and only if either the current of next state are longer than
% 'parameter.sup_dur'
% Result is stored in the field 'list_state.markov_chain_forward'

n_state=size(parameter.state,1);
list_state_number_inter=list_state.state_number;
list_state_number_inter=list_state_number_inter(:,1:parameter.markov_chain_Nmax);
M=size(list_state_number_inter,1);
Mat=zeros(n_state,n_state);
for m=1:M
    [state_vector_unique] = make_state_unique_nrep(list_state_number_inter(m,:));
    if size(state_vector_unique)>=2
        for n=2:size(state_vector_unique,2)
            from=state_vector_unique(1,n-1);
            to=state_vector_unique(1,n);
            
            %% if either before or after state are longer tha parameter.sup_dur
            if isnan(from)==0 && isnan(to)==0 ...
                    && ((state_vector_unique(2,n)>=parameter.sup_dur || state_vector_unique(2,n-1)>=parameter.sup_dur)...
                    ||(from ==1 && state_vector_unique(2,n)>=parameter.sup_dur)...
                    ||(to==1 && state_vector_unique(2,n-1)>=parameter.sup_dur))
                
                Mat(from,from)=Mat(from,from)+state_vector_unique(2,n-1);
                Mat(from,to)=Mat(from,to)+1;
            end
        end
    end
end
list_state.markov_chain_forward_supdur=Mat;
end
