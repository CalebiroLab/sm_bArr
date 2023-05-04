% Copyright (c) Lanoiselée and Calebiro 2023
% Please cite:
% Plasma membrane preassociation drives beta-arrestin coupling to receptors and activation
% Grimes, Koszegi, Lanoiselée et al.
% Cell 186, 1-18 (2023)
% https://doi.org/10.1016/j.cell.2023.04.018
% 
% This script is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY;
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
function [state_vector_unique] = make_state_unique_nrep(state_vector)
% Transform a row vector of state numbers in a two rows vector
% first vector contains states within the vector wihtout consecutive
% repetitions (except for nan)
% second vector is the number of repetition of the given state

index=1;
state_vector_unique=nan(2,numel(state_vector)); % pre-allocate memory for the worst case possible (no consecutive states)
cp=0;
while index<=size(state_vector,2)
    cp=cp+1;
    [pos_up] = find_index_equal_wrap_right( index,state_vector',state_vector(1,index));
    state_vector_unique(:,cp)=[state_vector(1,index);pos_up-index+1];
    index=pos_up+1;
end

if cp<numel(state_vector)
    state_vector_unique(:,cp+1:end)=[]; %% remove entries that were not used.
end

    function [pos_up] = find_index_equal_wrap_right( index,Vector_to_test,quantity)
        % Return min and max index equal of result where vector_to_test==quantity
        % quantity default value is 1
        % 0 0 0    1      1   1   1 1   1 0  0 0
        %                   index    pos_up
        %Testing continuous path
        N=size(Vector_to_test,1);
        
        stop_up=0;
        pos_up=index;
        while stop_up==0
            if stop_up==0
                if pos_up+1<=N
                    if Vector_to_test(pos_up+1,1)==quantity
                        pos_up=pos_up+1;
                    else
                        stop_up=1;
                    end
                else
                    stop_up=1;
                end
            end
            
        end
    end
end

