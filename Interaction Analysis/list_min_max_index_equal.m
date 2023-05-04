% Copyright (c) Lanoiselée and Calebiro 2023
% Please cite:
% Plasma membrane preassociation drives beta-arrestin coupling to receptors and activation
% Grimes, Koszegi, Lanoiselée et al.
% Cell 186, 1-18 (2023)
% https://doi.org/10.1016/j.cell.2023.04.018
% 
% This script is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY;
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
function [ List_min_max ] = list_min_max_index_equal( Vector,quantity )
% Return return a 2d vector with in first column all min and in the secon
% all max where vector_to_test consecutively equals quantity
% quantity default value is 1
if nargin==1
    quantity=1;
end
N_vec=size(Vector,1);
index=1;
List_min_max=[];
stop=0;
while stop==0
    if Vector(index,1)==quantity
        [pos_down,pos_up] = find_min_max_index_equal( index,Vector,quantity);
        List_min_max=[List_min_max;pos_down,pos_up];
        index=pos_up+1;
    else
        index=index+1;
    end
    if index>N_vec
        stop=1;
    end
end

    function [pos_down,pos_up] = find_min_max_index_equal( index,Vector_to_test,quantity)
        % Return min and max index equal of result where vector_to_test==quantity
        % quantity default value is 1
        % 0 0 0    1      1   1   1 1   1 0  0 0
        %       pos_down    index    pos_up
        N=size(Vector_to_test,1);
        if nargin==2
            quantity=1;
        end
        stop_up=0;
        stop_down=0;
        pos_up=index;
        pos_down=index;
        for z=1:N
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
            if stop_down==0
                if pos_down-1>=1
                    if Vector_to_test(pos_down-1,1)==quantity
                        pos_down=pos_down-1;
                    else
                        stop_down=1;
                    end
                else
                    stop_down=1;
                end
            end
            if stop_up==1 && stop_down==1
                break
            end
        end
    end
end

