% Copyright (c) Lanoiselée and Calebiro 2023
% Please cite:
% Plasma membrane preassociation drives beta-arrestin coupling to receptors and activation
% Grimes, Koszegi, Lanoiselée et al.
% Cell 186, 1-18 (2023)
% https://doi.org/10.1016/j.cell.2023.04.018
% 
% This script is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY;
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
function [X,Y]=gap_close(X,Y,method)
if nargin==2
    method='linear_interpolation'; % default is linear interpolation
end
for channel=1:numel(X)
    M=size(X{channel},1);
    list_exist=~isnan(X{channel});
    for m=1:M
        pos_beginning=find(list_exist(m,:),1);
        if ~isempty(pos_beginning)
            [ List_min_max ] = list_min_max_index_equal( list_exist(m,pos_beginning:end)',1 );
            if size(List_min_max,1)>1
                for n_gap=1:size(List_min_max,1)-1
                    gap_fstart=List_min_max(n_gap,2)+1+pos_beginning-1;
                    gap_fend=List_min_max(n_gap+1,1)-1+pos_beginning-1;
                    if strcmp(method,'linear_interpolation')
                        X{channel}(m,gap_fstart-1:gap_fend+1)=linspace(X{channel}(m,gap_fstart-1),X{channel}(m,gap_fend+1),gap_fend-gap_fstart+3);
                        Y{channel}(m,gap_fstart-1:gap_fend+1)=linspace(Y{channel}(m,gap_fstart-1),Y{channel}(m,gap_fend+1),gap_fend-gap_fstart+3);
                    elseif strcmp(method,'brownian_bridge')
                        driftX=linspace(X{channel}(m,gap_fstart-1),X{channel}(m,gap_fend+1),gap_fend-gap_fstart+3);
                        driftY=linspace(Y{channel}(m,gap_fstart-1),Y{channel}(m,gap_fend+1),gap_fend-gap_fstart+3);
                        sigx=std(diff(X{channel}(m,:)),'omitnan');
                        sigy=std(diff(Y{channel}(m,:)),'omitnan');
                        N=gap_fend-gap_fstart+2;
                        W=[0,cumsum(sigx*randn(1,N))];
                        bridgex=W-(0:N)/N.*W+driftX;
                        W=[0,cumsum(sigy*randn(1,N))];
                        bridgey=W-(0:N)/N.*W+driftY;
                        X{channel}(m,gap_fstart-1:gap_fend+1)=bridgex;
                        Y{channel}(m,gap_fstart-1:gap_fend+1)=bridgey;
                    else
                        disp('Wrong method name')
                        break
                    end
                end
            end
        end
    end
end
end


