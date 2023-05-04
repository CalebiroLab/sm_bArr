% Copyright (c) Lanoiselée and Calebiro 2023
% Please cite:
% Plasma membrane preassociation drives beta-arrestin coupling to receptors and activation
% Grimes, Koszegi, Lanoiselée et al.
% Cell 186, 1-18 (2023)
% https://doi.org/10.1016/j.cell.2023.04.018
% 
% This script is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY;
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
function [int_linking_matrix,int_linking_matrix_durations,int_families] = interaction_linking(int_matrix,gap_wind)
%from interaction matrix generates int_linking_matrix,int_linking_matrix_durations,int_families
int_matrix=sparse(int_matrix);
%% computation of int_linking_matrix
int_linking_matrix=(((int_matrix(:,4)==int_matrix(:,4)')|((int_matrix(:,8)==int_matrix(:,8)')))&((int_matrix(:,1)>= (int_matrix(:,2)-gap_wind+1)')&(int_matrix(:,1)<= (int_matrix(:,2)+gap_wind+1)')&(int_matrix(:,2)>int_matrix(:,2)')));

diag_fill=zeros(size(int_matrix,1),1);
int_linking_matrix=spdiags(diag_fill,0,int_linking_matrix);
[int1,int2]=find(int_linking_matrix==1); % find connected interaction pairs

%% compute duration of combined pairs of interactions
int_linking_matrix_durations=sparse(size(int_linking_matrix,1),size(int_linking_matrix,2));% init int_linking_duration
int_linking_matrix_durations(sub2ind(size(int_linking_matrix_durations),int1,int2))=(int_matrix(int1,2)-int_matrix(int1,1)+1)+(int_matrix(int2,2)-int_matrix(int2,1)+1);% int_linking_durations are duration of interaction 1 + duration of interaction int2

%% find interaction families (connected subgraphs) in the  int_linking_matrix (symetrised transition matrix with added diagonal terms)
connectivity_matrix=int_linking_matrix+int_linking_matrix'+speye(size(int_linking_matrix,1));
idx_int_family = conncomp(graph(connectivity_matrix));%find all disconencted components
list_fam_unique=unique(idx_int_family); % unique family number occurences

list_int_idx=1:numel(idx_int_family ); % all interaction IDs
int_families=[]; % init int_families
%% add interaction families
for n_fam=list_fam_unique
int_families(n_fam).ints=sort(list_int_idx(idx_int_family==n_fam));
end
end

