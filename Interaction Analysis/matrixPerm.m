% Copyright (c) Lanoiselée and Calebiro 2023
% Please cite:
% Plasma membrane preassociation drives beta-arrestin coupling to receptors and activation
% Grimes, Koszegi, Lanoiselée et al.
% Cell 186, 1-18 (2023)
% https://doi.org/10.1016/j.cell.2023.04.018
% 
% This script is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY;
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
function [matrix_out] = matrixPerm(matrix_in_1,matrix_in_2,biunivocal)
%this routine replaces one value of matrix_in_1 with one random value taken
%from matix_in_2. If biunivocal==1, it reiterates until a set is found
%for which each raw and each column contain no more than one TRUE (i.e. each
%link is either biunivocal or the int links to none).
found=0;
[row,col]=find(matrix_in_2>0);
n_row=size(row,1);
while found==0 && n_row>0

n=size(matrix_in_1,2); %the number of columns

pos=randi(n_row);
j=row(pos);
i=col(pos);

matrix_out=matrix_in_1;
matrix_out(:,i)=0;
matrix_out(j,:)=0;
matrix_out(j,i)=1;

if biunivocal==0
found=1;
else
if any(sum(matrix_out,1)>1)==0 && any(sum(matrix_out,2)>1)==0 %i.e. each link is either biunivocal or the int links to none
found=1;
else
row(pos)=[];
col(pos)=[];
n_row=n_row-1;
end
end

end

end

