% Copyright (c) Lanoiselée and Calebiro 2023
% Please cite:
% Plasma membrane preassociation drives beta-arrestin coupling to receptors and activation
% Grimes, Koszegi, Lanoiselée et al.
% Cell 186, 1-18 (2023)
% https://doi.org/10.1016/j.cell.2023.04.018
% 
% This script is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY;
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
function [list_state] = add_list_inmask(X,Y,stack_binary,state_name,list_state)

% colocalization of trajectories with the binary mask each frame
% result is appended in the field list_state.(state_name)
% binary mask can either be 1 frames or a stack with one image per frame
%% round trajectory coordinates
int_x=floor(X);
int_y=floor(Y);

%% Pre-Allocate memory for list in Mask
list_in_msk=nan(size(list_state.exist));
%% loop for aligning Mask image and calculating colocalisation at each frame
if size(stack_binary,3)==1
list_in_msk(list_state.exist(:))=double(stack_binary(sub2ind(size(stack_binary),int_y(list_state.exist(:)),int_x(list_state.exist(:)))));
list_in_msk(~list_state.exist==1)=nan;
else
size_im=size(stack_binary(:,:,1));
for frame=1:size(int_x,2)
im_to_check=stack_binary(:,:,frame);
%% localisation on CCP: 1 if on CCP, 0 if outside
list_in_msk(list_state.exist(:,frame),frame)=double(im_to_check(sub2ind(size(im_to_check),int_y(list_state.exist(:,frame),frame),int_x(list_state.exist(:,frame),frame))));
end
end
list_state.(state_name)=list_in_msk;
end

