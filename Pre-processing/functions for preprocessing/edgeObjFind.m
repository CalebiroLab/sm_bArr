
function [exclobj,OBJ] = edgeObjFind (OBJ, IFO, im, plotting)
% To be used with Till's array.
% This function identify the objects at the edges.
% The output is an array with all the excluded objects. And plot the position of the objects. 


a=[];
k= ((IFO.roiSize-1)/2)+(1);
for n= 1: size(OBJ.mxR,1);
   if round(min(OBJ.mxR(n,:))) - k <=0 || round(max(OBJ.mxR(n,:))) + k>=size(im,2) ||... 
      round(min(OBJ.myR(n,:))) - k<=0 || round(max(OBJ.myR(n,:))) + k>=size(im,1);
       a=[a,n];
   end
end
exclobj=a;
for i = 1:numel(exclobj);
    OBJ.trjRB(exclobj(i),:) = NaN;
end

%%
if plotting ==1
 figure ('name', 'excluded objects');
imagesc(im);
colormap(gray);
hold on 
for i=1:numel(exclobj);
    plot(OBJ.mxR(exclobj(i),1),OBJ.myR(exclobj(i),1),'oy');
    text(OBJ.mxR(exclobj(i),1)+3,OBJ.myR(exclobj(i),1), ['\color{yellow}' num2str(exclobj(i))])
end
end
