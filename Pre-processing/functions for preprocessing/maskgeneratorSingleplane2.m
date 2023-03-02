% creation of the mask according to object positions and a given radius.
function [mask3,im] = maskgeneratorSingleplane2(im, xarray, yarray,iarray, radius, plotting)
% version where the image is already loaded
% correction to use images that are rectangular.

im = double(im);
mask = zeros(size(im,1), size(im,2));

for n=1:size(xarray,1);
    if ~isnan(iarray(n,1)) && xarray(n,1)>0 && yarray(n,1)>0;
    mask(yarray(n,1),xarray(n,1)) = 1;
    end
end
SE = strel('disk', radius);
mask2 = imdilate(mask,SE);

mask3 = abs(mask2-1);           %background areas are set to 1 and the areas of objects to 0

if plotting==1;
figure(1)
subplot(1,4,1), imagesc(im);
colormap 'gray';
hold on
for n=1:size(xarray,1);
    if ~isnan(iarray(n,1));
    plot(xarray(n,1),yarray(n,1), 'oy');
    end
end
hold off
title(['image ' num2str(1)]);

subplot(1,4,2), imagesc(mask2);
colormap 'gray'
title('Mask');

subplot(1,4,3), imagesc(mask2.*im);
colormap 'gray'
hold on
for n=1:size(xarray,1);
    if ~isnan(iarray(n,1));
    plot(xarray(n,1),yarray(n,1), 'oy');
    end
end
hold off
title(['Area set to 0 by strel = ' num2str(radius) ' in mask']);

subplot(1,4,4), imagesc(mask3);
colormap 'gray'
title('Mask used to calculate background');
end