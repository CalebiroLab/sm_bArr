function [minpix maxpix] = minmaxpix(path, moviefile)

% This function determines the minimum and maximum pixel intensities in an image file.
%
% path: path for the image file
% moviefile: name of the image file,
% minpix: minimum intensity value
% maxpix: maximum inetnsity value
IFO.fpath = path;
cd(path);
info = imfinfo(moviefile, 'tif');
numFrames = size(info,2);
for n = 1 : numFrames
    iC=imread(moviefile,n);
    mini(n)=min(min(iC));
    maxi(n)=max(max(iC));
end
minpix = min(mini);
maxpix = max(maxi);
end
