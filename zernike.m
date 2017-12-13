function [ z ] = zernike(varargin)
%ZERNIKE z = zernike(mode,imgwidth, discwidth)
%   returns a 2d array matching the zernike function of the required index,
%   assuming a noll indexing convention. image width defaults to 256, and
%   disc width defaults to image width.
%   Detailed explanation goes here
zern = varargin{1};
if nargin == 1
    imgwidth = 256;
    discwidth = imgwidth;
else
    imgwidth = varargin{2};
end
if length(size(imgwidth))>1
    imgwidth = imgwidth(1);
end
if nargin == 2
    discwidth = imgwidth;
elseif nargin == 3
    discwidth = varargin{3};
    if length(size(discwidth))>1
        discwidth = discwidth(1);
    end
end
x = linspace(-1,1,discwidth);
[X,Y] = meshgrid(x,x);
[theta,r] = cart2pol(X,Y);
tmp = zeros(max(zern(:)),1);
tmp(zern) = 1;
z = getZernSum(tmp,r,theta); 
try
    z = padarray(z,[1 1]*(imgwidth-discwidth)/2);
catch
    error('Disc diameter padding error. Choose even square image dimensions.')
end

