function z = getAp(pix,varargin)
if length(size(pix))>1
    pix = pix(1);
end
if nargin>1
    D = varargin{1};
else
    D = pix;
end
if length(size(D))>1
    D = D(1);
end
z = zeros(pix);
x = linspace(-(pix-1)/2,(pix-1)/2,pix);
[X,Y] = meshgrid(x,x);
[~,r] = cart2pol(X,Y);
z(r<=D/2)=1;