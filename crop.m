function Z = crop(X,cropto)
% takes the centre cropto(1) by cropto(2) rectangle of points from X

[m,n] = size(X);
Z = X(fix(m/2-cropto(1)/2+1:m/2+cropto(1)/2),fix(n/2-cropto(2)/2+1:n/2+cropto(2)/2));