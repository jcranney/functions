function z = getZernSum(x,r,theta)
nollmax = length(x);
global marray narray
if isempty(narray)
    [marray, narray] = genNollindices;
end
[a,b] = size(r);
z = zeros(nollmax*a,b);
nolls = 1:nollmax;
idx = r<=1;
idxbig = repmat(idx,1,nollmax);
eyescaled = kron(x',eye(a));
try
    z(idxbig) = zernfun(narray((1:length(nolls))),marray((1:length(nolls))),r(idx),theta(idx));
    z = reshape(z,b,nollmax*a);
    z = z*eyescaled';
catch
    error('Need to include genNollIndicies.m in initialisation to get marray and narray globals')
end
