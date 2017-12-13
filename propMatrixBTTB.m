function A = propMatrixBTTB(n,a)
%PROPMATRIX(N,A) takes size N and scalars stored in array A to determine
%a propagation matrix.
if length(n)~=1
    error('"n" must be a constant')
end
if length(a)==1
    a_ = zeros((2*n-1)^2,1);
    a_(sub2ind([2*n-1,2*n-1],n,n)) = a;
    a = a_;
elseif length(a)==5 && numel(a)==5
    a_ = zeros((2*n-1)^2,1);
    idx = sub2ind([2*n-1,2*n-1],[n n n n-1 n+1],[n n-1 n+1 n n]);
    a_(idx) = a(:);
    a = a_;
else
    a = a(:);
end
A = BTTB(reshape(a,2*n-1,2*n-1));