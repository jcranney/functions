function [marray, narray] = genNollindices
nmax = 20;
narray = zeros(nchoosek(nmax+1,2),1);
marray = narray;
idx1 = 1;
for a = 1:nmax
    idx2 = nchoosek((a+1),2);
    narray(idx1:idx2) = a-1;
    msub = round(linspace(-(a-1),(a-1),a));
    msub_ = zeros(a,1);
    for b = 1:(a)
        msub_temp = (-1)^(b+mod(idx1,2)+1)*round(min(abs(msub+(mod(idx1,2)==0)*0.001)));
        msub_(b) = msub_temp(1);
        msub(msub==msub_(b))=[];
    end
    marray(idx1:idx2) = msub_;
    idx1 = idx2+1;
end