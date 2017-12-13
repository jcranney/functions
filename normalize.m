function [ A_ ] = normalize( A )
%NORMALIZE This is probably already a function but I can't find it
%immediately so I'll just make it myself. Take in a vector or matrix of any
%dimension and scale all of its values linearly to be within 0 and 1
mx = max(A(:));
mn = min(A(:));
A_ = (A-mn)/(mx-mn);
end

