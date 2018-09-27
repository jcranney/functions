function Z = sqrvec(x)
% takes a column vector and reshapes it to be square, just because im sick of
% going to vectors so easily and then doing this nonsense to go back. Only
% works for square matrices, not general oblongs dimensions.

l = size(x,1);
if l <= 1
    error('x must be a column vector');
end
Z = reshape(x,sqrt(l),sqrt(l));
