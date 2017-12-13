function F = idft2(m)
% get the dft matrix to be applied to a vector of length m

w = exp(2*pi*1i/m);
a = (0:(m-1))';
b = kron(a,a');
F = w*ones(m);
F = F.^b;
F = kron(F,F)/m^2;