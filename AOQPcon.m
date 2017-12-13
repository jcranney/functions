function [h,J] = AOQPcon(cube,Omega)
% AOQP(OMEGA,Y) takes the propogated samples in OMEGA and the actual
% outputs in Y, and determines the best set of H to minimize the difference
% using a Quadratic Program (QP).
% N is the length of H
% M is the length of Y
thresh = cube.thresh;
y = cube.getY;
options = optimoptions('fmincon');
options.Algorithm = 'sqp';
N = size(Omega,2);
M = size(Omega,1);
if length(y)~=M
    error('"Omega" and "y" must be the same height');
end
% We require that the rows of A be positive, which is the following
% constraint:
tmp = [ones(1,cube.n) zeros(1,cube.n-1)];
w = [];
for ni = 1:cube.n
    w = [w;tmp];
    tmp = circshift(tmp,[0 1]);
end
W = kron(w,w);
W(:,cube.getBad) = [];
% We require that all h be positive, and that the sum of h be less than 1.
% These restrictions can be loosened but this is a decent start.
% Let:
A = [-eye(N) ; W];
b = [zeros(N,1) ; ones(cube.n^2,1)];
alpha = 1e12;
[h,J] = quadprog(alpha*2*(Omega'*Omega),alpha*(-2*Omega'*y),A,b);
%[h,J] = fmincon(@(h) 1e10*(Omega*h-y)'*(Omega*h-y),ones(N,1),A,b,[],[],[],[],[],options);
h(h<thresh) = 0;