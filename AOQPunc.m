function [h,J] = AOQPunc(cube,Omega)
% AOQP(OMEGA,Y) takes the propogated samples in OMEGA and the actual
% outputs in Y, and determines the best set of H to minimize the difference
% using a Quadratic Program (QP).
% N is the length of H
% M is the length of Y
y = cube.getY;
options = optimoptions('fmincon');
options.Algorithm = 'sqp';
N = size(Omega,2);
M = size(Omega,1);
if length(y)~=M
    error('"Omega" and "y" must be the same height');
end
% We require that all h be positive, and that the sum of h be less than 1.
% These restrictions can be loosened but this is a decent start.
% Let:


alpha = 1e5;
[h,J] = quadprog(alpha*2*(Omega'*Omega),alpha*(-2*Omega'*y));