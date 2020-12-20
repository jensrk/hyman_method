clear all; clc

%% Tridiagonal Toeplitz Matrix
N = 20;
e = ones(N,1);
A = full(spdiags([-e,2*e,-e],[-1,0,1],N,N));

% Eigenwerte sind bekannt in geschlossener Form:
% k = 1:n;
% Lambda = 2 - 2*cos(k*pi/(n+1));
eigenvalues_exact = 2 - 2*cos(pi*(N:-1:1)'/(N+1));

%% evaluiere hyman_newton() Funktion
tic
eigenvalues = hyman_newton(A, 1e-10);
toc
disp("Maximum error to exact values is")
disp(max(abs(eigenvalues - eigenvalues_exact)))
