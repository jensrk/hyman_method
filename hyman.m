function [one_over_xN, one_over_xN_prime, c] = hyman(B,mu)
%HYMAN evaluates det(B - x*I) with the Hyman Method
%   B has to be of Hessenberg type
%   x can be scalar or array of scalar

N = size(B,1);
%   For (B - mu*I)*x = e_1 we have
%   p(mu) = det(B - mu*I) = c/xN with xN = x(N)
%   c = (-1)^(N-1)*b_21*b_32*...*b_N,N-1
%   p'(mu) = c*(1/xN)'

% product of subdiagonal
c = (-1)^(N-1)*prod(diag(B(2:N,1:N-1)));
if c == 0
    error('Subdiagonal of Hessenbergmatrix has zero')
end

one_over_xN = zeros(size(mu));
one_over_xN_prime = zeros(size(mu));
% for every entry of mu evaluate p(mu) with Hyman
for i = 1:length(mu)
    % M := (B - mu*I)
    M = (B - mu(i)*eye(N));
    
    % solve Rows 2 to N
    v = M(2:N,1:N-1)\(-M(2:N,N));
    % 1/xN is computed from first row
    one_over_xN(i) = M(1,:)*[v; 1];
    
    % solve Rows 2 to N
    z = M(2:N,1:N-1)\[v(2:end); 1];
    % (1/xN)' is computed from first row
    one_over_xN_prime(i) = (M(1,1:N-1)*z - v(1));
end
end

