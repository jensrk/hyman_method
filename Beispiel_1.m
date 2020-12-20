clear all; close all; clc

%% Beispiel Matrix in Hessenbergform
A = [-1     1     0;
      2     1    -2;
      0     2    -2;];
N = 3;

%% Auswertung des charakteristischen Polynoms
% über Werte zwischen -3 und 2
x = linspace(-3,2,40);

% mit der MATLAB Funktion det()
p = zeros(size(x));
for i = 1:length(p)
    p(i) = det(A - x(i)*eye(N));
end

% mit der Funktion hyman()
[one_over_xN, one_over_xN_prime, c] = hyman(A,x);

%% plot...
plot(x, p, "-.*"); hold on
plot(x, one_over_xN*c, "--o")
plot(x, one_over_xN_prime*c)
for lambda = eig(A)
    plot(real(lambda),imag(lambda),'rx')
end
drawaxis(gca, 'x', 0, 'movelabel', 1)
legend(["p(x)","p(x) Hyman","p'(x) Hyman","Eigenwerte"])
