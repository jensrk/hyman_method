clear all; close all; clc

%% Beispiel Matrix in Hessenbergform
A = [-1     1     0;
      2     1    -2;
      0     2    -2;];
N = 3;

%% Auswertung des charakteristischen Polynoms
% über Werte zwischen -3 und 2
x = linspace(-3,2,40);

% mit der Funktion hyman()
[one_over_xN, one_over_xN_prime, c] = hyman(A,x);

% plot...
plot(x, one_over_xN*c); hold on
drawaxis(gca, 'x', 0, 'movelabel', 1)

%% Suche den größten Eigenwert mit dem Newton-Verfahren
% finde obere Schranke für größten EW
[~,x_init] = gershgorin_bound(A)
p = inf;
x = x_init;
store_x = [x_init];
while abs(p) > 1e-6
    [one_over_xN, one_over_xN_prime, c] = hyman(A,x);
    p = one_over_xN*c;
    x = x - one_over_xN/one_over_xN_prime;
    store_x = [store_x; x];
end
lambda_max = x
plot(x,0,"rx")

%% Überprüfe Konvergenzordnung
figure(2)
eigenvalues = eig(A);
diffs = abs(store_x - eigenvalues(1));
loglog(diffs(1:end-1),diffs(2:end),"--x"); hold on
x = logspace(-10,2,50);
loglog(x, x.^1)
loglog(x, x.^2)
loglog(x, x.^3)
legend(["Newton Iteration","O(|x-x*|)","O(|x-x*|^2)","O(|x-x*|^3)"])
xlabel("|x_n - x*|")
ylabel("|x_{n+1} - x*|")