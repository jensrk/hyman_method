clear all; close all; clc

%% 4x4 Matrix nicht Hessenberg
A = [ ...
    5.0    0.5      0    0.5; ...
    0.5    3.0      0      0; ...
    0.5      0    1.0    0.5; ...
    0.5      0      0    6.0  ...
    ];
N = 4;

% Ähnlichkeitstransformation in Hessenbergform
A = hess(A);

%% plot des char. Polynoms
% Vektor von x Werten von 0 bis 7
x_vec = linspace(0,7,100)';

% Auswertung mit der Funktion hyman()
[one_over_xN, one_over_xN_prime, c] = hyman(A,x_vec);
p_of_x = c*one_over_xN;

% plot ...
plot(x_vec, p_of_x); hold on
axis([0 7 -10 10])
drawaxis(gca, 'x', 0, 'movelabel', 1)
pause

%% Beginne Newton Verfahren
% finde obere Schranke für größten EW
[~,x_init] = gershgorin_bound(A);

% leerer Vektor für EW
lambda = zeros(4,1);

% initialisiere Variablen für Newton Verfahren
err = 1e-6;

% N-fache Nullstellensuche mit dem Newton Verfahren
for i = 1:N   
    x = [x_init];
    p = [inf];
    while abs(p(end)) > err
        [one_over_xN, one_over_xN_prime, c] = hyman(A,x(end));
        if i == 1
            p = [p; one_over_xN*c];
            x = [x; x(end) - one_over_xN/one_over_xN_prime];
        else
            [q, q_prime] = divider(x(end),lambda(1:i));
            p = [p; one_over_xN*c/q];
            x = [x; x(end) - 1/(one_over_xN_prime/one_over_xN - q_prime/q)];
        end
    end
    % speichere und plotte gefundenen EW
    lambda(i) = x(end);
    plot(lambda(i),0,"rx"); pause
    
    % definiere Nullstellenpolynom zur Deflation von p
    [q, q_prime] = divider(x_vec, lambda(1:i));
    
    % plotte defaltioniertes Polynom
    plot(x_vec, p_of_x./q); pause;
    
    % nächster Startwert x_init > lambda(i)
    x_init = x(end) + 1e-1;
end

function [q, q_prime] = divider(x, eigenvalues)
    N = length(eigenvalues);
    K = length(x);
    A = zeros(K,N);
    for i = 1:N
        A(:,i) = (x - eigenvalues(i));
    end
    q = prod(A, 2);
    q_prime = zeros(size(q));
    for i = 1:N
        I = [1:i-1 i+1:N];
        q_prime = q_prime + prod(A(:,I), 2);
    end
end
