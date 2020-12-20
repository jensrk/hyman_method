function lambda = hyman_newton(A, err)
%HYMAN_NEWTON Summary of this function goes here
%   Detailed explanation goes here

B = hess(A);
N = size(A,1);

% finde obere Schranke für größten EW
[lower, x_init] = gershgorin_bound(B);

% leerer Vektor für EW
lambda = zeros(N,1);

% N-fache Nullstellensuche mit dem Newton Verfahren
for i = 1:N   
    x = [x_init];
    p = [inf];
    while abs(p(end)) > err
        [one_over_xN, one_over_xN_prime, c] = hyman(B,x(end));
        if i == 1 || i == N
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
    if i < N-1
        x_init = x(end) + 1e-1;
    else
        x_init = lower;
    end
    
    fprintf( ...
        "Found eigenvalue: %.4f in %i iterations\n", ...
        lambda(i), length(x)-1)
end

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


