function [lower, upper] = gershgorin_bound(A)
N = size(A,1);
x_intersection = zeros(2*N,1);
for i = 1:N
    I = [1:(i-1) (i+1):N];
    radius = sum(abs(A(i,I)));
    x_intersection(2*i-1) = A(i,i) - radius;
    x_intersection(2*i)   = A(i,i) + radius;
end
lower = min(x_intersection);
upper = max(x_intersection);
end

