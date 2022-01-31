%A m × n orthogonal matrix Q and a N × N upper-triangular matrix R.
%/* The pseudo-code begins below */
function [Q, R] = simple_QR(A)
[m,n] = size(A);
Im = eye(m);
Qtranspose = Im; %/* Qtranspose is initially the m × m identity matrix */
R = A;
for j = 1:n
    %clear x, v, w; %/* The length of these vectors change as j changes. */
    x = R(j:m,j);
     e1 = [1; zeros(m-j, 1)];
    if x(1) >= 0 
        u = x-norm(x)*e1;
    else
        u = x+norm(x)*e1;
    end
    if norm(u) < 1e-10
        if j < m
            w = [1; zeros(m-j, 1)];
        else
            w = 1;
        end
    else
        w = u/norm(u);
    end
    if j > 1 
        v = [zeros(j-1,1); w];
    else
        v = w;
    end
    P = Im-2*v*v';
    R = P*R;
    Qtranspose = P*Qtranspose;
end
Q = Qtranspose';
end