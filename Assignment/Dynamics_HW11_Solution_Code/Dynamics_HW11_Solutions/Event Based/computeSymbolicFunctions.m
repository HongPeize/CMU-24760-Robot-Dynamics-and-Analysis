syms x y dx dy real
q = [x;y];
dq = [dx;dy];

m = 1;
g = 9.8;

a1 = 2*y - x;
a2 = 2*y + x;

a = [a1;a2];
A = jacobian(a,q);
dA = sym(zeros(size(A)));
for i = 1:length(q)
    dA = dA + diff(A, q(i))*dq(i);
end

matlabFunction(a,'File', 'compute_a', 'Vars', {q});
matlabFunction(A,'File', 'computeA', 'Vars', {q});
matlabFunction(dA,'File', 'computedA', 'Vars', {q,dq});