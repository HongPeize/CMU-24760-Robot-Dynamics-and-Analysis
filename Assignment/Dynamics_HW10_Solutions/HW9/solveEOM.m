function [ddq,lambda] = solveEOM(x,contactMode)
% Split up state vector into generalized coordinates and velocities
q = x(1:2);
dq = x(3:4);

% Compute A and dA at the current state, only select rows for the current
% contact mode
A = computeA(q);
dA = computedA(q,dq);
A = A(contactMode,:);
dA = dA(contactMode,:);

% Extract the number of constraints
c = size(A,1);

% Compute EOM matrices
N = [0;9.8];
C = [0, 0;0, 0];
Y = [0;0];

% Compute block matrix inverse
blockMatrixInv = computeBlockMatrixInverse(x,contactMode);

% Solve EOM
sol = blockMatrixInv*[Y - N;zeros(c,1)] - blockMatrixInv*[C;dA]*dq;

% Extract accelerations and contact forces
ddq = sol(1:2);
if length(sol)>=3
    lambda = sol(3:end);
else
    lambda = [];
end
end