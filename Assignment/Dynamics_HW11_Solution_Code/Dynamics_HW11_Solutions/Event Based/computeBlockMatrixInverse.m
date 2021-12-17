function blockMatrixInv = computeBlockMatrixInverse(x,contactMode)
% Split up state vector into generalized coordinates and velocities
q = x(1:2);
dq = x(3:4);

% Compute A at the current state, only select rows for the current
% contact mode
A = computeA(q);
A = A(contactMode,:);

% Extract the number of constraints
c = size(A,1);

% Compute mass matrix
m = 1;
M = m*eye(2);

% Compute block matrix inverse
blockMatrixInv = inv([M A';A, zeros(c,c)]);