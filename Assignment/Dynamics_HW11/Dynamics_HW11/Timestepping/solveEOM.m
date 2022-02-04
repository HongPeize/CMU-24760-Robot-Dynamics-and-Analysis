function [qn, dqn, lambda] = solveEOM(q, dq, h)

% Declare symbolic variables
syms xn yn dxn dyn real
syms lambda_1 lambda_2 real

% Group into vectors
qn = [xn;yn];
dqn = [dxn;dyn];
lambda = [lambda_1; lambda_2];

% Construct M and A matrices
M = eye(2);
a = compute_a(qn);
A = computeA(qn);

% Equalities
eq1 = M*(dqn - dq) - (-A'*lambda + h*[0;-9.8]);
eq2 = (qn - q) - h*dqn;
eq3 = a.*lambda;

% Inequalities 
in1 = a>=0;
in2 = lambda<=0;

% Group equations, inequalities, and variables
eqns = [eq1;eq2;eq3;in1;in2];
vars = [qn;dqn;lambda];

% Solve and save into symbolic struct
sol = solve(eqns, vars);

% Extract doubles from solve output
qn = double([sol.xn(1);sol.yn(1)]);
dqn = double([sol.dxn(1);sol.dyn(1)]);
lambda = double([sol.lambda_1;sol.lambda_2]);