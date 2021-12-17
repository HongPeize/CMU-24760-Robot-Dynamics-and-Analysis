function dx = dynamics(t,x,contactMode)
% Initialize time derivative of state vector as column vector
dx = zeros(length(x),1);

% Split up state vector into generalized coordinates and velocities
q = x(1:2);
dq = x(3:4);

% Solve the equations of motion
[ddq, lambda] = solveEOM(x,contactMode);

% Time derivative of generalized positions are the generalized velocities
dx(1:2) = dq;

% Time derivative of generalized velocities are accelerations (from EOM)
dx(3:4) = ddq;

end