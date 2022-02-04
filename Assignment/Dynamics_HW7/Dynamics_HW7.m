%% Initialization

clear;clc;close all

% States
syms theta1 theta2 o_x o_y o_theta dtheta1 dtheta2 do_x do_y do_theta ddtheta1 ddtheta2 ddo_x ddo_y ddo_theta real

% Physics
syms m_l L I_l w m_o I_o g real

% Inputs
syms tau_1 tau_2 real

% Constraint forces
syms lambda_1 lambda_2
lambda = [lambda_1; lambda_2];

% Finger coordinates
theta = [theta1;theta2];
dtheta = [dtheta1;dtheta2];
ddtheta = [ddtheta1;ddtheta2];

% Body coordinates
x = [o_x; o_y; o_theta];
dx = [do_x; do_y; do_theta];
ddx = [ddo_x; ddo_y; ddo_theta];

% Local coordinates
q = [theta;x];
dq = [dtheta; dx];
ddq = [ddtheta; ddx];

%% Problem 1.1

% Compute generalized inertia matrices for each link and object
M1 = [m_l 0 0;
    0 m_l 0;
    0 0 I_l];
M2 = M1;
M_o =  [m_o 0 0;
    0 m_o 0;
    0 0 I_o];

% Compute body Jacobians for each link
Jb_sL1 = [0 0;
    L/2 0
    1 0];
Jb_sL2 = [L*sin(theta2) 0
    L*cos(theta2)+L/2 L/2
    1 1];

% Compute manipulator inertia tensor
Mhat = [Jb_sL1'*M1*Jb_sL1 + Jb_sL2'*M2*Jb_sL2,  zeros(length(dtheta),length(dx))
    zeros(length(dtheta),length(dx))',    M_o];

Jb_po = [cos(o_theta) sin(o_theta) 0
    -sin(o_theta) cos(o_theta) 0
    0 0 1];

% Use body jacobians to put in terms of local coordinates rather than body
% twist (go from Mhat to Mbar)
Mbar = [eye(length(dtheta)), zeros(length(dtheta),length(dx))
    zeros(length(dx),length(dtheta)), Jb_po']*Mhat*[eye(length(dtheta)), zeros(length(dtheta),length(dx))
    zeros(length(dx),length(dtheta)), Jb_po];

% Compute Coriolis matrix (Uncomment below to evaluate, takes some time)
Cbar  = sym(zeros(length(q),length(q)));
for ii = 1:length(q)
    for jj = 1:length(q)
        for kk = 1:length(q)
            Cbar(ii,jj) = Cbar(ii,jj) + 1/2*(diff(Mbar(ii,jj),q(kk)) + diff(Mbar(ii,kk),q(jj)) - diff(Mbar(jj,kk),q(ii)))*dq(kk);
        end
    end
end

% Compute nonlinear and applied force terms
% If we assume g is negative, we should compute a negative V
V = m_l*g*(w/2+L/2*sin(theta1)) + m_l*g*(w/2 + L*sin(theta1) + L/2*sin(theta1+theta2)) + m_o*g*o_y;
Nbar = jacobian(V, q)';
Y = [tau_1;tau_2;0;0;0];

% Contact wrench basis
Bc = eye(6,2);

% Compute transform from O to C
g_oc = [1 0 0 0;
    0 1 0 -w/2;
    0 0 1 0;
    0 0 0 1];

% Compute Grasp Map
G = [tform2adjoint(inv(g_oc))'*Bc];

% Switch to planar
G = G([1,2,6],:);

% Compute body jacobian of object wrt palm
Jb_po = [cos(o_theta) sin(o_theta) 0
    -sin(o_theta) cos(o_theta) 0
    0 0 1];

% Compute grasp map in local coordinates
GbarT = G'*Jb_po;

% Compute hand Jacobian
g_sc = [cos(o_theta) -sin(o_theta) 0 L*cos(theta1) + L*cos(theta1+theta2) % o_theta = -o_theta from previous since o_theta
    sin(o_theta) cos(o_theta) 0 L*sin(theta1) + L*sin(theta1+theta2)      % rotates the s frame, not the c frame
    0 0 1 0
    0 0 0 1];

Js_sf = [0 sin(theta1)*L
    0 -cos(theta1)*L
    0 0
    0 0
    0 0
    1 1];

Jh = Bc'*inv(tform2adjoint(g_sc))*Js_sf;

% Compute A matrix
A = [-Jh GbarT];

% Compute dA matrix
dA = sym(zeros(size(A)));
for i = 1:length(A)
    dA = dA + diff(A, q(i))*dq(i);  % Chain rule
end
c = size(A,1); % number of constraints

% Compute equations of motion
EOM = Mbar*ddq + Cbar*dq + Nbar + A'*lambda - Y;

%% Problem 1.2

% Calculate block matrix inverse
blockMinv = inv([Mbar A';A zeros(c,c)]);

% Rename blocks of inverse
M_dag = blockMinv(1:length(q), 1:length(q));
A_dag = blockMinv(length(q)+1:end, 1:length(q));
Lambda = blockMinv(length(q)+1:end, length(q)+1:end);

% Compute acceleration
ddq_massive = M_dag*(Y - Cbar*dq - Nbar) - A_dag'*dA*dq;

% Numerical evaluation
symbols = [tau_1; tau_2; L; m_l; I_l; w; m_o; I_o; g; q; dtheta1; dtheta2; do_x; do_y; do_theta];
values = [200; 0; 0.1; 1; 8.33e-4; 0.2; 24; 0.16; 9.81; pi/2; -pi/2; 0.1; 0.2; 0; 0; 0; 0; 0; 0];

ddq_eval_massive = subs(ddq_massive, symbols, values);

%% Problem 1.3

% Compute massless EOM
EOM_massless = subs(EOM, [m_l, I_l], [0,0]);

%% Problem 1.4

% Numerical evaluation
symbols = [tau_1; tau_2; L; m_l; I_l; w; m_o; I_o; g; q; dtheta1; dtheta2; do_x; do_y; do_theta];
values = [200; 0; 0.1; 0; 0; 0.2; 24; 0.16; 9.81; pi/2; -pi/2; 0.1; 0.3; 0; 0; 0; 0; 0; 0];

% In the massless case, Mbar is not invertable so you have to do the block matrix inverse approach
ddq_eval_massless = subs(ddq_massive, symbols, values);

% Compute error percentage
error_from_massless = abs(1 - ddq_eval_massless./ddq_eval_massive)*100;

%% Problem 1.5

% Compute constraints
lambda_frictionless = lambda_2;
c_frictionless = length(lambda_frictionless); % number of constraints, should be same as length of lambda

% Compute contact velocity map A and its time derivative dA
A_frictionless = A(2, :);
dA_frictionless = sym(zeros(size(A_frictionless)));
for i = 1:length(A_frictionless)
    dA_frictionless = dA_frictionless + diff(A_frictionless, q(i))*dq(i);  % Chain rule
end

% Compute equations of motion and check against previously computed
EOM_frictionless = Mbar*ddq + Cbar*dq + Nbar + A_frictionless'*lambda_frictionless - Y;

%% Problem 1.6

% Compute massless and frictionless EOM
EOM_massless_frictionless = subs(EOM_frictionless, [m_l, I_l], [0,0]);

% Compute rank
rank_block = rank(subs([Mbar A_frictionless';A_frictionless zeros(c_frictionless,c_frictionless)], [m_l, I_l], [0,0]));
