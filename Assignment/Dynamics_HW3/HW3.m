clear all; clc; close all;

%% Problem 2.1
gst0 = [1, 0, 0, 1407;
    0, 1, 0, 0;
    0, 0, 1, 1855;
    0, 0, 0, 1];

%% Problem 2.2

% Define symbolic variables
syms q1 q2 q3 q4 q5 q6 real;
q = [q1 q2 q3 q4 q5 q6]';

% Define known frame offsets
lx2 = 320; lx5 = 887; lx6 = 200; lz2 = 680; lz3 = 975; lz4 =200;

% Define rigid body transformations between successive links

% All frames are same as the stationary base frame
gs1 = [cos(q1) -sin(q1) 0 0;sin(q1) cos(q1) 0 0;0 0 1 0;0 0 0 1];
g12 = [cos(q2) 0 sin(q2) lx2;0 1 0 0;-sin(q2) 0 cos(q2) lz2;0 0 0 1];
g23 = [cos(q3) 0 sin(q3) 0;0 1 0 0;-sin(q3) 0 cos(q3) lz3;0 0 0 1];
g34 = [1 0 0 0;0 cos(q4) -sin(q4) 0;0 sin(q4) cos(q4) lz4;0 0 0 1];
g45 = [cos(q5) 0 sin(q5) lx5;0 1 0 0;-sin(q5) 0 cos(q5) 0;0 0 0 1];
g5t = [1 0 0 lx6;0 cos(q6) -sin(q6) 0;0 sin(q6) cos(q6) 0;0 0 0 1];

% Another frame configuration
% gs1 = [cos(q1) -sin(q1) 0 0;sin(q1) cos(q1) 0 0;0 0 1 0;0 0 0 1];
% g12 = [sin(q2) cos(q2) 0 lx2;0 0 1 0;cos(q2) -sin(q2) 0 lz2;0 0 0 1];
% g23 = [cos(q3) -sin(q3) 0 lz3;sin(q3) cos(q3) 0 0;0 0 1 0;0 0 0 1];
% g34 = [cos(q4) -sin(q4) 0 lz4;0 0 1 0;-sin(q4) -cos(q4) 0 0;0 0 0 1];
% g45 = [cos(q5) -sin(q5) 0 0;0 0  -1 0;sin(q5) cos(q5) 0 lx5;0 0 0 1];
% g5t = [0 sin(q6) cos(q6) 0;1 0 0 lx6;0 cos(q6) -sin(q6) 0;0 0 0 1];

% Compute forward  kinematics
gst = simplify(gs1*g12*g23*g34*g45*g5t);

% Compute gst(0) with assigning all symbolic variables as 0
gst0_sym = double(subs(gst, q, zeros(6, 1)));

% Compare gst0 to gst0_sym
disp('Difference between gst0 and gst0_sym: ')
disp(gst0 - gst0_sym)

%% Problem 2.3

% Define joint twists in initial configuration
xi1 = [0;0;0;0;0;1];
xi2 = [-lz2;0;lx2;0;1;0];
xi3 = [-(lz2+lz3);0;lx2;0;1;0];
xi4 = [0;lz2+lz3+lz4;0;1;0;0];
xi5 = [-(lz2+lz3+lz4);0;lx2+lx5;0;1;0];
xi6 = [0;lz2+lz3+lz4;0;1;0;0];

% Compute product of exponentials
gst_exp = expm(twist2rbvel(xi1)*q1)*expm(twist2rbvel(xi2)*q2)*expm(twist2rbvel(xi3)*q3)*expm(twist2rbvel(xi4)*q4)...
    *expm(twist2rbvel(xi5)*q5)*expm(twist2rbvel(xi6)*q6)*gst0;
gst_exp = simplify(gst_exp);

% Compare to 2.2
disp('Difference between gst and gst_exp:')
disp(simplify(gst_exp - gst))

%% Problem 2.4

% Define goal position as initial configuration +100mm in +y direction 
gDes = gst0;
gDes(2, 4) = gDes(2, 4) + 100;

% Define an optimization problem
F = sum(sum((gst-gDes).^2));
fun = matlabFunction(F, 'var', {q}); 

% Call fminunc to solve Inverse Kinematics
options = optimoptions(@fminunc, 'Display', 'iter');
[q_sol, fval] = fminunc(fun, zeros(6, 1), options);
gAchieved = double(subs(gst, q, q_sol));

% Compare gAchieved with gDes
disp('Difference between gst and gst_exp:')
disp(gAchieved - gDes)
disp('Optimization tolerance:')
disp(fval)

%% Problem 2.5

% Calculate the inverse of Foward Kinematics transformation
gst_inv = inv(gst);

% Calculate each portion of the spatial Jacobian and get Js
J1 = rbvel2twist(diff(gst , q1)*gst_inv);
J2 = rbvel2twist(diff(gst , q2)*gst_inv);
J3 = rbvel2twist(diff(gst , q3)*gst_inv);
J4 = rbvel2twist(diff(gst , q4)*gst_inv);
J5 = rbvel2twist(diff(gst , q5)*gst_inv);
J6 = rbvel2twist(diff(gst , q6)*gst_inv);

Js = simplify([J1,J2,J3,J4,J5,J6]);

%% Problem 2.6

% Calculate adjoint of product of exponential map
xi2_prime = tform2adjoint(expm(twist2rbvel(xi1)*q1))*xi2;
xi3_prime = tform2adjoint(expm(twist2rbvel(xi1)*q1)*expm(twist2rbvel(xi2)*q2))*xi3;
xi4_prime = tform2adjoint(expm(twist2rbvel(xi1)*q1)*expm(twist2rbvel(xi2)*q2)*expm(twist2rbvel(xi3)*q3))*xi4;
xi5_prime = tform2adjoint(expm(twist2rbvel(xi1)*q1)*expm(twist2rbvel(xi2)*q2)*expm(twist2rbvel(xi3)*q3)...
    *expm(twist2rbvel(xi4)*q4))*xi5;
xi6_prime = tform2adjoint(expm(twist2rbvel(xi1)*q1)*expm(twist2rbvel(xi2)*q2)*expm(twist2rbvel(xi3)*q3)...
    *expm(twist2rbvel(xi4)*q4)*expm(twist2rbvel(xi5)*q5))*xi6;

% Calculate and compare spatial Jacobian
Js_exp = [xi1 xi2_prime xi3_prime xi4_prime xi5_prime xi6_prime];

% Compare to 2.6
disp('Difference between gst and gst_exp:')
disp(simplify(Js - Js_exp))

%% Problem 2.7

% Define body twist in initial configuration
Vb = [0;1;0;0;0;0];

% Convert to spatial twist through adjoint, at initial configuration
Vs = tform2adjoint(gst)*Vb;
Vs = double(subs(Vs, q, zeros(6, 1)));

% Compute rank of spatial jacobian (singular if rank < 6)
Js0 = subs(Js, q, zeros(6, 1));
rank_Js0 = rank(Js0);
disp('Rank of Js0 is:')
disp(rank_Js0)

% We cannot achieve the velocity as spatial Jacobian is singular.
