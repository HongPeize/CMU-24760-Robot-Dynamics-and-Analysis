clear; clc; close all;

syms theta_1 theta_2 l_1 l_2 m g real
q=[theta_1; theta_2];

%% 1.1
g_st_0=[eye(3), [l_1+l_2; 0; 0]; 0, 0, 0, 1];
xi_1=[0; 0; 0; 0; 0; 1];
xi_2=[-angvel2skew([0; 0; 1])*[l_1; 0; 0]; 0; 0; 1];
g_st=expm(twist2rbvel(xi_1)*q(1))*expm(twist2rbvel(xi_2)*q(2))*g_st_0;

R_st=g_st(1:3, 1:3);
F_t=[R_st'*[0; -m*g; 0]; 0; 0; 0];
F_t=simplify(F_t);

%% 1.2

J_b=sym(zeros(6, 2));
for i=1:2
    J_b(:, i)=rbvel2twist(g_st\diff(g_st, q(i)));
end

tau_pm=J_b'*(-F_t);
tau_pm=simplify(tau_pm);

%% 1.3
syms m_1 m_2 real

g_s1_0=[eye(3), [l_1; 0; 0]; 0, 0, 0, 1];
g_s1=expm(twist2rbvel(xi_1)*q(1))*g_s1_0;

g_1=[eye(3), [l_1/2; 0; 0]; 0, 0, 0, 1];
g_2=[eye(3), [l_2/2; 0; 0]; 0, 0, 0, 1];

R_1=g_s1(1:3, 1:3);
R_2=g_st(1:3, 1:3);
F_1_s=[R_1'*[0; -m_1*g; 0]; 0; 0; 0];
F_2_s=[R_2'*[0; -m_2*g; 0]; 0; 0; 0];

F_1_b=tform2adjoint(g_1)'*F_1_s;
F_2_b=tform2adjoint(g_2)'*F_2_s;

J_b1=sym(zeros(6, 2));
for i=1:2
    J_b1(:, i)=rbvel2twist(g_s1\diff(g_s1, q(i)));
end

tau_lm=-J_b1'*F_1_b-J_b'*F_2_b;
tau_lm=simplify(tau_lm);

%% 1.4
energy=tau_lm'*tau_lm;
energy=simplify(energy);
