function a = compute_a(in1)
%COMPUTE_A
%    A = COMPUTE_A(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    19-Nov-2019 19:27:13

x = in1(1,:);
y = in1(2,:);
t2 = y.*2.0;
a = [t2-x;t2+x];