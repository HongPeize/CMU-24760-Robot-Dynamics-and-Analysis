function A = computeA(in1)
%COMPUTEA
%    A = COMPUTEA(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    04-Nov-2019 15:37:45

x = in1(1,:);
y = in1(2,:);
A = reshape([0.0,1.0,x.*2.0-4.0,1.0,1.0,y.*2.0-2.0],[3,2]);
