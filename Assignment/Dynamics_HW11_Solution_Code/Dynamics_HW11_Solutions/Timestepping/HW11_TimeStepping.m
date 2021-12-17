clear;clc;close all;
set(0, 'DefaultAxesFontSize', 24);
set(0, 'DefaultLineMarkerSize', 10);
set(0,'defaultfigurecolor',[1 1 1]);
set(0,'DefaultAxesXGrid','off','DefaultAxesYGrid','off')
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
%% Initialization

% Define start and stop times, set a dt to keep outputs at constant time
% interval
t = 0;
tfinal = 3;
h = 0.01;
tspan = t:h:tfinal;

% Initialize state and contact mode
q = [0.2;1];
dq = [0;0];

% Initialize arrays for logging data
xout = [];
lambdaout = [];

% Initialize contact mode
oldContactMode = zeros(0,1);
disp(['Initialize in mode {', num2str(oldContactMode'), '}.']);
%% Main Loop

for i =  1:length(tspan)

    % Solve for new state and contact forces
    [qn, dqn, lambda] = solveEOM(q, dq, h);
    
    % Log state and contact forces
    xout = [xout; qn' dqn'];
    lambdaout = [lambdaout; lambda'];
    
    % Check new contact mode and determine if there has been a transition
    contactMode = find([(-lambda(1))>0;(-lambda(2))>0]);
    if ~isequal(contactMode,oldContactMode)
        transition = {oldContactMode, contactMode};
        disp(['Transition from mode {', num2str(oldContactMode'), '} to mode {', num2str(contactMode'), '} at t = ', num2str(tspan(i+1)), '.']);
    end
    
    % Reset data
    oldContactMode = contactMode;
    q = qn;
    dq = dqn;
end

disp(['Terminate in mode {', num2str(oldContactMode'), '} at t = ', num2str(tfinal), '.']);

% This function shows animation for the simulation, don't modify it
animateHW11(xout, h);

figure
plot(tspan, -lambdaout);
xlabel('t (s)')
ylabel('$U(\lambda)$')
legend('$U(\lambda_1)$','$U(\lambda_2)$',  'location', 'northeast')