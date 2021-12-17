clear;clc;close all;
%% Initialization (Do Not Change This Section)

% Define start and stop times, set a dt to keep outputs at constant time
% interval
tstart = 0;
tfinal = 3;
dt = 0.01;

% Initialize state and contact mode
q0 = [0.2;1];
dq0 = [0;0];
x0 = [q0;dq0];
contactMode = [];
%% Main Loop

% Tell ode45 what event function to use and set max step size to make sure
% we don't miss a zero crossing
options = odeset('Events', @guardFunctions,'MaxStep',0.01);

% Initialize output arrays
tout = tstart;
xout = x0.'; % collection of state vectors at each time
teout = []; % collection of ode45 output
xeout = []; % collection of ode45 output
ieout = []; % collection of ode45 output
while tstart < tfinal
    % Initialize simulation time vector
    tspan = [tstart:dt:tfinal];
    
    % Simulate
    [t,x,te,xe,ie] = ode45(@dynamics,tspan,x0,options,contactMode);
    
    % Sometimes the events function will record a nonterminal event if the
    % initial condition is a zero. We want to ignore this, so we will only
    % use the last row in the terminal state, time, and index.
    if ~isempty(ie)
        te = te(end,:);
        xe = xe(end,:);
        ie = ie(end,:);
    end
    
    % Log output
    nt = length(t);
    tout = [tout; t(2:nt)];
    xout = [xout; x(2:nt,:)];
    teout = [teout; te];
    xeout = [xeout; xe];
    ieout = [ieout; ie];
    
    % Quit if simulation completes
    if isempty(ie) 
        disp('Final time reached');
        break; % abort if simulation has completed
    end
    
    % If flag was caused by a_i < 0 (i not in contact mode), compute the
    % proper contact mode via IV complemetarity
    if any(ie <= length(setdiff([1:2], contactMode)))
        % Determine next contact mode
        contactMode = compIV(xe');
        
        % Compute the reset map
        [dq_p, p_hat] = computeResetMap(xe',contactMode);
        xe = [xe(1:2),dq_p'];
        
        % Report contact mode transition
        disp(['Transition to contact mode {', num2str(contactMode'), '} at time t = ', num2str(te), ' s.'])
    end
    
    % Check to see if there should be liftoff (positive lambda), if so
    % compute the proper contact mode via FA complementarity
    [ddq,lambda] = solveEOM(xe',contactMode);
    if (-lambda)<=0
        contactMode = compFA(xe');
        
        % Report contact mode transition
        disp(['Transition to contact mode {', num2str(contactMode'), '} at time t = ', num2str(te), ' s.'])
    end
    
    % Update initial conditions for next iteration
    x0 = xe';
    tstart = t(end);
    
    % Stop if the particle comes to a rest
    if all(abs(dq_p)<1e-6)
        break;
    end
end

% This function shows animation for the simulation, don't modify it
animateHW11(xout, dt);

