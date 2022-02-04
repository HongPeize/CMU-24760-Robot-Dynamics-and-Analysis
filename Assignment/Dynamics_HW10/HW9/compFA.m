function contactMode = compFA(x)
% Split up state vector into generalized coordinates and velocities
q = x(1:2);
dq = x(3:4);

% Define J, set of contact modes
J = {[], [1], [2], [3], [1;2], [1;3]};

% Define constraint functions
a = compute_a(q);

% Identify possible contact modes where constraint function is 0
possibleModes = [ []; find(abs(a)<1e-6)];

% Test each contact mode
for k = 1:length(J)
    % Define current contact mode
    mode = J{k};
    notMode = setdiff(possibleModes, mode);
    
    % Check that current mode is included in possibleModes, otherwise pass
    if all(ismember(mode, possibleModes))
        
        % Compute reset map into possible modes, reformat the forces
        [ddq,lambda] = solveEOM(x,mode);
        [ddq_union, lambda_union] = solveEOM(x,possibleModes);
        lambda_union(possibleModes) = lambda_union;
        
        % Generate complementarity conditions
        cond_1 = all(-lambda>=0);
        cond_2 = all(-lambda_union(notMode)<=0);
        
        % If both complementarity conditions are satisfied, switch modes
        if cond_1 && cond_2
            contactMode = mode;
            return;
        end
    end
end

end