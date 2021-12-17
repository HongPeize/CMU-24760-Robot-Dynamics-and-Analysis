function contactMode = compIV(x)
% Split up state vector into generalized coordinates and velocities
q = x(1:2);
dq = x(3:4);

% Define constraint functions
a = compute_a(q);

% Identify possible contact modes where constraint function is 0
possibleModes = find(abs(a)<1e-6);

% Define J, set of contact modes
J = {[1], [2], [1;2]};

% Test each contact mode
for k = 1:length(J)
    % Define current contact mode
    mode = J{k};
    notMode = setdiff(possibleModes, mode);
    
    % Check that current mode is included in possibleModes, otherwise pass
    if all(ismember(mode, possibleModes))
        % Compute reset map into mode
        [dq_p, p_hat] = computeResetMap(x,mode);
        
        % Compute reset map into possible modes, reformat the impulses
        [dq_p_union, p_hat_union] = computeResetMap(x,possibleModes);        
        p_hat_union(possibleModes) = p_hat_union;
        
        % Generate complementarity conditions
        cond_1 = all(-p_hat>=0);
        cond_2 = all(-p_hat_union(notMode)<0);
        
        % If both complementarity conditions are satisfied, switch modes
        if cond_1 && cond_2
            contactMode = mode;
        end
    end
end
end