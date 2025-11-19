function [best_pos, best_val] = SimplePSO(obj_func, n_vars, lb, ub, options)
    % SimplePSO - Basic Particle Swarm Optimization implementation
    %
    % Inputs:
    %   obj_func: Function handle for objective function
    %   n_vars: Number of variables
    %   lb: Lower bounds (1 x n_vars)
    %   ub: Upper bounds (1 x n_vars)
    %   options: Struct with fields 'MaxIter', 'SwarmSize', 'Display' (optional)
    %
    % Outputs:
    %   best_pos: Best position found
    %   best_val: Best objective value found

    % Default options
    if nargin < 5; options = struct(); end
    if ~isfield(options, 'MaxIter'); options.MaxIter = 100; end
    if ~isfield(options, 'SwarmSize'); options.SwarmSize = 50; end
    if ~isfield(options, 'Display'); options.Display = 'iter'; end
    
    w = 0.7; % Inertia weight
    c1 = 1.5; % Cognitive weight
    c2 = 1.5; % Social weight
    
    n_particles = options.SwarmSize;
    
    % Initialize particles
    positions = repmat(lb, n_particles, 1) + rand(n_particles, n_vars) .* repmat(ub - lb, n_particles, 1);
    velocities = zeros(n_particles, n_vars);
    
    personal_best_pos = positions;
    personal_best_val = inf(n_particles, 1);
    
    global_best_pos = positions(1, :);
    global_best_val = inf;
    
    % Initial evaluation
    for i = 1:n_particles
        val = obj_func(positions(i, :));
        personal_best_val(i) = val;
        if val < global_best_val
            global_best_val = val;
            global_best_pos = positions(i, :);
        end
    end
    
    % Main loop
    for iter = 1:options.MaxIter
        % Update velocities and positions
        r1 = rand(n_particles, n_vars);
        r2 = rand(n_particles, n_vars);
        
        velocities = w * velocities + ...
                     c1 * r1 .* (personal_best_pos - positions) + ...
                     c2 * r2 .* (repmat(global_best_pos, n_particles, 1) - positions);
                 
        positions = positions + velocities;
        
        % Boundary handling (clamp)
        positions = max(positions, repmat(lb, n_particles, 1));
        positions = min(positions, repmat(ub, n_particles, 1));
        
        % Evaluation
        for i = 1:n_particles
            val = obj_func(positions(i, :));
            
            if val < personal_best_val(i)
                personal_best_pos(i, :) = positions(i, :);
                personal_best_val(i) = val;
                
                if val < global_best_val
                    global_best_val = val;
                    global_best_pos = positions(i, :);
                end
            end
        end
        
        if strcmp(options.Display, 'iter')
            fprintf('Iteration %d: Best Value = %.6f\n', iter, global_best_val);
        end
    end
    
    best_pos = global_best_pos;
    best_val = global_best_val;
end

