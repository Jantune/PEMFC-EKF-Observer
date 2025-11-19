function identified_params = Identify_Thermal_Parameters()
    % Identify_Thermal_Parameters
    % Identifies thermal model parameters using Global PSO on sparse snapshots.
    
    current_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile(current_dir, '../Data'));
    addpath(fullfile(current_dir, '../Models'));
    
    data = LoadData();
    stack_config = LoadStackConfig();
    config = stack_config;
    
    if isempty(data.time)
        fprintf('No data available. Returning dummy thermal parameters.\n');
        identified_params = config.defaults;
        return;
    end
    
    % Determine N
    if isfield(data, 'V_cell') && ~isempty(data.V_cell)
        config.N = size(data.V_cell, 2);
    else
        config.N = 180; % Example
    end
    
    % Determine dt
    if length(data.time) > 1
        config.dt = mean(diff(data.time));
    else
        config.dt = 1.0;
    end
    
    % Prepare inputs for simulation
    inputs.time = data.time;
    inputs.T_c_in = data.T_in;
    inputs.I_st = data.I_st;
    inputs.v_c_matrix = data.V_cell;
    % Placeholder for flow distribution
    inputs.q_c_matrix = repmat(data.flow_rate / config.N, 1, config.N); 
    
    % Optimization
    % Params: [C_ed, C_cell, R_cell_ed, R_cell, hconv]
    lb = [1000, 100, 0.001, 0.01, 5];
    ub = [5000, 1000, 0.1, 0.5, 100];
    
    options.MaxIter = 50;
    options.SwarmSize = 30;
    options.Display = 'iter';
    
    obj_func = @(p) cost_function(p, data, inputs, config);
    
    fprintf('Starting Thermal Parameter Identification via PSO...\n');
    [best_p, best_cost] = SimplePSO(obj_func, 5, lb, ub, options);
    
    identified_params.C_ed = best_p(1);
    identified_params.C_cell = best_p(2);
    identified_params.R_cell_ed = best_p(3);
    identified_params.R_cell = best_p(4);
    identified_params.hconv = best_p(5);
    
    fprintf('Identification Complete.\n');
    disp(identified_params);
end

function cost = cost_function(p, data, inputs, config)
    % Pack params
    params.C_ed = p(1);
    params.C_cell = p(2);
    params.R_cell_ed = p(3);
    params.R_cell = p(4);
    params.hconv = p(5);
    
    % Run Simulation
    % Initial state
    if isfield(data, 'T_dist') && ~isempty(data.T_dist)
        initial_state = [data.T_dist(1, 1); data.T_dist(1, :)'; data.T_dist(1, end)]; 
    else
        initial_state = ones(config.N + 2, 1) * config.Tenv;
    end
    
    [T_cell_sim, ~, ~] = ThermalModel(params, initial_state, inputs, config);
    
    % Cost Calculation
    if isfield(data, 'T_dist') && ~isempty(data.T_dist)
        diff = T_cell_sim - data.T_dist;
        cost = sum(diff(:).^2);
    else
        cost = 1e9; 
    end
end
