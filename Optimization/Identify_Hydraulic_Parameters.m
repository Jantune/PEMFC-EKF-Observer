function identified_params = Identify_Hydraulic_Parameters()
    % Identify_Hydraulic_Parameters
    % Identifies alpha and beta for the hydraulic model using Global PSO.
    
    % Add paths
    current_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile(current_dir, '../Data'));
    addpath(fullfile(current_dir, '../Models'));
    
    % Load Data
    data = LoadData();
    
    if isempty(data.flow_rate)
        fprintf('No data available for identification. Please populate Data/LoadData.m\n');
        identified_params = [0, 0];
        return;
    end
    
    % Load Config
    stack_config = LoadStackConfig();
    
    % Optimization configuration
    lb = [1e4, 1e3]; % Lower bounds for alpha, beta
    ub = [1e8, 1e6]; % Upper bounds for alpha, beta
    
    options.MaxIter = 50;
    options.SwarmSize = 30;
    options.Display = 'iter';
    
    % Objective Function
    obj_func = @(params) cost_function(params, data, stack_config);
    
    % Run PSO
    fprintf('Starting Hydraulic Parameter Identification via PSO...\n');
    [best_params, best_cost] = SimplePSO(obj_func, 2, lb, ub, options);
    
    alpha = best_params(1);
    beta = best_params(2);
    
    fprintf('Identification Complete.\n');
    fprintf('Alpha: %.4f\n', alpha);
    fprintf('Beta: %.4f\n', beta);
    fprintf('Final Cost (SSE): %.4f\n', best_cost);
    
    identified_params = best_params;
end

function cost = cost_function(params, data, config)
    % Calculate Sum of Squared Errors between Measured and Predicted Total Flow
    
    q_meas = data.flow_rate;
    P_in = data.P_in;
    P_out = data.P_out;
    
    % Determine N dynamically if possible
    if isfield(data, 'V_cell') && ~isempty(data.V_cell)
        N = size(data.V_cell, 2);
    else
        N = 180; % Fallback
    end
    
    n_points = length(q_meas);
    total_error = 0;
    
    for i = 1:n_points
        if q_meas(i) > 0 % Only consider non-zero flow points
            [q_pred, ~] = HydraulicModel(params, P_in(i), P_out(i), N, config);
            total_error = total_error + (q_meas(i) - q_pred)^2;
        end
    end
    
    cost = total_error;
end
