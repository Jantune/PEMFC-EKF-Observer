function optimized_qr = Optimize_QR_Matrices()
    % Optimize_QR_Matrices
    % Optimizes Q and R matrices for EKF using Global PSO.
    
    current_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile(current_dir, '../Data'));
    addpath(fullfile(current_dir, '../Models'));
    addpath(fullfile(current_dir, '../EKF'));
    
    data = LoadData();
    stack_config = LoadStackConfig();
    config = stack_config;
    
    if isempty(data.time)
        fprintf('No data available. Returning default QR params.\n');
        optimized_qr = config.defaults;
        return;
    end
    
    % Determine N
    if isfield(data, 'V_cell') && ~isempty(data.V_cell)
        config.N = size(data.V_cell, 2);
    else
        config.N = 180;
    end
    
    % Determine dt
    if length(data.time) > 1
        config.dt = mean(diff(data.time));
    else
        config.dt = 1.0;
    end
    
    params = config.defaults;
    
    % Bounds for [q_std_cell, q_std_ed, r_var_meas]
    lb = [0.1, 0.1, 0.001];
    ub = [10, 10, 1.0];
    
    options.MaxIter = 30;
    options.SwarmSize = 20;
    options.Display = 'iter';
    
    obj_func = @(p) cost_function(p, data, params, config);
    
    fprintf('Starting QR Optimization via PSO...\n');
    [best_p, best_cost] = SimplePSO(obj_func, 3, lb, ub, options);
    
    optimized_qr.q_std_cell = best_p(1);
    optimized_qr.q_std_ed = best_p(2);
    optimized_qr.r_var_measurement = best_p(3);
    
    fprintf('QR Optimization Complete.\n');
    disp(optimized_qr);
end

function cost = cost_function(p, data, params, config)
    q_std_cell = p(1);
    q_std_ed = p(2);
    r_var = p(3);
    
    % Construct Matrices
    Q = diag([q_std_ed^2; repmat(q_std_cell^2, config.N, 1); q_std_ed^2]);
    R = r_var;
    
    % Initialize EKF
    if isfield(data, 'T_dist') && ~isempty(data.T_dist)
        initial_X = [data.T_dist(1, 1); data.T_dist(1, :)'; data.T_dist(1, end)];
    else
        initial_X = ones(config.N+2, 1) * config.Tenv;
    end
    
    initial_P = eye(config.N+2);
    
    ekf = ExtendedKalmanFilter(initial_X, initial_P, Q, R, params, config);
    
    % Run EKF
    n_steps = length(data.time);
    rmse_accum = 0;
    count = 0;
    
    q_total_vec = data.flow_rate;
    q_dist_mat = repmat(q_total_vec / config.N, 1, config.N);
    
    for k = 1:n_steps
        inputs.q_c = q_dist_mat(k, :)';
        inputs.v_c = data.V_cell(k, :)';
        inputs.T_c_in = data.T_in(k);
        inputs.I_st = data.I_st(k);
        
        ekf.Predict(inputs);
        
        y_meas = data.T_out(k);
        ekf.Update(y_meas, q_total_vec(k), inputs.q_c);
        
        if isfield(data, 'T_dist') && ~isempty(data.T_dist)
             est_T = ekf.X(2:end-1);
             true_T = data.T_dist(k, :)';
             rmse_accum = rmse_accum + mean((est_T - true_T).^2);
             count = count + 1;
        end
    end
    
    if count > 0
        cost = sqrt(rmse_accum / count);
    else
        cost = 1e9;
    end
end
