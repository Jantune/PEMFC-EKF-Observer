function Compare_ML_Models()
    % Compare_ML_Models
    % Compares Physics-Based EKF with Data-Driven Hybrid Models.
    
    current_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile(current_dir, '../Data'));
    addpath(fullfile(current_dir, '../Models'));
    addpath(fullfile(current_dir, '../EKF'));
    
    data = LoadData();
    if isempty(data.time)
        fprintf('No data. Exiting.\n');
        return;
    end
    
    stack_config = LoadStackConfig();
    config = stack_config;
    
    % Determine N and dt
    if isfield(data, 'V_cell') && ~isempty(data.V_cell)
        config.N = size(data.V_cell, 2);
    else
        config.N = 180;
    end
    
    if length(data.time) > 1
        config.dt = mean(diff(data.time));
    else
        config.dt = 1.0;
    end
    
    fprintf('=== ML Baseline Comparison ===\n');
    
    % Determine validation snapshot index
    % Try paper-specific time (1708s) first, else use last 20% of data
    val_idx = find(data.time == 1708, 1);
    if isempty(val_idx)
        n_total = length(data.time);
        val_idx = floor(0.9 * n_total); % Use 90% point
        fprintf('Using validation at time index %d (%.1fs)\n', val_idx, data.time(val_idx));
    else
        fprintf('Using paper validation time: 1708s\n');
    end
    
    % 1. Train and Evaluate Ridge Regression Hybrid
    fprintf('\n--- Ridge Regression Hybrid ---\n');
    ridge_results = TrainResidualBaseline(data, 'ridge');
    
    % 2. Train and Evaluate GPR Hybrid
    fprintf('\n--- GPR Hybrid ---\n');
    gpr_results = TrainResidualBaseline(data, 'gpr');
    
    % 3. Evaluate EKF
    fprintf('\n--- EKF Observer ---\n');
    
    params = config.defaults;
    
    q_std_cell = config.defaults.q_std_cell;
    q_std_ed = config.defaults.q_std_ed;
    r_var = config.defaults.r_var;
    Q = diag([q_std_ed^2; repmat(q_std_cell^2, config.N, 1); q_std_ed^2]);
    R = r_var;
    
    initial_X = ones(config.N+2, 1) * config.Tenv;
    ekf = ExtendedKalmanFilter(initial_X, eye(config.N+2), Q, R, params, config);
    
    % Run simulation up to validation point
    n_steps = length(data.time);
    q_dist_mat = repmat(data.flow_rate / config.N, 1, config.N);
    
    target_idx = min(val_idx, n_steps);
    
    for k = 1:target_idx
        inputs.q_c = q_dist_mat(k, :)';
        inputs.v_c = data.V_cell(k, :)';
        inputs.T_c_in = data.T_in(k);
        inputs.I_st = data.I_st(k);
        
        ekf.Predict(inputs);
        ekf.Update(data.T_out(k), data.flow_rate(k), inputs.q_c);
    end
    
    X_ekf = ekf.X(2:end-1); % Cells
    
    % Compare
    if isfield(data, 'T_dist') && ~isempty(data.T_dist)
         T_true = data.T_dist(target_idx, :)';
         err = X_ekf - T_true;
         rmse_ekf = sqrt(mean(err.^2));
         fprintf('EKF RMSE (time idx %d): %.4f\n', target_idx, rmse_ekf);
    else
         rmse_ekf = NaN;
         fprintf('Ground truth not available for EKF validation.\n');
    end
    
    % Summary Table
    fprintf('\n=== Summary ===\n');
    fprintf('Method\t\tRMSE (degC)\n');
    fprintf('Ridge\t\t%.4f\n', ridge_results.metrics.rmse);
    fprintf('GPR\t\t%.4f\n', gpr_results.metrics.rmse);
    fprintf('EKF\t\t%.4f\n', rmse_ekf);
end
