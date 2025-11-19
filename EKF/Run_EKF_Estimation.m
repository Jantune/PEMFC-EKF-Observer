function Run_EKF_Estimation()
    % Run_EKF_Estimation
    % Main script to run the Temperature Estimation Observer.
    
    current_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile(current_dir, '../Data'));
    addpath(fullfile(current_dir, '../Models'));
    addpath(fullfile(current_dir, '../EKF'));
    addpath(fullfile(current_dir, '../Optimization'));
    
    % 1. Load Data
    data = LoadData();
    if isempty(data.time)
        fprintf('No data. Exiting.\n');
        return;
    end
    
    % 2. Load Stack Configuration
    stack_config = LoadStackConfig();
    
    % Dynamic Configuration from Data
    config = stack_config; % Inherit constants
    
    % Determine N
    if isfield(data, 'V_cell') && ~isempty(data.V_cell)
        config.N = size(data.V_cell, 2);
    else
        config.N = 180; % Fallback
        fprintf('Warning: N could not be determined from V_cell. Using default %d.\n', config.N);
    end
    
    % Determine dt
    if length(data.time) > 1
        config.dt = mean(diff(data.time));
    else
        config.dt = 1.0; % Fallback
    end
    
    % 3. Parameters (Thermal/Hydraulic)
    % Load identified if available, else use defaults from config
    params = config.defaults;
    
    % 4. QR Matrices (Optimized or Defaults)
    q_std_cell = config.defaults.q_std_cell;
    q_std_ed = config.defaults.q_std_ed;
    r_var = config.defaults.r_var;
    
    Q = diag([q_std_ed^2; repmat(q_std_cell^2, config.N, 1); q_std_ed^2]);
    R = r_var;
    
    % 5. Initialize EKF
    % Initial state: Ambient temp or from data
    initial_X = ones(config.N+2, 1) * config.Tenv; 
    initial_P = eye(config.N+2);
    
    ekf = ExtendedKalmanFilter(initial_X, initial_P, Q, R, params, config);
    
    % 6. Run Loop
    n_steps = length(data.time);
    results.X_est = zeros(n_steps, config.N+2);
    
    fprintf('Running EKF Estimation (N=%d, dt=%.2fs)...\n', config.N, config.dt);
    
    % Hydraulic Pre-calculation (Using Hydraulic Model + Defaults)
    % Ideally this should be dynamic if pressures vary significantly
    hydraulic_params = [params.alpha, params.beta];
    
    % Use P_in/P_out if available, else use flow-based lookup or simple distribution
    % For speed in this loop, we often use pre-calculated distribution weights
    % Since HydraulicModel is available, use it if pressure data is available.
    % Note: HydraulicModel is an iterative solver, might be slow inside the loop.
    % Standard EKF approach: Calculate distribution weights W based on Flow.
    % If flow changes, distribution changes slightly.
    % Simplified: Assumed uniform or pre-calibrated map.
    % Let's use the HydraulicModel once to get a distribution profile or assume uniform if P not available.
    
    use_hydraulic_model = isfield(data, 'P_in') && ~isempty(data.P_in);
    
    if use_hydraulic_model
         fprintf('Using Hydraulic Model for flow distribution...\n');
    else
         fprintf('Pressure data missing. Using Uniform Flow Distribution.\n');
    end
    
    tic;
    for k = 1:n_steps
        % Current inputs
        total_flow = data.flow_rate(k);
        
        if use_hydraulic_model && total_flow > 0
            % This might be computationally expensive for every step. 
            % Optimization: Interpolation table or only update when flow changes significantly.
            % For now, we call it directly to ensure physics fidelity.
            [~, q_dist] = HydraulicModel(hydraulic_params, data.P_in(k), data.P_out(k), config.N, stack_config);
        else
            q_dist = ones(config.N, 1) * (total_flow / config.N);
        end
        
        inputs.q_c = q_dist;
        inputs.v_c = data.V_cell(k, :)';
        inputs.T_c_in = data.T_in(k);
        inputs.I_st = data.I_st(k);
        
        % Predict
        ekf.Predict(inputs);
        
        % Update
        y_meas = data.T_out(k);
        ekf.Update(y_meas, total_flow, inputs.q_c);
        
        results.X_est(k, :) = ekf.X';
    end
    t_run = toc;
    
    fprintf('Estimation Complete. Time: %.2f s\n', t_run);
end
