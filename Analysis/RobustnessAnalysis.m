function RobustnessAnalysis()
    % RobustnessAnalysis
    % Analyzes EKF sensitivity to parameter uncertainties.
    
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
    
    % Nominal Params
    nominal_params = config.defaults;
    
    q_std_cell = config.defaults.q_std_cell;
    q_std_ed = config.defaults.q_std_ed;
    r_var = config.defaults.r_var;
    
    param_names = {'C_ed', 'C_cell', 'R_cell_ed', 'R_cell', 'hconv'};
    perturbations = [-0.2, -0.1, 0.1, 0.2]; % +/- 10%, 20%
    
    fprintf('=== Robustness Analysis ===\n');
    
    % Baseline RMSE
    rmse_base = RunEKF_GetRMSE(nominal_params, config, data, q_std_cell, q_std_ed, r_var);
    fprintf('Baseline RMSE: %.4f\n', rmse_base);
    
    for p_idx = 1:length(param_names)
        p_name = param_names{p_idx};
        fprintf('\nPerturbing %s:\n', p_name);
        
        for pert = perturbations
            % Perturb
            params_pert = nominal_params;
            val = params_pert.(p_name);
            new_val = val * (1 + pert);
            params_pert.(p_name) = new_val;
            
            % Run
            rmse = RunEKF_GetRMSE(params_pert, config, data, q_std_cell, q_std_ed, r_var);
            
            % Report
            change = (rmse - rmse_base) / rmse_base * 100;
            fprintf('  %+d%% (%.4f) -> RMSE: %.4f (%.2f%% change)\n', ...
                    pert*100, new_val, rmse, change);
        end
    end
end

function rmse = RunEKF_GetRMSE(params, config, data, q_std_cell, q_std_ed, r_var)
    % Run EKF and calculate RMSE
    
    Q = diag([q_std_ed^2; repmat(q_std_cell^2, config.N, 1); q_std_ed^2]);
    R = r_var;
    
    initial_X = ones(config.N+2, 1) * config.Tenv;
    ekf = ExtendedKalmanFilter(initial_X, eye(config.N+2), Q, R, params, config);
    
    n_steps = length(data.time);
    q_dist_mat = repmat(data.flow_rate / config.N, 1, config.N);
    
    rmse_sum = 0;
    count = 0;
    
    for k = 1:n_steps
        inputs.q_c = q_dist_mat(k, :)';
        inputs.v_c = data.V_cell(k, :)';
        inputs.T_c_in = data.T_in(k);
        inputs.I_st = data.I_st(k);
        
        ekf.Predict(inputs);
        ekf.Update(data.T_out(k), data.flow_rate(k), inputs.q_c);
        
        if isfield(data, 'T_dist') && ~isempty(data.T_dist)
            est = ekf.X(2:end-1);
            truth = data.T_dist(k, :)';
            rmse_sum = rmse_sum + mean((est - truth).^2);
            count = count + 1;
        end
    end
    
    if count > 0
        rmse = sqrt(rmse_sum / count);
    else
        rmse = NaN;
    end
end
