function results = TrainResidualBaseline(data, method, train_time_indices, val_time_index)
    % TrainResidualBaseline
    % Trains a hybrid physics-ML model (Residual Correction).
    %
    % Inputs:
    %   data: Data struct from LoadData
    %   method: 'ridge' or 'gpr'
    %   train_time_indices: (Optional) Vector of time indices for training snapshots
    %   val_time_index: (Optional) Time index for validation snapshot
    %
    % If train/val indices not provided, attempts to use specific times from paper
    % (635s, 1349s for training, 1708s for validation) or falls back to generic split.
    
    if nargin < 2; method = 'ridge'; end
    
    % Load Stack Config
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
    
    % Parameters
    params = config.defaults;
    
    inputs.time = data.time;
    inputs.T_c_in = data.T_in;
    inputs.I_st = data.I_st;
    inputs.v_c_matrix = data.V_cell;
    inputs.q_c_matrix = repmat(data.flow_rate / config.N, 1, config.N); 
    
    initial_state = ones(config.N+2, 1) * config.Tenv;
    
    [T_phys_sim, ~, ~] = ThermalModel(params, initial_state, inputs, config);
    
    % Determine Training and Validation Indices
    if nargin < 3 || isempty(train_time_indices)
        % Try paper-specific times first
        train_t1 = find(data.time == 635, 1);
        train_t2 = find(data.time == 1349, 1);
        
        if ~isempty(train_t1) && ~isempty(train_t2)
            train_indices = [train_t1, train_t2];
        else
            % Generic fallback: Use sparse snapshots from first 70% of data
            n_total = length(data.time);
            train_region = 1:floor(0.7 * n_total);
            % Select 2-3 evenly spaced snapshots
            if length(train_region) >= 3
                step = floor(length(train_region) / 3);
                train_indices = train_region([step, 2*step]);
            else
                train_indices = train_region;
            end
            fprintf('Paper snapshot times not found. Using generic training indices: [%s]\n', num2str(train_indices));
        end
    else
        train_indices = train_time_indices;
    end
    
    if nargin < 4 || isempty(val_time_index)
        % Try paper-specific validation time
        val_t = find(data.time == 1708, 1);
        
        if ~isempty(val_t)
            val_index = val_t;
        else
            % Generic fallback: Use a snapshot from last 30% of data
            n_total = length(data.time);
            val_region = floor(0.7 * n_total):n_total;
            val_index = val_region(round(length(val_region)/2));
            fprintf('Paper validation time not found. Using generic validation index: %d\n', val_index);
        end
    else
        val_index = val_time_index;
    end
    
    % Check if ground truth available
    if ~isfield(data, 'T_dist') || isempty(data.T_dist)
        fprintf('Ground truth T_dist not available. Cannot train ML baseline.\n');
        results.metrics.rmse = NaN;
        results.metrics.mae = NaN;
        results.metrics.max_rel_pct = NaN;
        return;
    end
    
    X_train = [];
    y_train = [];
    
    % Feature Engineering
    for t_idx = train_indices
        if t_idx < 1 || t_idx > size(T_phys_sim, 1)
            warning('Training index %d out of bounds. Skipping.', t_idx);
            continue;
        end
        
        T_phys_t = T_phys_sim(t_idx, :)';
        T_exp_t = data.T_dist(t_idx, :)';
        
        for i = 1:config.N
            feat = [T_phys_t(i), data.I_st(t_idx), data.T_in(t_idx), data.V_cell(t_idx, i), i/config.N];
            res = T_exp_t(i) - T_phys_t(i);
            
            X_train = [X_train; feat];
            y_train = [y_train; res];
        end
    end
    
    if isempty(X_train)
        fprintf('No valid training data collected.\n');
        results.metrics.rmse = NaN;
        return;
    end
    
    % Standardization
    mu = mean(X_train);
    sig = std(X_train);
    sig(sig==0) = 1;
    X_train_std = (X_train - mu) ./ sig;
    
    % Train Model
    if strcmpi(method, 'ridge')
        lam = 1.0;
        X_aug = [ones(size(X_train_std, 1), 1), X_train_std];
        w = (X_aug' * X_aug + lam * eye(size(X_aug, 2))) \ (X_aug' * y_train);
        
        predict_func = @(X) [ones(size(X,1), 1), X] * w;
        
    elseif strcmpi(method, 'gpr')
        try
            mdl = fitrgp(X_train_std, y_train, 'KernelFunction', 'squaredexponential');
            predict_func = @(X) predict(mdl, X);
        catch
            fprintf('GPR training failed (toolbox missing?). Using Ridge fallback.\n');
            lam = 1.0;
            X_aug = [ones(size(X_train_std, 1), 1), X_train_std];
            w = (X_aug' * X_aug + lam * eye(size(X_aug, 2))) \ (X_aug' * y_train);
            predict_func = @(X) [ones(size(X,1), 1), X] * w;
        end
    end
    
    % Validation
    if val_index < 1 || val_index > size(T_phys_sim, 1)
        warning('Validation index %d out of bounds.', val_index);
        results.metrics.rmse = NaN;
        return;
    end
    
    T_phys_t = T_phys_sim(val_index, :)';
    T_exp_t = data.T_dist(val_index, :)';
    
    X_val = [];
    for i = 1:config.N
        feat = [T_phys_t(i), data.I_st(val_index), data.T_in(val_index), data.V_cell(val_index, i), i/config.N];
        X_val = [X_val; feat];
    end
    X_val_std = (X_val - mu) ./ sig;
    
    pred_res = predict_func(X_val_std);
    T_final_pred = T_phys_t + pred_res;
    
    % Metrics
    err = T_final_pred - T_exp_t;
    results.metrics.rmse = sqrt(mean(err.^2));
    results.metrics.mae = mean(abs(err));
    results.metrics.max_rel_pct = max(abs(err) ./ max(abs(T_exp_t), eps)) * 100;
    
    fprintf('Training on %d snapshots, validation at time index %d\n', length(train_indices), val_index);
    fprintf('RMSE: %.4f, MAE: %.4f\n', results.metrics.rmse, results.metrics.mae);
end
