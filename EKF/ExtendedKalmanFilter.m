classdef ExtendedKalmanFilter < handle
    % ExtendedKalmanFilter Implementation of EKF for PEMFC Temp Estimation
    % Optimized with Analytical Tridiagonal Jacobian.
    
    properties
        N       % Number of cells
        dt      % Time step
        X       % State vector [T_ed1; T_cell(1:N); T_ed2] (N+2 x 1)
        P       % Covariance matrix
        Q       % Process noise covariance
        R       % Measurement noise covariance
        
        % Model parameters
        params
        config
    end
    
    methods
        function obj = ExtendedKalmanFilter(initial_X, initial_P, Q, R, params, config)
            obj.X = initial_X;
            obj.P = initial_P;
            obj.Q = Q;
            obj.R = R;
            obj.params = params;
            obj.config = config;
            obj.N = config.N;
            obj.dt = config.dt;
        end
        
        function Predict(obj, inputs)
            % Predict step
            % inputs: struct with q_c (vector), v_c (vector), T_c_in, I_st
            
            % Unpack for readability
            C_ed = obj.params.C_ed;
            C_cell = obj.params.C_cell;
            R_cell_ed = obj.params.R_cell_ed;
            R_cell = obj.params.R_cell;
            hconv = obj.params.hconv;
            
            T_ed1 = obj.X(1);
            T_cell = obj.X(2:obj.N+1);
            T_ed2 = obj.X(obj.N+2);
            
            q_c = inputs.q_c;
            v_c = inputs.v_c;
            T_in = inputs.T_c_in;
            I = inputs.I_st;
            
            E0 = obj.config.E0;
            Tenv = obj.config.Tenv;
            Aed = obj.config.Aed;
            dt = obj.dt;
            
            % --- Physics Update ---
            % Coolant props
            T_prop = T_in;
            rho = -1.5629e-5 * T_prop^3 + 0.011778 * T_prop^2 - 3.0726 * T_prop + 1227.8;
            cp = (1.1105e-5 * T_prop^3 - 3.1078e-3 * T_prop^2 - 1.4787 * T_prop + 4631.9) / 1000;
            factor = rho * cp * 1000 * (0.001 / 60); % q(L/min) -> Q(W) factor
            
            T_next = zeros(obj.N+2, 1);
            
            % Cell 1
            Q_cool_1 = q_c(1) * factor * (T_cell(1) - T_in);
            Q_gen_1 = I * (E0 - v_c(1));
            dT1 = ((T_ed1 - T_cell(1))/R_cell_ed + (T_cell(2) - T_cell(1))/R_cell + Q_gen_1 - Q_cool_1) / C_cell;
            T_next(2) = T_cell(1) + dT1 * dt;
            
            % Mid Cells
            if obj.N > 2
                Q_cool = q_c(2:end-1) .* factor .* (T_cell(2:end-1) - T_in);
                Q_gen = I .* (E0 - v_c(2:end-1));
                dT_cond = (T_cell(1:end-2) - T_cell(2:end-1))/R_cell + (T_cell(3:end) - T_cell(2:end-1))/R_cell;
                dT = (dT_cond + Q_gen - Q_cool) / C_cell;
                T_next(3:obj.N) = T_cell(2:end-1) + dT * dt;
            end
            
            % Cell N
            Q_cool_N = q_c(obj.N) * factor * (T_cell(obj.N) - T_in);
            Q_gen_N = I * (E0 - v_c(obj.N));
            dTN = ((T_ed2 - T_cell(obj.N))/R_cell_ed + (T_cell(obj.N-1) - T_cell(obj.N))/R_cell + Q_gen_N - Q_cool_N) / C_cell;
            T_next(obj.N+1) = T_cell(obj.N) + dTN * dt;
            
            % End Plates
            dT_ed1 = ((T_cell(1) - T_ed1)/R_cell_ed - hconv * Aed * (T_ed1 - Tenv)) / C_ed;
            T_next(1) = T_ed1 + dT_ed1 * dt;
            
            dT_ed2 = ((T_cell(obj.N) - T_ed2)/R_cell_ed - hconv * Aed * (T_ed2 - Tenv)) / C_ed;
            T_next(obj.N+2) = T_ed2 + dT_ed2 * dt;
            
            X_pred = T_next;
            
            % 2. Jacobian Calculation (Analytical)
            F = obj.CalculateJacobian(q_c, rho, cp, factor);
            
            % 3. Covariance Prediction
            obj.P = F * obj.P * F' + obj.Q;
            obj.X = X_pred;
        end
        
        function Update(obj, y_meas, q_c_total, q_c_dist)
            % Update step
            % y_meas: Measured stack outlet temperature
            
            if q_c_total > 1e-6
                w = q_c_dist / q_c_total;
            else
                w = ones(obj.N, 1) / obj.N; % Fallback
            end
            
            H = zeros(1, obj.N+2);
            H(2:obj.N+1) = w';
            
            % Innovation
            y_pred = H * obj.X;
            res = y_meas - y_pred;
            
            % Gain
            S = H * obj.P * H' + obj.R;
            K = (obj.P * H') / S;
            
            % State Update
            obj.X = obj.X + K * res;
            
            % Covariance Update
            obj.P = (eye(obj.N+2) - K * H) * obj.P;
        end
        
        function F = CalculateJacobian(obj, q_c, rho, cp, factor)
            % Analytical Tridiagonal Jacobian
            
            N_states = obj.N + 2;
            
            % Diagonals
            d_main = zeros(N_states, 1);
            d_upper = zeros(N_states-1, 1);
            d_lower = zeros(N_states-1, 1);
            
            C_ed = obj.params.C_ed;
            C_cell = obj.params.C_cell;
            R_cell_ed = obj.params.R_cell_ed;
            R_cell = obj.params.R_cell;
            hconv = obj.params.hconv;
            Aed = obj.config.Aed;
            dt = obj.dt;
            
            gamma = q_c * factor;
            
            % End Plate 1 (Index 1)
            d_main(1) = 1 + dt * (-1/R_cell_ed - hconv*Aed) / C_ed;
            d_upper(1) = dt * (1/R_cell_ed) / C_ed;
            
            % Cell 1 (Index 2)
            d_lower(1) = dt * (1/R_cell_ed) / C_cell;
            d_main(2) = 1 + dt * (-1/R_cell_ed - 1/R_cell - gamma(1)) / C_cell;
            d_upper(2) = dt * (1/R_cell) / C_cell;
            
            % Middle Cells (Indices 3 to N) -> Cell i corresponds to State i+1
            for i = 2:obj.N-1
                idx = i + 1;
                d_lower(idx-1) = dt * (1/R_cell) / C_cell;
                d_main(idx) = 1 + dt * (-2/R_cell - gamma(i)) / C_cell;
                d_upper(idx) = dt * (1/R_cell) / C_cell;
            end
            
            % Cell N (Index N+1)
            idx = obj.N + 1;
            d_lower(idx-1) = dt * (1/R_cell) / C_cell;
            d_main(idx) = 1 + dt * (-1/R_cell_ed - 1/R_cell - gamma(obj.N)) / C_cell;
            d_upper(idx) = dt * (1/R_cell_ed) / C_cell;
            
            % End Plate 2 (Index N+2)
            idx = obj.N + 2;
            d_lower(idx-1) = dt * (1/R_cell_ed) / C_ed;
            d_main(idx) = 1 + dt * (-1/R_cell_ed - hconv*Aed) / C_ed;
            
            % Construct Sparse Matrix
            F = spdiags([ [d_lower; 0], d_main, [0; d_upper] ], -1:1, N_states, N_states);
        end
    end
end
