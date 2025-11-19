function [T_cell_sim, T_ed1_sim, T_ed2_sim] = ThermalModel(params, initial_state, inputs, config)
    % ThermalModel Runs the 1D thermal simulation over a time horizon.
    %
    % Inputs:
    %   params: Struct with fields C_ed, C_cell, R_cell_ed, R_cell, hconv
    %   initial_state: [T_ed1; T_cell(1:N); T_ed2] (N+2 x 1)
    %   inputs: Struct with fields:
    %           time (T x 1)
    %           q_c_matrix (T x N) - Coolant flow distribution
    %           v_c_matrix (T x N) - Cell voltages
    %           T_c_in (T x 1)
    %           I_st (T x 1)
    %   config: Struct with fields:
    %           N (int)
    %           E0, Tenv, Aed (constants)
    %           dt (step size)
    %
    % Outputs:
    %   T_cell_sim: Simulated cell temperatures (T x N)
    %   T_ed1_sim: Simulated End Plate 1 temperatures (T x 1)
    %   T_ed2_sim: Simulated End Plate 2 temperatures (T x 1)

    N = config.N;
    dt = config.dt;
    E0 = config.E0;
    Tenv = config.Tenv;
    Aed = config.Aed;
    
    % Unpack parameters
    C_ed = params.C_ed;
    C_cell = params.C_cell;
    R_cell_ed = params.R_cell_ed;
    R_cell = params.R_cell;
    hconv = params.hconv;
    
    % Unpack inputs
    steps = length(inputs.time);
    
    % Preallocate outputs
    T_cell_sim = zeros(steps, N);
    T_ed1_sim = zeros(steps, 1);
    T_ed2_sim = zeros(steps, 1);
    
    % Initialize state
    X_current = initial_state;
    
    T_ed1_sim(1) = X_current(1);
    T_cell_sim(1, :) = X_current(2:N+1)';
    T_ed2_sim(1) = X_current(N+2);
    
    % Simulation Loop
    for k = 1:steps-1
        % Get inputs for current step
        q_c_t = inputs.q_c_matrix(k, :)'; % column vector
        v_c_t = inputs.v_c_matrix(k, :)';
        T_c_in_t = inputs.T_c_in(k);
        I_st_t = inputs.I_st(k);
        
        % Run one step
        [X_next, ~, ~] = RunOneStep(X_current, C_ed, C_cell, R_cell_ed, R_cell, hconv, ...
                                    N, q_c_t, v_c_t, T_c_in_t, I_st_t, E0, Tenv, Aed, dt);
                                
        % Update state
        X_current = X_next;
        
        % Store results
        T_ed1_sim(k+1) = X_current(1);
        T_cell_sim(k+1, :) = X_current(2:N+1)';
        T_ed2_sim(k+1) = X_current(N+2);
    end
end

function [X_next, T_ed1_next, T_ed2_next] = RunOneStep(X_current, C_ed, C_cell, R_cell_ed, R_cell, hconv, N, q_c_t, v_c_t, T_c_in_t, I_st_t, E0, Tenv, Aed, dt)
    % Logic adapted from original run_physical_model_onestep.m
    
    T_ed1_current = X_current(1);
    T_cell_current = X_current(2:N+1); 
    T_ed2_current = X_current(N+2);

    T_cell_next_val = zeros(N, 1);

    % Coolant properties (Water/Glycol mixture typical, or pure water approximation)
    % Original code used polynomials based on T_c_in_t
    current_temp = T_c_in_t;
    
    % Density (kg/m^3)
    rou = -1.5629e-5 * current_temp^3 + 0.011778 * current_temp^2 - 3.0726 * current_temp + 1227.8;
    
    % Specific Heat (J/kgK)
    % Original polynomial resulted in kJ/kgK approx, then multiply by 1000 to get J/kgK
    % capacity calculation from original code:
    capacity_kJ = (1.1105e-5 * current_temp^3 - 3.1078e-3 * current_temp^2 - 1.4787 * current_temp + 4631.9) / 1000;
    capacity = capacity_kJ; % Keep it as is, handled in Q calc
    
    % Factor for heat removal: q(L/min) * (1e-3/60) -> m^3/s * rho * cp * 1000 -> W
    % Original code: (q_c_t * 0.001 / 60) * rou * capacity * (deltaT) * 10^3;
    % The 10^3 suggests capacity was in kJ/kgK, converting final result to W.
    
    factor = (0.001 / 60) * rou * capacity * 1000;

    % Cell 1
    Q_coolant_1 = q_c_t(1) * factor * (T_cell_current(1) - T_c_in_t);
    Q1 = I_st_t * (E0 - v_c_t(1)) - Q_coolant_1;
    dT_cell_1 = ((T_ed1_current - T_cell_current(1)) / R_cell_ed + (T_cell_current(2) - T_cell_current(1)) / R_cell + Q1) / C_cell;
    T_cell_next_val(1) = T_cell_current(1) + dT_cell_1 * dt;
    
    % Cells 2 to N-1
    % Vectorized
    Q_coolant = q_c_t(2:N-1) .* factor .* (T_cell_current(2:N-1) - T_c_in_t);
    Q_gen = I_st_t * (E0 - v_c_t(2:N-1));
    Q_net = Q_gen - Q_coolant;
    
    dT_conduction = (T_cell_current(1:N-2) - T_cell_current(2:N-1)) / R_cell + ...
                    (T_cell_current(3:N) - T_cell_current(2:N-1)) / R_cell;
                    
    dT_cell = (dT_conduction + Q_net) / C_cell;
    T_cell_next_val(2:N-1) = T_cell_current(2:N-1) + dT_cell * dt;

    % Cell N
    Q_coolant_N = q_c_t(N) * factor * (T_cell_current(N) - T_c_in_t);
    QN = I_st_t * (E0 - v_c_t(N)) - Q_coolant_N;
    dT_cell_N = ((T_ed2_current - T_cell_current(N)) / R_cell_ed + (T_cell_current(N-1) - T_cell_current(N)) / R_cell + QN) / C_cell;
    T_cell_next_val(N) = T_cell_current(N) + dT_cell_N * dt;

    % End Plates
    dT_ed_1 = ((T_cell_current(1) - T_ed1_current) / R_cell_ed - hconv * Aed * (T_ed1_current - Tenv)) / C_ed;
    T_ed1_next = T_ed1_current + dT_ed_1 * dt;
    
    dT_ed_2 = ((T_cell_current(N) - T_ed2_current) / R_cell_ed - hconv * Aed * (T_ed2_current - Tenv)) / C_ed;
    T_ed2_next = T_ed2_current + dT_ed_2 * dt;
    
    X_next = [T_ed1_next; T_cell_next_val; T_ed2_next];
end
