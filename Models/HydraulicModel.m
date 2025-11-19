function [q_total_pred, q_distribution] = HydraulicModel(params, P_in, P_out, N, config)
    % HydraulicModel calculates total coolant flow and distribution
    % based on inlet/outlet pressures and model parameters.
    %
    % Inputs:
    %   params: [alpha, beta] - Hydraulic resistance coefficients
    %   P_in:  Inlet pressure (kPa)
    %   P_out: Outlet pressure (kPa)
    %   N:     Number of cells
    %   config: Stack configuration struct (LoadStackConfig)
    %
    % Outputs:
    %   q_total_pred: Predicted total flow rate (L/min)
    %   q_distribution: Flow rate for each cell (L/min) (N x 1)

    alpha = params(1);
    beta = params(2);
    
    % Constants from config
    L = config.L;
    Am = config.Am;
    z = config.z;
    Ach = config.Ach;
    Mu = config.Mu;
    rou = config.rou;
    D_hm = config.D_hm;
    K_fm_in = config.K_fm_in;
    K_fm_out = config.K_fm_out;
    R_c_in = config.R_c_in;
    R_c_out = config.R_c_out;

    % Derived constants for equation: A*q^2 + B*q + C = 0
    AA = (R_c_in + R_c_out + 2*alpha/rou) * rou / (2 * z^2 * Ach^2);
    BB = beta / (z * Ach);
    
    delta_P_target = (P_in - P_out) * 1e3; % Convert to Pa
    
    % Iterative solver to find q_c(N) that matches total pressure drop
    
    function P_drop_calc = calculate_pressure_drop(q_N_guess)
        q_c = zeros(N, 1);
        q_m = zeros(N, 1);
        deltP_m_in = zeros(N, 1);
        deltP_m_out = zeros(N, 1);
        
        q_c(N) = q_N_guess;
        
        CC = 0;
        for k = N:-1:2
             % Pressure drop in cell k
             deltP_cell_k = (R_c_in + R_c_out) * rou * q_c(k)^2 / (2 * z^2 * Ach^2) + ...
                            alpha * (q_c(k) / (z * Ach))^2 + ...
                            beta * q_c(k) / (z * Ach);
            
             % Flow in manifold segments
             if k == N
                 q_m(k) = q_c(k); 
             else
                 q_m(k) = q_m(k+1) + q_c(k);
             end
             
             % Reynolds number
             Re_k = (rou * D_hm * (q_m(k) / Am)) / Mu;
             if Re_k < 2100
                 fm_k = 64 / max(Re_k, 1e-5); 
             else
                 fm_k = 0.316 / max(Re_k, 1e-5)^0.25;
             end
             
             % Manifold pressure drops
             deltP_m_in(k) = (fm_k * L / D_hm + K_fm_in) * rou * q_m(k)^2 / (2 * Am^2);
             deltP_m_out(k) = (fm_k * L / D_hm + K_fm_out) * rou * q_m(k)^2 / (2 * Am^2);
             
             if k == N
                 CC = deltP_m_in(k) + deltP_m_out(k);
             else
                 CC = deltP_m_in(k) + deltP_m_out(k) - deltP_m_in(k+1) - deltP_m_out(k+1); 
             end
             
             % Solve for q_c(k-1)
             % AA*q^2 + BB*q - (CC + deltP_cell_k) = 0
             q_c(k-1) = (-BB + sqrt(BB^2 + 4 * AA * (CC + deltP_cell_k))) / (2*AA);
        end
        
        % Calculate total pressure drop for cell 1 (loop ends at k=2, computing q_c(1))
        q_m(1) = q_m(2) + q_c(1);
        Re_1 = (rou * D_hm * (q_m(1) / Am)) / Mu;
        if Re_1 < 2100
            fm_1 = 64 / max(Re_1, 1e-5);
        else
             fm_1 = 0.316 / max(Re_1, 1e-5)^0.25;
        end
        
        deltP_cell_1 = (R_c_in + R_c_out) * rou * q_c(1)^2 / (2 * z^2 * Ach^2) + ...
                       alpha * (q_c(1) / (z * Ach))^2 + ...
                       beta * q_c(1) / (z * Ach);
                       
        deltP_m_in_1 = (fm_1 * L / D_hm + K_fm_in) * rou * q_m(1)^2 / (2 * Am^2);
        deltP_m_out_1 = (fm_1 * L / D_hm + K_fm_out) * rou * q_m(1)^2 / (2 * Am^2);
        
        P_drop_calc = deltP_cell_1 + deltP_m_in_1 + deltP_m_out_1;
    end

    % Simple bisection or fzero to find q_N
    err_func = @(q_val) calculate_pressure_drop(q_val) - delta_P_target;
    
    options = optimset('TolX', 1e-8, 'Display', 'off');
    try
        q_N_sol = fzero(err_func, [1e-7, 5e-4], options);
    catch
        q_N_sol = fzero(err_func, 1e-5, options);
    end
    
    % Re-run to get distribution (identical logic)
    q_c_final = zeros(N, 1);
    q_c_final(N) = q_N_sol;
    
    q_m = zeros(N, 1);
    deltP_m_in = zeros(N, 1);
    deltP_m_out = zeros(N, 1);
    CC = 0;
    for k = N:-1:2
         deltP_cell_k = (R_c_in + R_c_out) * rou * q_c_final(k)^2 / (2 * z^2 * Ach^2) + ...
                        alpha * (q_c_final(k) / (z * Ach))^2 + ...
                        beta * q_c_final(k) / (z * Ach);
         if k == N
             q_m(k) = q_c_final(k);
         else
             q_m(k) = q_m(k+1) + q_c_final(k);
         end
         Re_k = (rou * D_hm * (q_m(k) / Am)) / Mu;
         if Re_k < 2100
             fm_k = 64 / max(Re_k, 1e-5);
         else
             fm_k = 0.316 / max(Re_k, 1e-5)^0.25;
         end
         deltP_m_in(k) = (fm_k * L / D_hm + K_fm_in) * rou * q_m(k)^2 / (2 * Am^2);
         deltP_m_out(k) = (fm_k * L / D_hm + K_fm_out) * rou * q_m(k)^2 / (2 * Am^2);
         if k == N
             CC = deltP_m_in(k) + deltP_m_out(k);
         else
             CC = deltP_m_in(k) + deltP_m_out(k) - deltP_m_in(k+1) - deltP_m_out(k+1); 
         end
         q_c_final(k-1) = (-BB + sqrt(BB^2 + 4 * AA * (CC + deltP_cell_k))) / (2*AA);
     end
     
     q_distribution = q_c_final * 60 * 1000;
     q_total_pred = sum(q_distribution);
end
