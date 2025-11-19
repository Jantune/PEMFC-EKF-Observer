% TimePerformanceAnalysis.m
% Benchmarking Jacobian calculation methods (Analytical vs Numerical)
% Part of PEMFC-EKF-Observer Project

clear; clc;

current_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(current_dir, '../EKF'));
addpath(fullfile(current_dir, '../Models'));
addpath(fullfile(current_dir, '../Data'));

fprintf('========================================\n');
fprintf('Jacobian Computation Performance Benchmark\n');
fprintf('========================================\n');

%% 1. Setup
% Load Configuration
stack_config = LoadStackConfig();

% Use default N for benchmark (can be overridden by loading data)
config = stack_config;
config.N = 180; % Benchmark N (representative commercial stack size)
config.dt = 1;

params = config.defaults;

% Initial State
N_states = config.N + 2;
X_test = config.Tenv + randn(N_states, 1) * 2;
P_test = eye(N_states);
Q = eye(N_states);
R = 1;

fprintf('Testing with N = %d cells (state dimension: %d)\n', config.N, N_states);

% Instantiate EKF
ekf = ExtendedKalmanFilter(X_test, P_test, Q, R, params, config);

% Test Inputs
inputs_test.q_c = ones(config.N, 1) * 0.001;
inputs_test.v_c = ones(config.N, 1) * 0.85;
inputs_test.T_c_in = 23;
inputs_test.I_st = 50;

fprintf('Test setup complete.\n\n');

%% 2. Analytical Method (Class Method)
fprintf('>>> Testing Analytical Method (Tridiagonal)...\n');
num_runs = 100;
time_analytical = zeros(num_runs, 1);

% Warmup
for i = 1:10
    T_prop = inputs_test.T_c_in;
    rho = -1.5629e-5 * T_prop^3 + 0.011778 * T_prop^2 - 3.0726 * T_prop + 1227.8;
    cp = (1.1105e-5 * T_prop^3 - 3.1078e-3 * T_prop^2 - 1.4787 * T_prop + 4631.9) / 1000;
    factor = rho * cp * 1000 * (0.001 / 60);
    
    F = ekf.CalculateJacobian(inputs_test.q_c, rho, cp, factor);
end

for i = 1:num_runs
    tic;
    T_prop = inputs_test.T_c_in;
    rho = -1.5629e-5 * T_prop^3 + 0.011778 * T_prop^2 - 3.0726 * T_prop + 1227.8;
    cp = (1.1105e-5 * T_prop^3 - 3.1078e-3 * T_prop^2 - 1.4787 * T_prop + 4631.9) / 1000;
    factor = rho * cp * 1000 * (0.001 / 60);
    F_analytical = ekf.CalculateJacobian(inputs_test.q_c, rho, cp, factor);
    time_analytical(i) = toc;
end

mean_analytical = mean(time_analytical) * 1000;
fprintf('  Mean Time: %.4f ms\n', mean_analytical);
fprintf('  Non-zero elements: %d / %d (density: %.2f%%)\n', ...
        nnz(F_analytical), N_states^2, nnz(F_analytical)/N_states^2*100);

%% 3. Numerical Method (Finite Difference)
fprintf('\n>>> Testing Numerical Method (Finite Difference)...\n');

% Helper for numerical Jacobian
function F_num = numerical_jacobian(ekf_obj, inputs)
    delta = 1e-6;
    N_s = length(ekf_obj.X);
    F_num = zeros(N_s, N_s);
    X_orig = ekf_obj.X;
    
    ekf_obj.X = X_orig;
    ekf_obj.Predict(inputs);
    X_next_base = ekf_obj.X;
    ekf_obj.X = X_orig;
    
    for j = 1:N_s
        X_pert = X_orig;
        X_pert(j) = X_pert(j) + delta;
        
        ekf_obj.X = X_pert;
        ekf_obj.Predict(inputs);
        X_next_pert = ekf_obj.X;
        ekf_obj.X = X_orig;
        
        F_num(:, j) = (X_next_pert - X_next_base) / delta;
    end
end

time_numerical = zeros(num_runs, 1);
% Warmup
for i = 1:5
    F_num = numerical_jacobian(ekf, inputs_test);
end

for i = 1:num_runs
    tic;
    F_numerical = numerical_jacobian(ekf, inputs_test);
    time_numerical(i) = toc;
end

mean_numerical = mean(time_numerical) * 1000;
fprintf('  Mean Time: %.4f ms\n', mean_numerical);
fprintf('  Speedup: %.2fx\n\n', mean_numerical / mean_analytical);

%% 4. ECU Projection
fprintf('>>> ECU Projection (Assuming 10-15x slowdown)...\n');
slowdown_factors = [10, 12, 15];

for sf = slowdown_factors
    ecu_analytical = mean_analytical * sf;
    fprintf('  %dx slower: Analytical = %.2f ms (%.2f%% of 1s sampling period)\n', ...
            sf, ecu_analytical, ecu_analytical / 10);
end

fprintf('\n========================================\n');
fprintf('Performance Summary\n');
fprintf('========================================\n');
fprintf('Analytical (Tridiagonal): %.4f ms\n', mean_analytical);
fprintf('Numerical (Dense):        %.4f ms\n', mean_numerical);
fprintf('Speedup Factor:           %.1fx\n', mean_numerical / mean_analytical);
fprintf('Memory Efficiency:        %.1fx (sparse vs dense)\n', N_states^2 / nnz(F_analytical));
fprintf('========================================\n');
