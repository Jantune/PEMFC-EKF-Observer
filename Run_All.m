% Run_All.m
% Orchestrates the PEMFC EKF Observer Project Execution

clc; clear; close all;

addpath('Data');
addpath('Models');
addpath('EKF');
addpath('Optimization');
addpath('Analysis');
addpath('ML_Baselines');

fprintf('======================================================\n');
fprintf('PEMFC Cell-by-Cell Temperature Observer Project\n');
fprintf('======================================================\n\n');

% 1. Data Check
fprintf('1. Checking Data...\n');
data = LoadData();
if isempty(data.time)
    fprintf('   [WARNING] Data is empty. Most scripts will exit early or use defaults.\n');
    fprintf('   Please populate PEMFC-EKF-Observer/Data/LoadData.m with experimental data.\n\n');
else
    % Display data info
    if isfield(data, 'V_cell') && ~isempty(data.V_cell)
        N_detected = size(data.V_cell, 2);
        fprintf('   Data loaded: %d cells, %d time steps\n', N_detected, length(data.time));
    end
    if length(data.time) > 1
        dt_detected = mean(diff(data.time));
        fprintf('   Time step: %.2f s\n', dt_detected);
    end
end

% 2. Parameter Identification
fprintf('\n2. Parameter Identification (PSO)...\n');
fprintf('   [SKIPPED] To run identification, uncomment lines in Run_All.m\n');
% Uncomment to run (might be slow without data)
% hydraulic_params = Identify_Hydraulic_Parameters();
% thermal_params = Identify_Thermal_Parameters();

% 3. QR Optimization
fprintf('\n3. QR Matrix Optimization (PSO)...\n');
fprintf('   [SKIPPED] Using defaults from LoadStackConfig.m\n');
% Uncomment to run optimization
% optimized_qr = Optimize_QR_Matrices();

% 4. EKF Estimation
fprintf('\n4. Running EKF Estimation...\n');
try
    Run_EKF_Estimation();
catch ME
    fprintf('   [ERROR] EKF failed: %s\n', ME.message);
end

% 5. Time Performance Analysis
fprintf('\n5. Time Performance Analysis...\n');
try
    run('Analysis/TimePerformanceAnalysis.m');
catch ME
    fprintf('   [ERROR] Performance analysis failed: %s\n', ME.message);
end

% 6. ML Baseline Comparison
fprintf('\n6. ML Baseline Comparison...\n');
try
    Compare_ML_Models();
catch ME
    fprintf('   [ERROR] ML comparison failed: %s\n', ME.message);
end

% 7. Robustness Analysis
fprintf('\n7. Robustness Analysis...\n');
try
    RobustnessAnalysis();
catch ME
    fprintf('   [ERROR] Robustness analysis failed: %s\n', ME.message);
end

fprintf('\n======================================================\n');
fprintf('All tasks completed.\n');
fprintf('======================================================\n');
