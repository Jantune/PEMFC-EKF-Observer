function GenerateSampleData()
    % GenerateSampleData
    % Generates a synthetic dataset for testing the PEMFC Temperature Observer.
    % This data mimics the structure and trends described in the paper 
    % "Cell-by-Cell Temperature Observer for Commercial PEMFC Stacks".
    %
    % The generated data includes:
    % - Time vector
    % - Stack Current (Dynamic profile)
    % - Flow Rate (Variable)
    % - Voltages (Cell-wise)
    % - Inlet/Outlet Temperatures
    % - Pressures
    %
    % Saves to 'PEMFC-EKF-Observer/Data/sample_data.mat'

    fprintf('Generating synthetic sample data...\n');

    % Configuration
    N = 180; % Number of cells
    dt = 1.0; % Time step (s)
    T_duration = 2000; % Total duration (s)
    time = (0:dt:T_duration-1)';
    
    % 1. Current Profile (Mimicking Fig 5a: Preheat, Startup, Dynamic Load)
    I_st = zeros(T_duration, 1);
    % Preheating (0-1060s): I=0
    % Startup (1060-1300s): Ramp to 50A
    % Dynamic (1300-end): Steps [100, 180, 240, 300]
    
    idx_start = 1240;
    I_st(idx_start:1300) = linspace(0, 50, 1300-idx_start+1);
    I_st(1301:1500) = 100;
    I_st(1501:1700) = 180;
    I_st(1701:1900) = 240;
    I_st(1901:end) = 300;
    
    % Add some noise
    I_st = I_st + randn(size(I_st)) * 0.5;
    I_st(I_st < 0) = 0;

    % 2. Flow Rate Profile (Mimicking Fig 5b)
    % Base flow related to current (stoichiometry) + Preheating flow
    % Preheat: 20 L/min
    % Interruption: 0 L/min at ~1100s
    
    flow_rate = zeros(T_duration, 1);
    flow_rate(1:1086) = 20;
    flow_rate(1087:1208) = 0; % Coolant interruption fault
    
    % Recovery and Load Following
    target_flow = 20 + (I_st / 300) * 40; % Rough scaling
    
    % Smooth transition for recovery
    flow_rate(1209:end) = target_flow(1209:end);
    
    % 3. Inlet Temperature (Mimicking Fig 5b)
    % Ramps up from 25C to ~60C
    T_in = zeros(T_duration, 1);
    T_in(1:1000) = linspace(25, 60, 1000);
    T_in(1001:end) = 60 + randn(1000, 1) * 0.2;
    
    % 4. Cell Voltages
    % V = E0 - R*I - eta (Simplified polarization)
    % N cells, slight variation
    V_cell = zeros(T_duration, N);
    E_ocv = 0.95;
    
    for i = 1:N
        % Slight parameter variation per cell
        R_ohm = 0.0003 * (1 + 0.1*randn()); 
        V_cell(:, i) = E_ocv - R_ohm * I_st - 0.05 * log(I_st + 1);
    end
    
    % 5. Ground Truth Temperature (Simulate using a "True" Model)
    % We use the same simplified physics model to generate "True" data,
    % but maybe with slightly different parameters to simulate model mismatch.
    
    true_params.C_ed = 1800;
    true_params.C_cell = 330;
    true_params.R_cell_ed = 0.0055;
    true_params.R_cell = 0.095;
    true_params.hconv = 22;
    
    config.N = N;
    config.E0 = 1.48;
    config.Tenv = 25;
    config.Aed = 0.1;
    config.dt = dt;
    
    % Run thermal simulation
    inputs.time = time;
    inputs.q_c_matrix = repmat(flow_rate / N, 1, N); % Uniform for generation
    inputs.v_c_matrix = V_cell;
    inputs.T_c_in = T_in;
    inputs.I_st = I_st;
    
    initial_state = ones(N+2, 1) * 25;
    
    % We need the ThermalModel function available. 
    % Assuming we are in Data/, need to add path to Models/
    if exist('../Models', 'dir')
        addpath('../Models');
    else
        % Try relative to file
        p_data = fileparts(mfilename('fullpath'));
        addpath(fullfile(p_data, '../Models'));
    end
    
    % Call Thermal Model
    try
        [T_cell_sim, ~, ~] = ThermalModel(true_params, initial_state, inputs, config);
        T_dist = T_cell_sim;
    catch
        warning('ThermalModel not found or failed. Generating random T_dist.');
        T_dist = 60 + randn(T_duration, N);
    end
    
    % 6. Outlet Temperature (Measured)
    % Weighted average of T_dist
    T_out = mean(T_dist, 2) + randn(T_duration, 1) * 0.1; % Add sensor noise
    
    % 7. Pressures (Based on Flow)
    % Simple quadratic relation P = k * q^2
    P_in = 100 + 0.1 * flow_rate.^2; 
    P_out = 100 * ones(T_duration, 1);
    
    % Save
    sample_data.time = time;
    sample_data.flow_rate = flow_rate;
    sample_data.P_in = P_in;
    sample_data.P_out = P_out;
    sample_data.T_in = T_in;
    sample_data.I_st = I_st;
    sample_data.V_cell = V_cell;
    sample_data.T_dist = T_dist; % Ground truth
    sample_data.T_out = T_out;   % Measurement
    
    save('sample_data.mat', 'sample_data');
    fprintf('Sample data saved to sample_data.mat\n');
end

