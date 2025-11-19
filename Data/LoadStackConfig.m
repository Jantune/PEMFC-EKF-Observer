function config = LoadStackConfig()
    % LoadStackConfig
    % Returns the geometric and physical configuration of the stack.
    % This allows modifying stack constants in one place.
    
    % Geometric / Hydraulic Constants (Used in HydraulicModel)
    config.L = 2e-3;            % Length of main channel segment per cell (m)
    config.Am = 15e-4;          % Main channel cross-sectional area (m^2)
    config.z = 10;              % Number of channels per cell
    config.Ach = 6e-6;          % Channel cross-sectional area (m^2)
    config.Mu = 464e-3;         % Dynamic viscosity (N*s/m^2)
    config.rou = 982;           % Density (kg/m^3)
    config.D_hm = 3.75e-2;      % Hydraulic diameter (m)
    config.K_fm_in = 0.1;       % Resistance coefficient (manifold inlet)
    config.K_fm_out = 0.5;      % Resistance coefficient (manifold outlet)
    config.R_c_in = 1.6;        % Resistance coefficient (cell inlet)
    config.R_c_out = 2;         % Resistance coefficient (cell outlet)
    
    % Thermal / Electrochemical Constants
    config.E0 = 1.48;           % Thermo-neutral voltage (V)
    config.Tenv = 25;           % Default Ambient Temperature (degC)
    config.Aed = 0.1;           % End plate area (m^2)
    
    % Default Initial Parameters (Identified)
    % These serve as defaults if identification is not run.
    config.defaults.C_ed = 1811.6;
    config.defaults.C_cell = 335.6;
    config.defaults.R_cell_ed = 0.005;
    config.defaults.R_cell = 0.09;
    config.defaults.hconv = 20;
    
    config.defaults.alpha = 1e5; % Hydraulic alpha
    config.defaults.beta = 1e4;  % Hydraulic beta
    
    % Default QR (for EKF)
    config.defaults.q_std_cell = 0.5;
    config.defaults.q_std_ed = 0.5;
    config.defaults.r_var = 0.05;
end

