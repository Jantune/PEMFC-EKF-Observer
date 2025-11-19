function data = LoadData()
    % LoadData Loads experimental data for the project.
    %
    % Returns a struct with fields:
    %   time: Time vector (s)
    %   flow_rate: Total coolant flow rate (L/min)
    %   P_in: Inlet pressure (kPa)
    %   P_out: Outlet pressure (kPa)
    %   T_in: Inlet temperature (degC)
    %   I_st: Stack current (A)
    %   V_cell: Cell voltages (matrix N x T)
    %   T_dist: Temperature distribution (matrix N x T or similar) - Optional (Ground Truth)
    %   T_out: Outlet temperature (degC) - Measured
    
    % Try to load sample data if it exists
    sample_file = fullfile(fileparts(mfilename('fullpath')), 'sample_data.mat');
    
    if exist(sample_file, 'file')
        loaded = load(sample_file);
        data = loaded.sample_data;
        fprintf('Loaded sample data from %s\n', sample_file);
    else
        fprintf('Sample data not found. Generating now...\n');
        try
            GenerateSampleData();
            loaded = load(sample_file);
            data = loaded.sample_data;
            fprintf('Generated and loaded sample data.\n');
        catch ME
             fprintf('Error generating data: %s\n', ME.message);
             data.time = [];
        end
    end
    
    % NOTE: To use your own data, comment out the above block and populate 'data' manually:
    % data.time = ...; 
    % ...
end
