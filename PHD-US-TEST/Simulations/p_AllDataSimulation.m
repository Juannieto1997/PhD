%% p_AllDataSimulation
% Runs the K-wave simulation in for all the frequencies and values in the
% the CSV file. 
% Stores a .mat file with a table containing the frequency, driving
% voltage, intensity and pressure. 
%% Prepare workspace 
% clear workspace 
clc; close all; clear all; beep off; 
% Add necesary folders 
dxf_domains_dir ='.\Domains2\'; % Data for the different medium
addpath('C:\Users\Juan\Documents\GitHub\PhD\FUNCTIONS') % Functions folder
addpath('C:\Users\Juan\Documents\MATLAB\k-Wave')
%% Create constant variables 
%CSV file location
str_CSV = 'C:\Users\Juan\Documents\SecondBrain\PhD\Resources\Transducer Characterization\20240729\ScaleMeassurements.csv';
str_SaveName = 'PressureSimulation.mat'; % name of the file to save 

% Properties 
FUS_on_time = 400e-6;       % (BHTE) set source on time -> Pulse duration [s]
FUS_off_time = 20;        % (BHTE) set source off time [s] -> Post-stim cooling simulated duration
axial_tx_offset = 0;        % m -> Vertical transducer offset
freefield = false;

% ---- Tx properties ----

tx_aperture = 15e-3;                        % Tx aperture [m]
tx_radius = 15e-3;                          % Spherical Tx radius of curvature [m]
tx_surface_area = pi * (tx_aperture / 2)^2; % Surface area [m2]

alpha_water = 0.0022;                       % Water attenuation [dB/(MHz^y cm)]
rho_water = 994.04;                         % Water density [kg.m-3]
c_water = 1482.3;                           % Speed of sound in water at 20°C [m.s-1]

% ---- Domain size ----

zdim = 20e-3;                               % axial domain dimension [m]
xdim = 10e-3;                               % radial domain dimension [m]

npts_per_lambda = 9.5;                      % N points per wavelenght
pml_size = 3;                              % PML size [grid points]
%% Load CSV and prepare data for analysis
m_Data = readmatrix(str_CSV);

v_Freq = m_Data(:,1); % Frequency tested 
v_Driv = m_Data(:,2); % Driving voltage
v_Inte = m_Data(:,11);% Intensity at the transducer 
v_Pres = zeros(1,length(v_Freq));

%% run the simulation for each value
for idx = 1:length(v_Freq)
    %% update properties of the simulation
    tx_freq = v_Freq(idx)*1000000;             % Hz
    tx_emitted_ac_power = v_Inte(idx);   % W -> I went with a ridiculously high value just so that you can see clear temperature elevation a 1MHz
    if tx_emitted_ac_power == 0
        v_Pres(idx) = 0; 
        continue 
    end 

    lambda_water = c_water / tx_freq;           % Ultrasound wavelength [m]
    % number of grid points in the axial (x) and radial (y) directions
    Nz = round(npts_per_lambda * zdim / lambda_water) + 2 * pml_size;
    Nx = round(npts_per_lambda * xdim / lambda_water) + pml_size;

    % grid point spacing along x and y [m]
    dz = lambda_water / npts_per_lambda;
    dx = lambda_water / npts_per_lambda;
    kgrid = kWaveGrid(Nz, dz, Nx, dx);
    % create plot axis
    z_vec_AS = (kgrid.x_vec - kgrid.x_vec(1)) - pml_size*dz;
    x_vec_AS = (kgrid.y_vec - kgrid.y_vec(1));

    % Medium properties
    medium.sound_speed = c_water * ones(Nz, Nx);        % [m/s]
    medium.density = rho_water * ones(Nz, Nx);          % [kg/m^3]
    alpha_coeff_map = alpha_water * ones(Nz, Nx);       % [dB/(MHz^y cm)]
    alpha_power_map = ones(Nz, Nx);
    if freefield ~= true
    
        % Air bottom domain
        air1_dxf = join([dxf_domains_dir 'Air1.dxf']);
        air1_bm = dxf2binary_mask(air1_dxf, Nz, Nx, dz, dx, pml_size*dz);
        medium.sound_speed(air1_bm) = 343;              % [m/s]
        medium.density(air1_bm) = 1.16;                 % [kg/m^3]
        alpha_coeff_map(air1_bm) = 0.0034;              % [dB/(MHz^y cm)]
        alpha_power_map(air1_bm) = 2;                   %

        % Air upper domain
        air2_dxf = join([dxf_domains_dir 'Air2.dxf']);
        air2_bm = dxf2binary_mask(air2_dxf, Nz, Nx, dz, dx, pml_size*dz);
        medium.sound_speed(air2_bm) = 343;              % [m/s]
        medium.density(air2_bm) = 1.16;                 % [kg/m^3]
        alpha_coeff_map(air2_bm) = 0.0034;              % [dB/(MHz^y cm)]
        alpha_power_map(air2_bm) = 2;                   %

        % Agarose coupling gel (Yang2017 / Culjat2010)
        agarose_gel_dxf = join([dxf_domains_dir 'AgaroseGel.dxf']);
        agarose_gel_bm = dxf2binary_mask(agarose_gel_dxf, Nz, Nx, dz, dx, pml_size*dz);
        medium.sound_speed(agarose_gel_bm) = 1503;      % [m/s]
        medium.density(agarose_gel_bm) = 1050;          % [kg/m^3]
        alpha_coeff_map(agarose_gel_bm) = 0.014;        % [dB/(MHz^y cm)]
        alpha_power_map(agarose_gel_bm) = 1;            %

        % Coupling cone (PLA parameters from Ma2022)
        coupling_cone_dxf = join([dxf_domains_dir 'CouplingCone.dxf']);
        coupling_cone_bm = dxf2binary_mask(coupling_cone_dxf, Nz, Nx, dz, dx, pml_size*dz);
        medium.sound_speed(coupling_cone_bm) = 1250;    % [m/s]
        medium.density(coupling_cone_bm) = 1240;        % [kg/m^3]

        % Petri dish coverglass
        coverglass_dxf = join([dxf_domains_dir 'Coverglass.dxf']);
        coverglass_bm = dxf2binary_mask(coverglass_dxf, Nz, Nx, dz, dx, pml_size*dz + axial_tx_offset);
        medium.sound_speed(coverglass_bm) = 5570;       % [m/s]
        medium.density(coverglass_bm) = 2230;           % [kg/m^3]
        alpha_coeff_map(coverglass_bm) = 5;             % [dB/(MHz^y cm)]        DOES NOT CONVERGE WITH 16.2
        alpha_power_map(coverglass_bm) = 1;             %

        % PDMS
        pdms_dxf = join([dxf_domains_dir 'Brainslice.dxf']);
        pdms_bm = dxf2binary_mask(pdms_dxf, Nz, Nx, dz, dx, pml_size*dz + axial_tx_offset);
        medium.sound_speed(pdms_bm) = 1054;             % [m/s]
        medium.density(pdms_bm) = 1038;                 % [kg/m^3]
        alpha_coeff_map(pdms_bm) = 1.907;               % [dB/(MHz^y cm)]
        alpha_power_map(pdms_bm) = 1.46;                %
    end

    medium.alpha_mode ='stokes';
    medium.alpha_power = 2;
    medium.alpha_coeff = alpha_coeff_map .* ((tx_freq*1e-6).^alpha_power_map / (tx_freq*1e-6).^medium.alpha_power); % [dB/(MHz^y cm)]

    %% ---- Pressure source geometrical definition ----

    % Convert Tx aperture in grid point
    if rem(round(tx_aperture / kgrid.dx), 2)
        tx_aperture_asgrid = round(tx_aperture / kgrid.dx);
    else
        tx_aperture_asgrid = round(tx_aperture / kgrid.dx) + 1;
    end

    % Tx geometry
    tx_pole_grid_loc = [pml_size+1, 1];
    tx_radius_asgrid = round(tx_radius / kgrid.dx);
    tx_focus_grid_loc = tx_pole_grid_loc + [tx_radius_asgrid, 0];
    % Define the position of the source (Tx) in the mask with a shape of an arc
    source.p_mask = makeArc([kgrid.Nx, kgrid.Ny], tx_pole_grid_loc, tx_radius_asgrid, tx_aperture_asgrid, tx_focus_grid_loc);

    % Cell location retreival
    if freefield == true
        % Takes focus location if freefield sim
        [M, cell_location_index] = min(abs(z_vec_AS - tx_radius), [], "all", "linear");
    else
        % Takes PDMS upper interface location if pdms_bm is defined
        [M, cell_location_index] = max(diff(pdms_bm(:, 1)), [], "all", "linear");
    end

    %% Domain vizualisation
    % [C, ia, ic] = unique(medium.sound_speed);
    % dom_viz = single(reshape(ic, size(medium.sound_speed))) .* ~single(source.p_mask);

    % figure
    % imagesc(x_vec_AS, z_vec_AS, dom_viz);
    % yline(z_vec_AS(cell_location_index), 'r', LineWidth=2)
    % xlabel('x (radial) position [m]');
    % ylabel('z (axial) position [m]');
    % axis image
    if freefield == true
        CFL = 0.14;
    else
        CFL = 0.04;
    end

    n_reflections = 2;
    t_end = (Nz * dz) * n_reflections / c_water;  % last time [s]
    % k-wave function to create the time array in the kgrid object
    kgrid.makeTime(max(medium.sound_speed), CFL, t_end);
    kgrid % Uncomment to preview grid infos
    % ---- Pressure source ----

    % Tx coupling medium -> grab most common values at the Tx location
    % (ignore values resulting of rasterization errors -> eg coupling cone PLA properties)
    c_tx_coupling_medium = mode(medium.sound_speed(logical(source.p_mask)));
    rho_tx_coupling_medium = mode(medium.density(logical(source.p_mask)));

    % Generate array of continuous wave (CW) signals from amplitude and phase
    tx_p_mag = sqrt(2) * sqrt((tx_emitted_ac_power * rho_tx_coupling_medium * c_tx_coupling_medium) / tx_surface_area); % Pressure [Pa]
    source.p = createCWSignals(kgrid.t_array, tx_freq, tx_p_mag, 0);

    % define all point og the computational grip as a sensor
    sensor.mask = ones(kgrid.Nx, kgrid.Ny);
    sensor.record = {'p', 'u'};
    % record the last cycles in steady state
    T = 1 / tx_freq;
    num_periods = 3;
    T_points = round(num_periods * T / kgrid.dt);
    sensor.record_start_index = kgrid.Nt - T_points + 1;
    % assign the input options
    input_args = { ...
        'PlotSim', true, ... % real time sim plot
        'DataCast', 'single', ...
        'DisplayMask', source.p_mask, ... % dom_viz, ...
        'PMLSize', pml_size, ...
        'PMLInside', false, ...
        'PlotPML', true, ...
        'PlotScale', 10 * [-1, 1] * tx_p_mag ...
        };

    sensor_data = kspaceFirstOrderAS(kgrid, medium, source, sensor, input_args{:});

    % extract the pressure magnitude at each position
    p_mag = extractAmpPhase(sensor_data.p, 1/kgrid.dt, tx_freq); % Get pressure modulus
    p_mag = reshape(p_mag, Nz, Nx);
    p_mag_mirrored = [p_mag(:, end:-1:2) p_mag];
    p_mag_z_profile = squeeze(p_mag(:, 1));
    x_vec_mirrored = [-x_vec_AS(end:-1:2); x_vec_AS];
    
    [~,loc]=min(abs(z_vec_AS - 0.0150));
       
    v_Pres(idx) = p_mag_z_profile(loc);


    % Axial pressure pressure profile preview (modulus)
    % figure
    % hold on
    % plot(z_vec_AS, p_mag_z_profile) % Modulus
    % xlabel('z (axial or acoustic axis) position [m]');
    % ylabel('Pressure magnitude [Pa]');

    %% Close all figures for next plot 
    close all;

end 


save(str_SaveName,'v_Freq','v_Driv','v_Inte','v_Pres')