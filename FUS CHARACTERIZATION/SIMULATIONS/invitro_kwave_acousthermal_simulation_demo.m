clearvars;
clear all;
close all;
beep off;
clc;

dxf_domains_dir ='.\Domains2\'; % Data for the different medium
addpath('C:\Users\Juan\Documents\GitHub\PhD\FUNCTIONS') % Functions folder
addpath('C:\Users\Juan\Documents\MATLAB\k-Wave')
% =========================================================================
% GENERAL SIMULATION PARAMETERS
% =========================================================================

tx_freq = 2.24e6;              % Hz
tx_emitted_ac_power = 14;   % W -> I went with a ridiculously high value just so that you can see clear temperature elevation a 1MHz
FUS_on_time = 400e-6;       % (BHTE) set source on time -> Pulse duration [s]
FUS_off_time = 1;        % (BHTE) set source off time [s] -> Post-stim cooling simulated duration
axial_tx_offset = 0;        % m -> Vertical transducer offset
freefield = true;

disp('=========================================================================');
disp(join(['Running sim at ' sprintf('%0.0e', tx_freq) 'Hz, ' sprintf('%0.0e', tx_emitted_ac_power) 'W acoustic power' newline 'Tx z offset: ' sprintf('%0.0e', axial_tx_offset) 'm'], ''))
disp('=========================================================================');

% =========================================================================
% ACOUSTIC SIMULATION
% =========================================================================

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
lambda_water = c_water / tx_freq;           % Ultrasound wavelength [m]
npts_per_lambda = 9.5;                      % N points per wavelenght
pml_size = 20;                              % PML size [grid points]

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
    medium.sound_speed(pdms_bm) = 1500;             % [m/s]
    medium.density(pdms_bm) = 1044.5;                 % [kg/m^3]
    alpha_coeff_map(pdms_bm) = 0.104;               % [dB/(MHz^y cm)]
    alpha_power_map(pdms_bm) = 0.104;                % [dB/(MHz^y cm)]
end

medium.alpha_mode ='stokes';
medium.alpha_power = 2;
medium.alpha_coeff = alpha_coeff_map .* ((tx_freq*1e-6).^alpha_power_map / (tx_freq*1e-6).^medium.alpha_power); % [dB/(MHz^y cm)]

% ---- Pressure source geometrical definition ----

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

% Domain vizualisation
[C, ia, ic] = unique(medium.sound_speed);
dom_viz = single(reshape(ic, size(medium.sound_speed))) .* ~single(source.p_mask);

figure
imagesc(x_vec_AS, z_vec_AS, dom_viz);
yline(z_vec_AS(cell_location_index), 'r', LineWidth=2)
xlabel('x (radial) position [m]');
ylabel('z (axial) position [m]');
axis image

% %% ---- Time array ----

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

% Axial pressure pressure profile preview (modulus)
figure
hold on
plot(z_vec_AS, p_mag_z_profile) % Modulus
xlabel('z (axial or acoustic axis) position [m]');
ylabel('Pressure magnitude [Pa]');

% Pressure field preview (modulus)
figure
imagesc(x_vec_mirrored, z_vec_AS, p_mag_mirrored);
yline(z_vec_AS(cell_location_index), 'r', LineWidth=2)
xlabel('x (radial) position [m]');
ylabel('z (axial or acoustic axis) position [m]');
axis image

% Focal Spot
% find line closer to the focal spot 
v_z = abs(z_vec_AS - 0.0150);
[~,s_idx] = min(v_z);
v_centerP = p_mag(s_idx,:);

n = length(v_centerP);
theta = linspace(0,2*pi,n).';
x0 = linspace(0,x_vec_AS(end),n);
y0 = v_centerP;
x = x0.*cos(theta);
z = x0.*sin(theta);
y = repmat(y0,[n,1]);

figure;
surf(x,z,y,'FaceAlpha',0.85, 'EdgeColor', 'none')



% =========================================================================
%% BHTE EVALUATION
% =========================================================================

% ---- BHTE Subdomain construction ----

% Temperature profile at cell location
bhte_domain_height = 6e-3; % m
bhte_domain_width = 8e-3; % m

% bhte_domain_height = 12e-3; % m
% bhte_domain_width = 10e-3; % m

cells_z_vec = z_vec_AS - z_vec_AS(cell_location_index);
bhte_domain_x_mask = abs(cells_z_vec) < bhte_domain_height / 2;
bhte_domain_y_mask = abs(x_vec_mirrored) < bhte_domain_width / 2;
bhte_Nx = sum(bhte_domain_x_mask);
bhte_Ny = sum(bhte_domain_y_mask);

% Pressure field preview (modulus)
figure
imagesc(x_vec_mirrored(bhte_domain_y_mask), cells_z_vec(bhte_domain_x_mask), p_mag_mirrored(bhte_domain_x_mask, bhte_domain_y_mask));
xlabel('x (radial) position [m]');
ylabel('z (axial) position [m]');
axis image

% ---- Domain size ----
bhte_kgrid = kWaveGrid(bhte_Nx, kgrid.dx, bhte_Ny, kgrid.dy);
y_mirrored_vec = bhte_kgrid.y_vec;

% --- Material properties definition ---

bhte_medium.density = [medium.density(:, end:-1:2) medium.density];                 % [kg/m^3]
bhte_medium.density = bhte_medium.density(bhte_domain_x_mask, bhte_domain_y_mask);  % BHTE subdomain

bhte_medium.sound_speed = [medium.sound_speed(:, end:-1:2) medium.sound_speed];             % [m/s]
bhte_medium.sound_speed = bhte_medium.sound_speed(bhte_domain_x_mask, bhte_domain_y_mask);  % BHTE subdomain

bhte_medium.alpha_coeff = [medium.alpha_coeff(:, end:-1:2) medium.alpha_coeff];
bhte_medium.alpha_coeff = bhte_medium.alpha_coeff(bhte_domain_x_mask, bhte_domain_y_mask);  % BHTE subdomain

bhte_medium.alpha_power = medium.alpha_power * ones(bhte_Nx, bhte_Ny);  % medium.alpha_power -> scalar (2 in AxiSym)
bhte_medium.thermal_conductivity = 0.6045 * ones(bhte_Nx, bhte_Ny);     % [W/(m.K)]
bhte_medium.specific_heat = 4178 * ones(bhte_Nx, bhte_Ny);              % [J/(kg.K)]

if freefield ~= true

    % Air bottom domain
    air1_dxf = join([dxf_domains_dir 'Air1.dxf']);
    air1_bm = dxf2binary_mask(air1_dxf, kgrid.Nx, kgrid.Ny, kgrid.dx, kgrid.dy, 0);
    mirored_air1_bm = [air1_bm(:, end:-1:2) air1_bm];
    mirored_air1_bm = mirored_air1_bm(bhte_domain_x_mask, bhte_domain_y_mask);  % BHTE subdomain
    bhte_medium.thermal_conductivity(mirored_air1_bm) = 0.027381779;            % [W/(m.K)]  ok
    bhte_medium.specific_heat(mirored_air1_bm) = 1003.666667;                   % [J/(kg.K)] ok

    % Air upper domain
    air2_dxf = join([dxf_domains_dir 'Air2.dxf']);
    air2_bm = dxf2binary_mask(air2_dxf, kgrid.Nx, kgrid.Ny, kgrid.dx, kgrid.dy, 0);
    mirored_air2_bm = [air2_bm(:, end:-1:2) air2_bm];
    mirored_air2_bm = mirored_air2_bm(bhte_domain_x_mask, bhte_domain_y_mask);  % BHTE subdomain
    bhte_medium.thermal_conductivity(mirored_air2_bm) = 0.027381779;            % [W/(m.K)]  ok
    bhte_medium.specific_heat(mirored_air2_bm) = 1003.666667;                   % [J/(kg.K)] ok

    % Agarose coupling gel (Yang2017 / Culjat2010)
    agarose_gel_dxf = join([dxf_domains_dir 'AgaroseGel.dxf']);
    agarose_gel_bm = dxf2binary_mask(agarose_gel_dxf, kgrid.Nx, kgrid.Ny, kgrid.dx, kgrid.dy, 0);
    mirored_agarose_gel_bm = [agarose_gel_bm(:, end:-1:2) agarose_gel_bm];
    mirored_agarose_gel_bm = mirored_agarose_gel_bm(bhte_domain_x_mask, bhte_domain_y_mask);  % BHTE subdomain
    bhte_medium.thermal_conductivity(mirored_agarose_gel_bm) = 0.6045;          % [W/(m.K)]  ok
    bhte_medium.specific_heat(mirored_agarose_gel_bm) = 4178;                   % [J/(kg.K)] ok

    % Petri dish coverglass
    coverglass_dxf = join([dxf_domains_dir 'Coverglass.dxf']);
    coverglass_bm = dxf2binary_mask(coverglass_dxf, kgrid.Nx, kgrid.Ny, kgrid.dx, kgrid.dy, axial_tx_offset);
    mirored_coverglass_bm = [coverglass_bm(:, end:-1:2) coverglass_bm];
    mirored_coverglass_bm = mirored_coverglass_bm(bhte_domain_x_mask, bhte_domain_y_mask);  % BHTE subdomain
    bhte_medium.thermal_conductivity(mirored_coverglass_bm) = 1.2;              % [W/(m.K)]  ok
    bhte_medium.specific_heat(mirored_coverglass_bm) = 900;                     % [J/(kg.K)] ok

    % PDMS
    pdms_dxf = join([dxf_domains_dir 'Brainslice.dxf']);
    pdms_bm = dxf2binary_mask(pdms_dxf, kgrid.Nx, kgrid.Ny, kgrid.dx, kgrid.dy, axial_tx_offset);
    mirored_pdms_bm = [pdms_bm(:, end:-1:2) pdms_bm];
    mirored_pdms_bm = mirored_pdms_bm(bhte_domain_x_mask, bhte_domain_y_mask);  % BHTE subdomain
    bhte_medium.thermal_conductivity(mirored_pdms_bm) = 0.15;                   % [W/(m.K)]  ok
    bhte_medium.specific_heat(mirored_pdms_bm) = 1460;                          % [J/(kg.K)] ok
end

bhte_source.T0 = 37;
alpha_np = db2neper(bhte_medium.alpha_coeff, bhte_medium.alpha_power) .* (2 * pi * tx_freq).^bhte_medium.alpha_power;

bhte_input_args = {
    'PlotSim', true, ... % real time sim plot
    'PlotScale', [37, 50], ... % Huge color scale to account for the huge input power
};

bhte_ppw = max(bhte_medium.sound_speed, [], "all") / (tx_freq * dx); % points per wavelength
bhte_cfl = 40;                          % cfl number - OPTIM: 40 ???
bhte_ppp = ceil(bhte_ppw / bhte_cfl);   % points per period
bhte_T = 1 / tx_freq;                   % period [s]
bhte_dt = bhte_T / bhte_ppp;            % time step [s]

% Run kWaveDiffusion (no continuous recording to prevent huge RAM usage)
bhte_sensor.mask = zeros(bhte_Nx, bhte_Ny);
bhte_sensor_width = 5e-3; % m
bhte_sensor_y_mask = abs(bhte_kgrid.y_vec) < bhte_sensor_width / 2;
[mm, water_pdms_interface_idx] = min(abs(cells_z_vec(bhte_domain_x_mask)));
bhte_sensor.mask(water_pdms_interface_idx:water_pdms_interface_idx+1, bhte_sensor_y_mask) = 1;
kdiff = kWaveDiffusion(bhte_kgrid, bhte_medium, bhte_source, bhte_sensor, bhte_input_args{:});

% FUS on
kdiff.Q = alpha_np .* p_mag_mirrored(bhte_domain_x_mask, bhte_domain_y_mask).^2 ./ (bhte_medium.sound_speed .* bhte_medium.density);
kdiff.takeTimeStep(round(FUS_on_time / bhte_dt), bhte_dt);
bhte_result.dT0 = kdiff.T - bhte_source.T0;

% FUS off
kdiff.Q = 0;
kdiff.takeTimeStep(round(FUS_off_time / bhte_dt), bhte_dt);
bhte_result.dT1 = kdiff.T - bhte_source.T0;

% %% Cell Temperature elevation retreival
% 
% temperature_elevation_map = reshape(bhte_result.dT0, bhte_Nx, bhte_Ny);
% 
% % Temperature elevation map preview
% figure
% imagesc(x_vec_mirrored(bhte_domain_y_mask), cells_z_vec(bhte_domain_x_mask), temperature_elevation_map);
% yline(0, 'r', LineWidth=2)
% xlabel('x (radial) position [m]');
% ylabel('z (axial) position [m]');
% title(join(['Temperature elevation map (acoustic domain subset)' newline sprintf('%0.0e', tx_freq) 'Hz | ' sprintf('%0.0e', tx_emitted_ac_power) 'W | Tx z offset: ' sprintf('%0.0e', axial_tx_offset) ')'], ''));
% axis image

% %% ---- Cell temperature temporal profiles ----
% 
% bhte_data = reshape(kdiff.sensor_data - bhte_source.T0, 2, sum(bhte_sensor_y_mask), length(kdiff.sensor_data));
% bhte_result.bhte_t_vec = [0:1:kdiff.time_steps_taken-1] .* bhte_dt;
% 
% % -3dB focal spot region
% p_mag_peak = max(p_mag_mirrored(cell_location_index, bhte_domain_y_mask));
% p3dB = 10^((20*log10(p_mag_peak) - 3) / 20);
% p3dB_mask = p_mag_mirrored(cell_location_index, bhte_domain_y_mask) > p3dB;
% 
% % Water & PDMS temp elevation at the interface
% bhte_result.water_side_temp = squeeze(bhte_data(1, :, :));
% bhte_result.pdms_side_temp = squeeze(bhte_data(2, :, :));
% 
% % Peak pressure loc temp retreival trial
% % [mm, ppeak_idx] = max(p_mag_mirrored(cell_location_index, bhte_sensor_y_mask));
% 
% bhte_result.water_temp_peak = max(bhte_result.water_side_temp, [], 1);
% bhte_result.pdms_temp_peak = max(bhte_result.pdms_side_temp, [], 1);
% bhte_result.water_temp_fspotavg = mean(bhte_result.water_side_temp(p3dB_mask(bhte_sensor_y_mask), :), 1);
% bhte_result.pdms_temp_fspotavg = mean(bhte_result.pdms_side_temp(p3dB_mask(bhte_sensor_y_mask), :), 1);
% bhte_result.sensor_x_vec = x_vec_mirrored(bhte_sensor_y_mask);
% bhte_result.sensor_z_vec = cells_z_vec(bhte_domain_x_mask);
% 
% disp(sprintf('Max cell temperature elevation: %.3f°C', max(bhte_result.water_side_temp, [], 'all')));
% 
% % Temperature elevation profile preview
% figure
% hold on
% plot(x_vec_mirrored(bhte_sensor_y_mask), max(bhte_result.water_side_temp, [], 2), 'LineWidth', 1.5, 'DisplayName', 'Water T0-1 max')
% plot(x_vec_mirrored(bhte_sensor_y_mask), max(bhte_result.pdms_side_temp, [], 2), 'LineWidth', 1.5, 'DisplayName', 'PDMS T0-1 max')
% title(join(['Water & PDMS temporal peak temperature elevation' newline sprintf('%0.0e', tx_freq) 'Hz | ' sprintf('%0.0e', tx_emitted_ac_power) 'W | Tx z offset: ' sprintf('%0.0e', axial_tx_offset) ')'], ''));
% legend
% 
% % PDMS temperature elevation time profile
% figure
% hold on
% plot(bhte_result.bhte_t_vec, bhte_result.pdms_temp_peak, ':', 'LineWidth', 2, 'DisplayName', 'PDMS peak')
% plot(bhte_result.bhte_t_vec, bhte_result.pdms_temp_fspotavg, ':', 'LineWidth', 2, 'DisplayName', 'PDMS -3dB fspot avg')
% title(join(['PDMS temperature elevation profile' newline sprintf('%0.0e', tx_freq) 'Hz | ' sprintf('%0.0e', tx_emitted_ac_power) 'W | Tx z offset: ' sprintf('%0.0e', axial_tx_offset) ')'], ''));
% legend('Location', 'northwest')
% 
% % Water temperature elevation time profile
% figure
% hold on
% plot(bhte_result.bhte_t_vec, bhte_result.water_temp_peak, '-', 'LineWidth', 2, 'DisplayName', 'Water peak')
% plot(bhte_result.bhte_t_vec, bhte_result.water_temp_fspotavg, '-', 'LineWidth', 2, 'DisplayName', 'Water -3dB fspot avg')
% title(join(['Cells temperature elevation profile' newline sprintf('%0.0e', tx_freq) 'Hz | ' sprintf('%0.0e', tx_emitted_ac_power) 'W | Tx z offset: ' sprintf('%0.0e', axial_tx_offset) ')'], ''));
% legend('Location', 'northwest')