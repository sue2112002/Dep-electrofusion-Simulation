clc; clear; close all;

%% ------------------ PARAMETERS ------------------
% Channel and electrode design
channel_length = 500e-6;  % Length (m)
num_electrodes = 10;      % Number of electrode pairs (IDEs)
spacing = channel_length / num_electrodes;  % IDE periodicity (m)

% Physical properties
epsilon_0 = 8.854e-12;     % Vacuum permittivity (F/m)
epsilon_r = 80;            % Relative permittivity (water)
epsilon_m = epsilon_0 * epsilon_r;
r_cell = 10e-6;            % Cell radius (m)
sigma_p = 0.5;             % Cell conductivity (S/m) - adjusted for HSPCs
sigma_m = 1e-3;            % Medium conductivity (S/m)
rho_cell = 1050;           % Cell density (kg/m³)
eta = 0.001;               % Dynamic viscosity (Pa·s)
u_flow = 1e-4;             % Flow velocity (m/s)

% Electric field (IDE model)
E0 = 1.1e5;                % Max E-field (V/m)
electroporation_thresh = 2.5e5;  % Threshold (V/m)
E_threshold = 0.8 * E0;    % Threshold for trapping (V/m)

%% ------------------ FIELD CALCULATION ------------------
x = linspace(0, channel_length, 1000);  % Spatial domain

% Interdigitated electrode (IDE) field model (simplified cosine)
E = E0 * cos(2 * pi * x / spacing);  
dEdx = gradient(E.^2, x);  % Gradient of E² for DEP force

% Clausius-Mossotti factor (updated for HSPCs)
Re_CM = real((sigma_p - sigma_m) / (sigma_p + 2*sigma_m));

% DEP force
F_DEP = 2 * pi * epsilon_m * r_cell^3 * Re_CM .* dEdx;

% Joule heating (simplified)
deltaT = (sigma_m * E.^2) ./ (1000 * 4186);  % ΔT in K (water: ρ=1000, cp=4186)

%% ------------------ MULTI-CELL TRACKING ------------------
dt = 1e-5;                 % Time step (s)
T_end = 0.05;              % Total simulation time (s)
time = 0:dt:T_end;
num_cells = 1000;          % Number of cells to simulate

% Preallocate
x_cells = zeros(num_cells, length(time));
v_cells = zeros(num_cells, length(time));
x_cells(:, 1) = rand(num_cells, 1) * channel_length;  % Random initial positions

mass = (4/3) * pi * r_cell^3 * rho_cell;
k_drag = 6 * pi * eta * r_cell;

for t = 2:length(time)
    for n = 1:num_cells
        % Find nearest spatial index for field lookup
        [~, idx] = min(abs(x - x_cells(n, t-1)));
        
        % DEP force + flow advection
        F_now = F_DEP(idx);
        a = (F_now - k_drag * v_cells(n, t-1)) / mass;
        v_cells(n, t) = v_cells(n, t-1) + a * dt;
        x_cells(n, t) = x_cells(n, t-1) + (v_cells(n, t) + u_flow) * dt;
        
        % Periodic boundary conditions
        x_cells(n, t) = mod(x_cells(n, t), channel_length);
    end
end

%% ------------------ STATISTICAL METRICS ------------------
% Preallocate
max_E_exposure = zeros(num_cells, 1);
high_field_time = zeros(num_cells, 1);
E_abs = abs(E);

% Vectorized field exposure calculation
for n = 1:num_cells
    % Interpolate E-field at cell positions
    x_pos = mod(x_cells(n, :), channel_length);
    E_at_pos = interp1(x, E_abs, x_pos, 'nearest');
    
    max_E_exposure(n) = max(E_at_pos);
    high_field_time(n) = sum(E_at_pos > E_threshold) * dt;
end

% Trapping efficiency (cells in high-field at end time)
trapped = max_E_exposure > E_threshold;
trapping_efficiency = mean(trapped) * 100;

% Velocity statistics
final_velocities = v_cells(:, end);
mean_velocity = mean(final_velocities);
std_velocity = std(final_velocities);

%% ------------------ VIABILITY CHECKS ------------------
if any(deltaT > 1)
    warning('⚠️ Joule heating exceeds 1 K. May affect cell viability.');
end

if any(E_abs > electroporation_thresh)
    warning('⚠️ Electric field exceeds electroporation threshold!');
end

%% ------------------ PLOTS ------------------
figure;
set(gcf, 'Position', [100, 100, 900, 700]);

% Field and force profiles
subplot(3,1,1);
plot(x * 1e6, E / 1e5, 'b', 'LineWidth', 2);
xlabel('Position (μm)'); ylabel('E (×10^5 V/m)');
title('Interdigitated Electrode Field Profile'); grid on;

subplot(3,1,2);
plot(x * 1e6, F_DEP / 1e-12, 'r', 'LineWidth', 2);
xlabel('Position (μm)'); ylabel('DEP Force (pN)');
title('DEP Force Profile'); grid on;

subplot(3,1,3);
plot(x * 1e6, deltaT, 'k', 'LineWidth', 2);
xlabel('Position (μm)'); ylabel('ΔT (K)');
title('Joule Heating Profile'); grid on;

% Trajectory plot (subset)
figure;
hold on;
for n = 1:min(20, num_cells)
    plot(time * 1e3, x_cells(n, :) * 1e6, 'LineWidth', 1);
end
xlabel('Time (ms)'); ylabel('Cell Position (μm)');
title('Cell Trajectories Under DEP + Flow');
grid on;

% Statistical plots
figure;
set(gcf, 'Position', [200, 200, 1000, 800]);

subplot(2,2,1);
histogram(max_E_exposure / 1e5, 20, 'FaceColor', [0.2 0.6 0.9]);
xlabel('Max E-field (×10^5 V/m)'); ylabel('Cell Count');
title('Field Exposure Distribution');

subplot(2,2,2);
histogram(final_velocities * 1e6, 20, 'FaceColor', [0.9 0.3 0.2]);
xlabel('Final Velocity (μm/s)'); ylabel('Cell Count');
title('Velocity Distribution');

subplot(2,2,3);
histogram(high_field_time * 1e3, 20, 'FaceColor', [0.3 0.8 0.4]);
xlabel('Time in High Field (ms)'); ylabel('Cell Count');
title('High-Field Residence Time');

subplot(2,2,4);
bar(trapping_efficiency, 'FaceColor', [0.7 0.2 0.7]);
ylabel('Efficiency (%)'); ylim([0 100]);
title(sprintf('Trapping Efficiency: %.1f%%', trapping_efficiency));
text(1, trapping_efficiency-10, sprintf('n = %d cells', num_cells), 'HorizontalAlignment', 'center');

%% ------------------ DISPLAY STATS ------------------
fprintf('\n===== DEP Device Performance Metrics =====\n');
fprintf('Trapping Efficiency: %.1f%%\n', trapping_efficiency);
fprintf('Mean Cell Velocity: %.2f ± %.2f μm/s\n', mean_velocity*1e6, std_velocity*1e6);
fprintf('Max E-field Exposure: %.1f - %.1f ×10^5 V/m\n', ...
    min(max_E_exposure)/1e5, max(max_E_exposure)/1e5);
fprintf('Average Time in High Field: %.2f ms\n', mean(high_field_time)*1e3);