clc; clear; close all;

%% ================== SYSTEM PARAMETERS ==================
% Channel and electrode design
channel_length = 500e-6;      % Length (m)
num_electrodes = 10;          % Number of electrode pairs (IDEs)
spacing = channel_length / num_electrodes;

% Physical properties
epsilon_0 = 8.854e-12;        % Vacuum permittivity (F/m)
epsilon_r = 80;               % Relative permittivity (water)
epsilon_m = epsilon_0 * epsilon_r;
r_cell = 10e-6;               % Cell radius (m)
sigma_p = 0.5;                % Cell conductivity (S/m)
sigma_m = 1e-3;               % Medium conductivity (S/m)
rho_cell = 1050;              % Cell density (kg/m³)
eta = 0.001;                  % Dynamic viscosity (Pa·s)
u_flow = 1e-4;                % Flow velocity (m/s)

% Electric field parameters
E0_DEP = 1.1e5;              % DEP field strength (V/m)
E_threshold = 0.8 * E0_DEP;   % Trapping threshold (V/m)
electroporation_thresh = 2.5e5; % V/m

% Electrofusion parameters
acFrequency = 1e6;            % 1 MHz
acVoltage = 10;               % Vpp
dcFieldStrength = 0.7e3;      % V/cm
pulseDuration = 100e-6;       % 100 µs
numPulses = 10;
fusionSuccessRate = 0.8;       % Probability of successful fusion

% Color codes for cell states
colors = struct('untrapped', [0.7 0.7 0.7], ...
                'trapped', [0 0.5 1], ...
                'paired', [1 0.6 0], ...
                'fused', [0 0.8 0], ...
                'failed', [1 0 0]);

%% ================== DEP SYSTEM SETUP ==================
% Field calculation
x = linspace(0, channel_length, 1000);
E_DEP = E0_DEP * cos(2 * pi * x / spacing);  
dEdx = gradient(E_DEP.^2, x);
Re_CM = real((sigma_p - sigma_m) / (sigma_p + 2*sigma_m));
F_DEP = 2 * pi * epsilon_m * r_cell^3 * Re_CM .* dEdx;

%% ================== CELL INITIALIZATION ==================
num_cells = 20;  % Reduced for visualization
cells = struct('id', {}, 'status', {}, 'partner_id', {}, 'x_pos', {}, 'y_pos', {}, 'velocity', {});

% Random initial positions along channel
for i = 1:num_cells
    cells(i).id = i;
    cells(i).status = 'untrapped';
    cells(i).partner_id = NaN;
    cells(i).x_pos = rand() * channel_length;
    cells(i).y_pos = 1 + 0.2*randn();  % Slight vertical variation
    cells(i).velocity = [0; 0];
end

%% ================== SIMULATION LOOP ==================
dt = 1e-5;                   % Time step (s)
T_total = 0.1;               % Total simulation time (s)
sim_steps = round(T_total/dt);
fusion_phase_start = 0.05;   % Start fusion after initial trapping

% Preallocate tracking
cell_trajectories = zeros(num_cells, sim_steps, 2);

for t = 1:sim_steps
    current_time = t*dt;
    
    %% DEP Trapping Phase
    if current_time < fusion_phase_start
        for i = 1:num_cells
            % Find nearest field point
            [~, idx] = min(abs(x - cells(i).x_pos));
            
            % Calculate forces
            F_drag = -6 * pi * eta * r_cell * cells(i).velocity;
            F_total = [F_DEP(idx); 0] + F_drag;
            
            % Update position and velocity
            cells(i).velocity = cells(i).velocity + (F_total / ((4/3)*pi*r_cell^3*rho_cell)) * dt;
            cells(i).x_pos = cells(i).x_pos + cells(i).velocity(1) * dt;
            cells(i).y_pos = cells(i).y_pos + cells(i).velocity(2) * dt;
            
            % Periodic boundary
            cells(i).x_pos = mod(cells(i).x_pos, channel_length);
            
            % Check trapping status
            if abs(E_DEP(idx)) > E_threshold
                cells(i).status = 'trapped';
                % Freeze trapped cells near electrode peaks
                cells(i).velocity = [0; 0];
                cells(i).x_pos = round(cells(i).x_pos/spacing)*spacing;
            end
            
            cell_trajectories(i,t,:) = [cells(i).x_pos, cells(i).y_pos];
        end
        
        % Pair trapped cells (every 2nd electrode)
        if current_time > 0.02
            for i = 1:2:num_cells-1
                if strcmp(cells(i).status, 'trapped') && strcmp(cells(i+1).status, 'trapped')
                    cells(i).status = 'paired';
                    cells(i+1).status = 'paired';
                    cells(i).partner_id = cells(i+1).id;
                    cells(i+1).partner_id = cells(i).id;
                    % Move paired cells together
                    mean_x = mean([cells(i).x_pos, cells(i+1).x_pos]);
                    cells(i).x_pos = mean_x - spacing/4;
                    cells(i+1).x_pos = mean_x + spacing/4;
                end
            end
        end
        
    %% Electrofusion Phase
    else
        % Apply AC field for alignment
        if current_time < fusion_phase_start + 0.01
            for i = 1:num_cells
                if strcmp(cells(i).status, 'paired')
                    % Additional alignment force
                    partner_idx = cells(i).partner_id;
                    dx = cells(partner_idx).x_pos - cells(i).x_pos;
                    alignment_force = 1e-12 * sign(dx - spacing/10);
                    cells(i).x_pos = cells(i).x_pos + alignment_force * dt;
                end
            end
        
        % Apply DC pulses
        elseif current_time > fusion_phase_start + 0.02
            for i = 1:2:num_cells-1
                if strcmp(cells(i).status, 'paired') && rand() <= fusionSuccessRate
                    cells(i).status = 'fused';
                    cells(i+1).status = 'fused';
                    % Fused cells become one (visualization only)
                    cells(i).x_pos = mean([cells(i).x_pos, cells(i+1).x_pos]);
                    cells(i+1).x_pos = cells(i).x_pos;
                elseif strcmp(cells(i).status, 'paired')
                    cells(i).status = 'failed';
                    cells(i+1).status = 'failed';
                end
            end
        end
    end
    
    % Visualize every 100 steps
    if mod(t,100) == 0
        visualizeSystem(cells, channel_length, colors, current_time);
    end
end

%% ================== ANALYSIS & VISUALIZATION ==================
% Plot field profiles
figure;
subplot(2,1,1);
plot(x*1e6, E_DEP/1e5, 'b', 'LineWidth', 2);
xlabel('Position (μm)'); ylabel('DEP Field (×10^5 V/m)');
title('DEP Field Profile'); grid on;

subplot(2,1,2);
plot(x*1e6, F_DEP/1e-12, 'r', 'LineWidth', 2);
xlabel('Position (μm)'); ylabel('DEP Force (pN)');
title('DEP Force Profile'); grid on;

% Plot trajectories
figure;
hold on;
for i = 1:num_cells
    plot(squeeze(cell_trajectories(i,:,1))*1e6, squeeze(cell_trajectories(i,:,2))*1e6, ...
        'Color', colors.(cells(i).status), 'LineWidth', 1);
end
xlabel('X Position (μm)'); ylabel('Y Position (μm)');
title('Cell Trajectories'); grid on;

% Fusion statistics
fused_pairs = sum(strcmp({cells.status}, 'fused'))/2;
fprintf('\n=== System Performance ===\n');
fprintf('Successfully fused pairs: %d/%d (%.1f%%)\n', ...
    fused_pairs, floor(num_cells/2), 100*fused_pairs/floor(num_cells/2));
fprintf('Average trapping time: %.2f ms\n', fusion_phase_start*1e3);

%% ================== STATISTICAL METRICS ==================
% Calculate trapping efficiency
trapped_cells = find(strcmp({cells.status}, 'trapped') | ...
                    strcmp({cells.status}, 'paired') | ...
                    strcmp({cells.status}, 'fused') | ...
                    strcmp({cells.status}, 'failed'));
trapping_efficiency = length(trapped_cells)/num_cells * 100;

% Calculate cell velocities (x-component only)
final_velocities = zeros(1, num_cells);
for i = 1:num_cells
    final_velocities(i) = cells(i).velocity(1);
end
mean_velocity = mean(final_velocities);
std_velocity = std(final_velocities);

% Calculate field exposure metrics
max_E_exposure = zeros(1, num_cells);
high_field_time = zeros(1, num_cells);
for i = 1:num_cells
    % Get all positions for this cell
    cell_x_positions = squeeze(cell_trajectories(i,:,1));
    
    % Interpolate E-field at all positions
    E_at_positions = interp1(x, abs(E_DEP), mod(cell_x_positions, channel_length));
    
    max_E_exposure(i) = max(E_at_positions);
    high_field_time(i) = sum(E_at_positions > E_threshold) * dt;
end

% Display the requested metrics
fprintf('\n===== DEP Device Performance Metrics =====\n');
fprintf('Trapping Efficiency: %.1f%%\n', trapping_efficiency);
fprintf('Mean Cell Velocity: %.2f ± %.2f μm/s\n', mean_velocity*1e6, std_velocity*1e6);
fprintf('Max E-field Exposure: %.1f - %.1f ×10^5 V/m\n', ...
    min(max_E_exposure)/1e5, max(max_E_exposure)/1e5);
fprintf('Average Time in High Field: %.2f ms\n', mean(high_field_time)*1e3);

%% Helper Function
function visualizeSystem(cells, channel_length, colors, current_time)
    clf;
    hold on;
    
    % Draw electrodes
    electrode_pos = linspace(0, channel_length, 10);
    for pos = electrode_pos
        rectangle('Position', [pos*1e6-5, 0.7*1e6, 10, 10], ...
                 'FaceColor', [0.3 0.3 0.3], 'EdgeColor', 'none');
    end
    
    % Plot cells
    for i = 1:length(cells)
        scatter(cells(i).x_pos*1e6, cells(i).y_pos*1e6, 200, ...
               'MarkerFaceColor', colors.(cells(i).status), ...
               'MarkerEdgeColor', 'k');
        text(cells(i).x_pos*1e6, cells(i).y_pos*1e6+3, num2str(cells(i).id), ...
            'HorizontalAlignment', 'center', 'FontSize', 8);
        
        % Draw connection lines for pairs
        if ~isnan(cells(i).partner_id)
            partner = cells(i).partner_id;
            plot([cells(i).x_pos, cells(partner).x_pos]*1e6, ...
                 [cells(i).y_pos, cells(partner).y_pos]*1e6, ...
                 'Color', [0.5 0.5 0.5], 'LineStyle', '--');
        end
    end
    
    xlim([0 channel_length*1e6]);
    ylim([0.5*1e6 1.5*1e6]);
    xlabel('Position (μm)'); 
    title(sprintf('DEP-EF System | Time: %.1f ms', current_time*1e3));
    drawnow;
end
