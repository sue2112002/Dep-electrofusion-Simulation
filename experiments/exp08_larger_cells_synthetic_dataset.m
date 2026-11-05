clc; clear; close all;

%% ================== SYSTEM PARAMETERS ==================
% Channel and electrode design
channel_length = 200e-6;      % Length (m), matches dataset range (~200 µm)
num_electrodes = 10;          % Number of electrode pairs
spacing = channel_length / num_electrodes;

% Physical properties
epsilon_0 = 8.854e-12;        % Vacuum permittivity (F/m)
epsilon_r = 80;               % Relative permittivity (water)
epsilon_m = epsilon_0 * epsilon_r;
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
E_dep_base = E0_DEP * cos(2 * pi * x / spacing);  % Base cosine field

% Enhance field at second (20e-6 m) and fourth (60e-6 m) electrode pairs
E_enhancement = zeros(size(x));
electrode_positions = [1*spacing, 3*spacing]; % Positions of 2nd and 4th electrodes
enhancement_factors = [0.25, 0.25]; % 25% increase at both
sigma = spacing / 10; % Width of enhancement

for i = 1:length(electrode_positions)
    pos = electrode_positions(i);
    factor = enhancement_factors(i);
    E_enhancement = E_enhancement + factor * exp(-((x - pos).^2) / (2 * sigma^2));
end
E_dep = E_dep_base .* (1 + E_enhancement); % Apply enhancement
dEdx = gradient(E_dep.^2, x);

%% ================== CELL INITIALIZATION ==================
% Generate synthetic dataset for 400 cells (200 PLT, 200 RBC)
num_cells = 400;
data = table((1:num_cells)', ... % ID
             [repmat("PLT", 200, 1); repmat("RBC", 200, 1)], ... % Type
             rand(num_cells, 1) * 200, ... % y (0–200 µm)
             [1.5 + rand(200, 1) * 1.5; 3 + rand(200, 1) * 4.5], ... % radius (PLT: 1.5–3 µm, RBC: 3–7.5 µm)
             [0.1 + rand(200, 1) * 0.2; 0 + rand(200, 1) * 1.5], ... % sigma_cyt
             50 + rand(num_cells, 1) * 30, ... % epsilon_cyt
             'VariableNames', {'ID', 'Type', 'y', 'radius', 'sigma_cyt', 'epsilon_cyt'});

% Initialize cells
cells = struct('id', {}, 'status', {}, 'partner_id', {}, 'x_pos', {}, 'y_pos', {}, 'velocity', {}, ...
               'radius', {}, 'sigma_p', {}, 'epsilon_p', {}, 'type', {});

for i = 1:num_cells
    cells(i).id = data.ID(i);
    cells(i).type = string(data.Type(i));
    cells(i).status = 'untrapped';
    cells(i).partner_id = NaN;
    cells(i).x_pos = data.y(i) * 1e-6; % Convert y from µm to m
    cells(i).y_pos = 1 + 0.2*randn();
    cells(i).velocity = [0; 0];
    cells(i).radius = data.radius(i) * 1e-6; % Convert radius from µm to m
    cells(i).sigma_p = data.sigma_cyt(i);
    cells(i).epsilon_p = data.epsilon_cyt(i) * epsilon_0;
end

% Compute representative F_dep for plotting
mean_radius = mean([cells.radius]); % Updated mean radius (~3.75e-6 m)
mean_sigma_p = mean([cells.sigma_p]);
Re_CM_mean = real((mean_sigma_p - sigma_m) / (mean_sigma_p + 2*sigma_m));
F_dep = 2 * pi * epsilon_m * mean_radius^3 * Re_CM_mean * dEdx; % Updated DEP force

% Precompute Re_CM for each cell
Re_CM = zeros(1, num_cells);
for i = 1:num_cells
    Re_CM(i) = real((cells(i).sigma_p - sigma_m) / (cells(i).sigma_p + 2*sigma_m));
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
            F_drag = -6 * pi * eta * cells(i).radius * cells(i).velocity;
            F_dep_i = 2 * pi * epsilon_m * cells(i).radius^3 * Re_CM(i) * dEdx(idx);
            F_total = [F_dep_i; 0] + F_drag;
            
            % Update position and velocity
            cells(i).velocity = cells(i).velocity + (F_total / ((4/3)*pi*cells(i).radius^3*rho_cell)) * dt;
            cells(i).x_pos = cells(i).x_pos + cells(i).velocity(1) * dt;
            cells(i).y_pos = cells(i).y_pos + cells(i).velocity(2) * dt;
            
            % Periodic boundary
            cells(i).x_pos = mod(cells(i).x_pos, channel_length);
            
            % Check trapping status
            if abs(E_dep(idx)) > E_threshold
                cells(i).status = 'trapped';
                cells(i).velocity = [0; 0];
                cells(i).x_pos = round(cells(i).x_pos/spacing)*spacing;
            end
            
            cell_trajectories(i,t,:) = [cells(i).x_pos, cells(i).y_pos];
        end
        
        % Pair trapped cells
        if current_time > 0.02
            trapped_indices = find(strcmp({cells.status}, 'trapped'));
            for i = 1:2:length(trapped_indices)-1
                idx1 = trapped_indices(i);
                idx2 = trapped_indices(i+1);
                
                if strcmp(cells(idx1).type, 'PLT') && strcmp(cells(idx2).type, 'RBC') || ...
                   strcmp(cells(idx1).type, 'RBC') && strcmp(cells(idx2).type, 'PLT') || ...
                   strcmp(cells(idx1).type, cells(idx2).type)
                    cells(idx1).status = 'paired';
                    cells(idx2).status = 'paired';
                    cells(idx1).partner_id = cells(idx2).id;
                    cells(idx2).partner_id = cells(idx1).id;
                    mean_x = mean([cells(idx1).x_pos, cells(idx2).x_pos]);
                    cells(idx1).x_pos = mean_x - spacing/4;
                    cells(idx2).x_pos = mean_x + spacing/4;
                end
            end
        end
        
    %% Electrofusion Phase
    else
        % Apply AC field for alignment
        if current_time < fusion_phase_start + 0.01
            for i = 1:num_cells
                if strcmp(cells(i).status, 'paired')
                    partner_idx = find([cells.id] == cells(i).partner_id);
                    dx = cells(partner_idx).x_pos - cells(i).x_pos;
                    alignment_force = 1e-12 * sign(dx - spacing/10);
                    cells(i).x_pos = cells(i).x_pos + alignment_force * dt;
                end
            end
        
        % Apply DC pulses
        elseif current_time > fusion_phase_start + 0.02
            for i = 1:num_cells
                if strcmp(cells(i).status, 'paired') && rand() <= fusionSuccessRate
                    partner_idx = find([cells.id] == cells(i).partner_id);
                    if ~isempty(partner_idx)
                        cells(i).status = 'fused';
                        cells(partner_idx).status = 'fused';
                        cells(i).x_pos = mean([cells(i).x_pos, cells(partner_idx).x_pos]);
                        cells(partner_idx).x_pos = cells(i).x_pos;
                    end
                elseif strcmp(cells(i).status, 'paired')
                    partner_idx = find([cells.id] == cells(i).partner_id);
                    if ~isempty(partner_idx)
                        cells(i).status = 'failed';
                        cells(partner_idx).status = 'failed';
                    end
                end
            end
        end
    end
    
    % Visualize every 1000 steps
    if mod(t,1000) == 0
        visualizeSystem(cells, channel_length, colors, current_time);
    end
end

%% ================== ANALYSIS & VISUALIZATION ==================
% Plot field profiles
figure;
subplot(2,1,1);
plot(x*1e6, E_dep/1e5, 'b', 'LineWidth', 2);
xlabel('Position (μm)'); ylabel('DEP Field (×10^5 V/m)');
title('DEP Field Profile'); grid on;

subplot(2,1,2);
plot(x*1e6, F_dep/1e-12, 'r', 'LineWidth', 2);
xlabel('Position (μm)'); ylabel('DEP Force (pN)');
title('DEP Force Profile'); grid on;

% Plot trajectories for 10 cells (5 PLT, 5 RBC)
figure;
hold on;
plt_indices = find(strcmp({cells.type}, 'PLT'), 5); % First 5 PLT
rbc_indices = find(strcmp({cells.type}, 'RBC'), 5); % First 5 RBC
for i = plt_indices
    plot(squeeze(cell_trajectories(i,:,1))*1e6, squeeze(cell_trajectories(i,:,2))*1e6, ...
         'Color', [0.2 0.6 1], 'LineWidth', 1.5, 'DisplayName', sprintf('PLT %d', cells(i).id));
end
for i = rbc_indices
    plot(squeeze(cell_trajectories(i,:,1))*1e6, squeeze(cell_trajectories(i,:,2))*1e6, ...
         'Color', [1 0.2 0.2], 'LineWidth', 1.5, 'DisplayName', sprintf('RBC %d', cells(i).id));
end
xlabel('X Position (μm)'); ylabel('Y Position (μm)');
title('Trajectories of Enlarged Cells (5 PLT, 5 RBC)');
grid on;
legend('show', 'Location', 'best');
xlim([0 channel_length*1e6]);
ylim([0.5*1e6 1.5*1e6]);

% Fusion statistics
fused_pairs = sum(strcmp({cells.status}, 'fused'))/2;
fprintf('\n=== System Performance ===\n');
fprintf('Successfully fused pairs: %d/%d (%.1f%%)\n', ...
    fused_pairs, floor(num_cells/2), 100*fused_pairs/floor(num_cells/2));
fprintf('Average trapping time: %.2f ms\n', fusion_phase_start*1e3);

% Additional statistics
plt_rbc_fused = 0;
same_type_fused = 0;
for i = 1:num_cells
    if strcmp(cells(i).status, 'fused')
        partner_idx = find([cells.id] == cells(i).partner_id);
        if ~isempty(partner_idx)
            if strcmp(cells(i).type, 'PLT') && strcmp(cells(partner_idx).type, 'RBC') || ...
               strcmp(cells(i).type, 'RBC') && strcmp(cells(partner_idx).type, 'PLT')
                plt_rbc_fused = plt_rbc_fused + 0.5;
            elseif strcmp(cells(i).type, cells(partner_idx).type)
                same_type_fused = same_type_fused + 0.5;
            end
        end
    end
end
fprintf('PLT-RBC fused pairs: %d (%.1f%%)\n', plt_rbc_fused, 100*plt_rbc_fused/fused_pairs);
fprintf('Same-type fused pairs: %d (%.1f%%)\n', same_type_fused, 100*same_type_fused/fused_pairs);

%% Helper Function
function visualizeSystem(cells, channel_length, colors, current_time)
    clf;
    hold on;
    
    % Draw parallel electrodes
    electrode_pos = linspace(0, channel_length, 10);
    electrode_width = 10;
    electrode_height = 10;
    for pos = electrode_pos
        rectangle('Position', [pos*1e6 - electrode_width/2, 1.3*1e6, electrode_width, electrode_height], ...
                 'Curvature', [1 0], 'FaceColor', [0.3 0.3 0.3], 'EdgeColor', 'none');
        rectangle('Position', [pos*1e6 - electrode_width/2, 0.7*1e6, electrode_width, electrode_height], ...
                 'Curvature', [1 0], 'FaceColor', [0.3 0.3 0.3], 'EdgeColor', 'none');
    end
    
    % Plot cells (up to 50)
    max_cells_plot = min(50, length(cells));
    for i = 1:max_cells_plot
        scatter(cells(i).x_pos*1e6, cells(i).y_pos*1e6, 60, ... % Increased marker size
               'MarkerFaceColor', colors.(cells(i).status), ...
               'MarkerEdgeColor', 'k');
        text(cells(i).x_pos*1e6, cells(i).y_pos*1e6+3, sprintf('%d (%s)', cells(i).id, cells(i).type), ...
            'HorizontalAlignment', 'center', 'FontSize', 6);
        
        % Draw connection lines for pairs
        if ~isnan(cells(i).partner_id)
            partner_idx = find([cells.id] == cells(i).partner_id);
            if ~isempty(partner_idx)
                plot([cells(i).x_pos, cells(partner_idx).x_pos]*1e6, ...
                     [cells(i).y_pos, cells(partner_idx).y_pos]*1e6, ...
                     'Color', [0.5 0.5 0.5], 'LineStyle', '--');
            end
        end
    end
    
    xlim([0 channel_length*1e6]);
    ylim([0.5*1e6 1.5*1e6]);
    xlabel('Position (μm)'); 
    title(sprintf('DEP-EF System with Enlarged Cells | Time: %.1f ms (Showing %d/%d Cells)', ...
        current_time*1e3, max_cells_plot, length(cells)));
    grid on;
    drawnow;
end