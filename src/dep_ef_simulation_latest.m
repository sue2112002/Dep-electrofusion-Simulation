clc; clear; close all;

%% ================== SYSTEM PARAMETERS ==================
% Channel and electrode design
channel_length = 200e-6;      % Length (m), 200 µm
num_electrodes = 10;          % Number of electrode pairs
spacing = channel_length / num_electrodes; % 20e-6 m (20 µm)

% Physical properties
epsilon_0 = 8.854e-12;        % Vacuum permittivity (F/m)
epsilon_r = 80;               % Relative permittivity (water)
epsilon_m = epsilon_0 * epsilon_r;
sigma_m = 1e-3;               % Medium conductivity (S/m)
rho_cell = 1050;              % Cell density (kg/m³)
eta = 0.001;                  % Dynamic viscosity (Pa·s)
u_flow = 1e-4;                % Flow velocity (m/s)

% Electric field parameters
E0_DEP = 1.5e5;              % DEP field strength (V/m), increased for faster trapping
E_threshold = 0.8 * E0_DEP;   % Trapping threshold (V/m), 1.2e5 V/m
electroporation_thresh = 2.5e5; % V/m

% Electrofusion parameters
acFrequency = 1e6;            % 1 MHz
acVoltage = 10;               % Vpp
dcFieldStrength = 0.7e3;      % V/cm
pulseDuration = 50e-6;        % 50 µs, shorter for sharper pulses
numPulses = 15;               % Increased for higher fusion yield
fusionSuccessRate = 0.95;      % Increased to 95% for better fusion

% Color codes for cell states
colors = struct('untrapped', [0.7 0.7 0.7], ...
                'trapped', [0 0.5 1], ...
                'paired', [1 0.6 0], ...
                'fused', [0 0.8 0], ...
                'failed', [1 0 0]);

%% ================== DEP SYSTEM SETUP ==================
% Field calculation
x = linspace(0, channel_length, 1000);
E_DEP_base = E0_DEP * cos(2 * pi * x / spacing);  % Base cosine field

% Enhance field at second (20e-6 m) and fourth (60e-6 m) electrode pairs
E_enhancement = zeros(size(x));
electrode_positions = [1*spacing, 3*spacing]; % Positions of 2nd and 4th electrodes (20e-6 m, 60e-6 m)
enhancement_factors = [0.5, 0.5]; % 50% increase for stronger trapping
sigma = spacing / 10; % Width of enhancement (2e-6 m)

for i = 1:length(electrode_positions)
    pos = electrode_positions(i);
    factor = enhancement_factors(i);
    E_enhancement = E_enhancement + factor * exp(-((x - pos).^2) / (2 * sigma^2));
end
E_DEP = E_DEP_base .* (1 + E_enhancement); % Apply enhancement to base field
dEdx = gradient(E_DEP.^2, x);

%% ================== CELL INITIALIZATION ==================
% Load dataset
data = readtable('input_cells_dataset.csv');

% Select all cells
num_cells = height(data); % Simulate all cells
cells = struct('id', {}, 'status', {}, 'partner_id', {}, 'x_pos', {}, 'y_pos', {}, 'velocity', {}, ...
               'radius', {}, 'sigma_p', {}, 'epsilon_p', {}, 'type', {});

% Initialize cells from dataset with scaled positions
for i = 1:num_cells
    cells(i).id = data.ID(i);
    cells(i).type = string(data.Type(i)); % Store cell type (PLT or RBC)
    cells(i).status = 'untrapped';
    cells(i).partner_id = NaN;
    % Scale initial positions to channel length
    cells(i).x_pos = (data.y(i) / max(data.y)) * channel_length; % Normalize to 0–200 µm
    cells(i).y_pos = 1 + 0.2*randn();  % Slight vertical variation
    cells(i).velocity = [0; 0];
    cells(i).radius = data.radius(i) * 1e-6; % Convert radius from µm to m
    cells(i).sigma_p = data.sigma_cyt(i); % Cell conductivity (S/m)
    cells(i).epsilon_p = data.epsilon_cyt(i) * epsilon_0; % Convert relative to absolute permittivity
end

% Check PLT/RBC balance
plt_count = sum(strcmp({cells.type}, 'PLT'));
rbc_count = sum(strcmp({cells.type}, 'RBC'));
if abs(plt_count - rbc_count) / num_cells > 0.2
    warning('PLT/RBC ratio is unbalanced (PLT: %d, RBC: %d). Optimal pairing may be limited.', ...
            plt_count, rbc_count);
end

% Compute representative F_DEP for plotting (using mean radius and Re_CM)
mean_radius = mean([cells.radius]); % Mean radius of cells (m)
mean_sigma_p = mean([cells.sigma_p]); % Mean conductivity of cells (S/m)
Re_CM_mean = real((mean_sigma_p - sigma_m) / (mean_sigma_p + 2*sigma_m));
F_DEP = 2 * pi * epsilon_m * mean_radius^3 * Re_CM_mean * dEdx; % Representative DEP force

% Precompute Re_CM for each cell
Re_CM = zeros(1, num_cells);
for i = 1:num_cells
    Re_CM(i) = real((cells(i).sigma_p - sigma_m) / (cells(i).sigma_p + 2*sigma_m));
end

%% ================== SIMULATION LOOP ==================
dt = 1e-5;                   % Time step (s)
T_total = 0.1;               % Total simulation time (s)
sim_steps = round(T_total/dt);
fusion_phase_start = 0.08;   % Start fusion after extended trapping (80 ms)

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
            F_DEP_i = 2 * pi * epsilon_m * cells(i).radius^3 * Re_CM(i) * dEdx(idx);
            F_total = [F_DEP_i; 0] + F_drag;
            
            % Update position and velocity
            cells(i).velocity = cells(i).velocity + (F_total / ((4/3)*pi*cells(i).radius^3*rho_cell)) * dt;
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
        
        % Reposition trapped cells to 2nd or 4th pairs at end of trapping phase
        if current_time >= fusion_phase_start - dt
            for i = 1:num_cells
                if strcmp(cells(i).status, 'trapped')
                    % Find closest of 2nd (20e-6 m) or 4th (60e-6 m) pair
                    target_positions = [20e-6, 60e-6];
                    [~, closest_idx] = min(abs(target_positions - cells(i).x_pos));
                    cells(i).x_pos = target_positions(closest_idx);
                end
            end
        end
        
        % Pair trapped cells (prioritize PLT-RBC pairs)
        if current_time > 0.02
            trapped_indices = find(strcmp({cells.status}, 'trapped'));
            plt_indices = trapped_indices(strcmp({cells(trapped_indices).type}, 'PLT'));
            rbc_indices = trapped_indices(strcmp({cells(trapped_indices).type}, 'RBC'));
            
            % Pair PLT-RBC first
            num_pairs = min(length(plt_indices), length(rbc_indices));
            for i = 1:num_pairs
                idx1 = plt_indices(i);
                idx2 = rbc_indices(i);
                cells(idx1).status = 'paired';
                cells(idx2).status = 'paired';
                cells(idx1).partner_id = cells(idx2).id;
                cells(idx2).partner_id = cells(idx1).id;
                % Move paired cells together, prefer 2nd or 4th electrode
                target_pos = electrode_positions(mod(i-1, 2) + 1); % Alternate between 20e-6, 60e-6
                cells(idx1).x_pos = target_pos - spacing/4;
                cells(idx2).x_pos = target_pos + spacing/4;
            end
            
            % Pair remaining cells (same-type)
            remaining_indices = trapped_indices(~strcmp({cells(trapped_indices).status}, 'paired'));
            for i = 1:2:length(remaining_indices)-1
                idx1 = remaining_indices(i);
                idx2 = remaining_indices(i+1);
                cells(idx1).status = 'paired';
                cells(idx2).status = 'paired';
                cells(idx1).partner_id = cells(idx2).id;
                cells(idx2).partner_id = cells(idx1).id;
                % Move paired cells together
                mean_x = mean([cells(idx1).x_pos, cells(idx2).x_pos]);
                cells(idx1).x_pos = mean_x - spacing/4;
                cells(idx2).x_pos = mean_x + spacing/4;
            end
        end
        
    %% Electrofusion Phase
    else
        % Apply AC field for alignment
        if current_time < fusion_phase_start + 0.01
            for i = 1:num_cells
                if strcmp(cells(i).status, 'paired')
                    % Additional alignment force
                    partner_idx = find([cells.id] == cells(i).partner_id);
                    if ~isempty(partner_idx)
                        dx = cells(partner_idx).x_pos - cells(i).x_pos;
                        alignment_force = 1e-12 * sign(dx - spacing/10);
                        cells(i).x_pos = cells(i).x_pos + alignment_force * dt;
                    end
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
                        % Fused cells become one (visualization only)
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
    
    % Visualize every 1000 steps to improve performance
    if mod(t,1000) == 0
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

% Plot trajectories (subset for clarity)
figure;
hold on;
max_trajectories = min(50, num_cells); % Plot up to 50 trajectories
for i = 1:max_trajectories
    plot(squeeze(cell_trajectories(i,:,1))*1e6, squeeze(cell_trajectories(i,:,2))*1e6, ...
        'Color', colors.(cells(i).status), 'LineWidth', 1);
end
xlabel('X Position (μm)'); ylabel('Y Position (μm)');
title('Cell Trajectories (First 50 Cells)'); grid on;

% Efficiency metrics
num_trapped = sum(strcmp({cells.status}, 'trapped') | strcmp({cells.status}, 'paired') | ...
                  strcmp({cells.status}, 'fused') | strcmp({cells.status}, 'failed'));
num_paired = sum(strcmp({cells.status}, 'paired') | strcmp({cells.status}, 'fused') | ...
                 strcmp({cells.status}, 'failed'));
fused_pairs = sum(strcmp({cells.status}, 'fused'))/2;
fprintf('\n=== System Performance ===\n');
fprintf('Trapping Efficiency: %d/%d cells trapped (%.1f%%)\n', num_trapped, num_cells, ...
        100*num_trapped/num_cells);
fprintf('Pairing Efficiency: %d/%d trapped cells paired (%.1f%%)\n', num_paired, num_trapped, ...
        100*num_paired/num_trapped);
fprintf('Fusion Efficiency: %d/%d possible pairs fused (%.1f%%)\n', fused_pairs, ...
        floor(num_cells/2), 100*fused_pairs/floor(num_cells/2));
fprintf('Average trapping time: %.2f ms\n', fusion_phase_start*1e3);

% Additional statistics: Count PLT-RBC vs. same-type fusions
plt_rbc_fused = 0;
same_type_fused = 0;
for i = 1:num_cells
    if strcmp(cells(i).status, 'fused')
        partner_idx = find([cells.id] == cells(i).partner_id);
        if ~isempty(partner_idx)
            if strcmp(cells(i).type, 'PLT') && strcmp(cells(partner_idx).type, 'RBC') || ...
               strcmp(cells(i).type, 'RBC') && strcmp(cells(partner_idx).type, 'PLT')
                plt_rbc_fused = plt_rbc_fused + 0.5; % Count pair once
            elseif strcmp(cells(i).type, cells(partner_idx).type)
                same_type_fused = same_type_fused + 0.5; % Count pair once
            end
        end
    end
end
fprintf('PLT-RBC fused pairs: %d (%.1f%% of fused pairs)\n', plt_rbc_fused, ...
        100*plt_rbc_fused/fused_pairs);
fprintf('Same-type fused pairs: %d (%.1f%% of fused pairs)\n', same_type_fused, ...
        100*same_type_fused/fused_pairs);

%% Helper Function
function visualizeSystem(cells, channel_length, colors, current_time)
    clf;
    hold on;
    
    % Draw parallel electrodes with semicircular tips
    electrode_pos = linspace(0, channel_length, 10);
    electrode_width = 15; % Width of electrode in µm, increased for visibility
    electrode_height = 15; % Height of electrode in µm, increased for visibility
    for pos = electrode_pos
        % Top electrode (at y = 1.3e6 µm)
        rectangle('Position', [pos*1e6 - electrode_width/2, 1.3*1e6, electrode_width, electrode_height], ...
                 'Curvature', [1 0], ... % Semicircular ends along x-axis
                 'FaceColor', [0.3 0.3 0.3], ...
                 'EdgeColor', 'none');
        % Bottom electrode (at y = 0.7e6 µm)
        rectangle('Position', [pos*1e6 - electrode_width/2, 0.7*1e6, electrode_width, electrode_height], ...
                 'Curvature', [1 0], ... % Semicircular ends along x-axis
                 'FaceColor', [0.3 0.3 0.3], ...
                 'EdgeColor', 'none');
    end
    
    % Plot cells (subset for clarity)
    max_cells_plot = min(50, length(cells)); % Plot up to 50 cells
    for i = 1:max_cells_plot
        scatter(cells(i).x_pos*1e6, cells(i).y_pos*1e6, 50, ...
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
    title(sprintf('DEP-EF System with Parallel Electrodes | Time: %.1f ms (Showing %d/%d Cells)', ...
        current_time*1e3, max_cells_plot, length(cells)));
    grid on;
    drawnow;
end