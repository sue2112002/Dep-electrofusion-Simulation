clc; clear; close all;

%% ================== SYSTEM PARAMETERS ==================
% Physical constants
epsilon_0 = 8.854e-12;        % Vacuum permittivity (F/m)
epsilon_r = 80;               % Relative permittivity (water)
epsilon_m = epsilon_0 * epsilon_r; % Base medium permittivity
sigma_m = 1e-3;               % Base medium conductivity (S/m)

% Simulation domain
domain_size = 100e-6;         % 100 μm (adjusted to match code domain)
nx = 100;                     % Grid points in x
ny = 100;                     % Grid points in y
dx = domain_size / (nx-1);    % Grid spacing
dy = domain_size / (ny-1);
x = linspace(0, domain_size, nx);
y = linspace(0, domain_size, ny);
[X, Y] = meshgrid(x, y);

% Electrode parameters
electrode_potential = 10;     % V
electrode_size = 20e-6;       % 20 μm (base width or diameter)

% Configurations
configs = {'triangular', 'rectangular', 'semicircular'};

%% ================== LOAD CSV DATASET ==================
data = readtable('input_cells_dataset.csv');
num_cells = height(data);

% Extract cell properties and scale y to fit 100 μm domain
cell_x = data.x * 1e-6;       % Convert x to meters (all 0)
cell_y = data.y * 1e-6 * (domain_size / max(data.y * 1e-6)); % Scale y to 0-100 μm
cell_radius = data.radius * 1e-6; % Convert radius to meters
cell_epsilon = data.epsilon_cyt * epsilon_0; % Cytoplasmic permittivity
cell_sigma = data.sigma_cyt;   % Cytoplasmic conductivity

for cfg = 1:length(configs)
    config = configs{cfg};
    
    %% ================== ELECTRODE MASK ==================
    V = zeros(nx, ny); % Potential grid
    mask = false(nx, ny); % Electrode mask
    
    center_x = domain_size / 2;
    center_y = domain_size / 2;
    
    switch config
        case 'triangular'
            % Isosceles triangle, base 20 μm, height 20 μm
            for i = 1:nx
                for j = 1:ny
                    if (Y(j,i) >= center_y - electrode_size/2) && (Y(j,i) <= center_y + electrode_size/2)
                        half_base = (electrode_size/2) * (1 - (Y(j,i) - (center_y - electrode_size/2)) / electrode_size);
                        if (X(j,i) >= center_x - half_base) && (X(j,i) <= center_x + half_base)
                            mask(j,i) = true;
                            V(j,i) = electrode_potential;
                        end
                    end
                end
            end
            
        case 'rectangular'
            % Rectangle, 20 μm × 20 μm
            mask = (X >= center_x - electrode_size/2) & (X <= center_x + electrode_size/2) & ...
                   (Y >= center_y - electrode_size/2) & (Y <= center_y + electrode_size/2);
            V(mask) = electrode_potential;
            
        case 'semicircular'
            % Semi-circle, radius 10 μm, flat side at y = center_y
            mask = ((X - center_x).^2 + (Y - center_y).^2 <= (electrode_size/2)^2) & (Y >= center_y);
            V(mask) = electrode_potential;
    end
    
    %% ================== CELL INFLUENCE ==================
    % Initialize grids with medium properties
    epsilon_grid = ones(nx, ny) * epsilon_m;
    sigma_grid = ones(nx, ny) * sigma_m;
    
    % Adjust permittivity and conductivity in cell regions
    for c = 1:num_cells
        % Find grid points within cell radius, scaled to domain
        dist = sqrt((X - cell_x(c)).^2 + (Y - cell_y(c)).^2);
        cell_mask = dist <= cell_radius(c);
        epsilon_grid(cell_mask) = cell_epsilon(c);
        sigma_grid(cell_mask) = cell_sigma(c);
    end
    
    %% ================== SOLVE LAPLACE'S EQUATION ==================
    % Finite difference method with variable properties (simplified)
    max_iter = 1000;
    tolerance = 1e-6;
    for iter = 1:max_iter
        V_old = V;
        for i = 2:nx-1
            for j = 2:ny-1
                if ~mask(j,i)
                    % Simplified averaging (uniform medium assumed for now)
                    V(j,i) = 0.25 * (V(j+1,i) + V(j-1,i) + V(j,i+1) + V(j,i-1));
                end
            end
        end
        % Boundary conditions (grounded)
        V(:,1) = 0; V(:,end) = 0; V(1,:) = 0; V(end,:) = 0;
        % Check convergence
        if max(abs(V(:) - V_old(:))) < tolerance
            break;
        end
    end
    
    %% ================== COMPUTE ELECTRIC FIELD ==================
    [Ex, Ey] = gradient(-V, dx, dy);
    E_magnitude = sqrt(Ex.^2 + Ey.^2);
    
    %% ================== EXPORT TO CSV ==================
    % Flatten grids for CSV
    data = [X(:), Y(:), E_magnitude(:)];
    csvwrite(sprintf('electric_field_%s.csv', config), data);
    
    %% ================== VISUALIZATION ==================
    figure;
    surf(X*1e6, Y*1e6, E_magnitude/1e5, 'EdgeColor', 'none');
    xlabel('X (μm)'); ylabel('Y (μm)'); zlabel('Electric Field (×10^5 V/m)');
    title(sprintf('Electric Field Distribution (%s)', config));
    colorbar;
    
    % Overlay cell positions
    hold on;
    scatter3(cell_x*1e6, cell_y*1e6, zeros(num_cells,1), 20, 'r', 'filled');
end