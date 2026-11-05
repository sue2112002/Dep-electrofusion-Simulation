function MicrofluidicElectroporationSystem
    % Create main figure window
    f = figure('Name', 'Microfluidic Electroporation System', 'NumberTitle', 'off', ...
               'Position', [100, 100, 1200, 800], 'Color', [0.95 0.95 0.95]);
    
    % Create tab group
    tabgroup = uitabgroup(f, 'Position', [0.3 0.05 0.68 0.9]);
    
    % Tab 1: Flow and Electroporation (Original)
    tab1 = uitab(tabgroup, 'Title', 'Flow & EP Monitoring');
    hFlowAxes = axes('Parent', tab1, 'Position', [0.1, 0.55, 0.35, 0.35]);
    hElectroAxes = axes('Parent', tab1, 'Position', [0.55, 0.55, 0.35, 0.35]);
    hPulseAxes = axes('Parent', tab1, 'Position', [0.1, 0.1, 0.35, 0.35]);
    hTransitAxes = axes('Parent', tab1, 'Position', [0.55, 0.1, 0.35, 0.35]);
    
    % Tab 2: Electric Field Analysis
    tab2 = uitab(tabgroup, 'Title', 'Field Analysis');
    hEFieldAxes = axes('Parent', tab2, 'Position', [0.1, 0.1, 0.35, 0.8]);
    hFFTAxes = axes('Parent', tab2, 'Position', [0.55, 0.1, 0.35, 0.8]);
    
    % Tab 3: Cell Viability
    tab3 = uitab(tabgroup, 'Title', 'Viability');
    hViabilityAxes = axes('Parent', tab3, 'Position', [0.2, 0.2, 0.6, 0.6]);
    
    % Parameters Panel (Left Side)
    uipanel(f, 'Title', 'Parameters', 'Position', [0.02, 0.05, 0.27, 0.9]);
    
    % Microfluidic Parameters
    uicontrol(f, 'Style', 'text', 'Position', [30, 750, 150, 20], 'String', 'Channel Length (µm)');
    uicontrol(f, 'Style', 'text', 'Position', [30, 720, 150, 20], 'String', 'Channel Width (µm)');
    uicontrol(f, 'Style', 'text', 'Position', [30, 690, 150, 20], 'String', 'Flow Rate (µL/min)');
    uicontrol(f, 'Style', 'text', 'Position', [30, 660, 150, 20], 'String', 'Viscosity (cP)');
    
    hLength = uicontrol(f, 'Style', 'edit', 'Position', [180, 750, 100, 20], 'String', '1000');
    hWidth = uicontrol(f, 'Style', 'edit', 'Position', [180, 720, 100, 20], 'String', '50');
    hFlow = uicontrol(f, 'Style', 'edit', 'Position', [180, 690, 100, 20], 'String', '10');
    hVisc = uicontrol(f, 'Style', 'edit', 'Position', [180, 660, 100, 20], 'String', '1');
    
    % Electroporation Parameters
    uicontrol(f, 'Style', 'text', 'Position', [30, 600, 150, 20], 'String', 'Applied Voltage (V)');
    uicontrol(f, 'Style', 'text', 'Position', [30, 570, 150, 20], 'String', 'Pulse Frequency (Hz)');
    uicontrol(f, 'Style', 'text', 'Position', [30, 540, 150, 20], 'String', 'Baseline Current (µA)');
    
    hVoltage = uicontrol(f, 'Style', 'edit', 'Position', [180, 600, 100, 20], 'String', '3');
    hFreq = uicontrol(f, 'Style', 'edit', 'Position', [180, 570, 100, 20], 'String', '1000');
    hCurrent = uicontrol(f, 'Style', 'edit', 'Position', [180, 540, 100, 20], 'String', '30');
    
    % Cell Parameters
    uicontrol(f, 'Style', 'text', 'Position', [30, 480, 150, 20], 'String', 'Cell Diameter (µm)');
    uicontrol(f, 'Style', 'text', 'Position', [30, 450, 150, 20], 'String', 'Cell Count');
    
    hCellDiam = uicontrol(f, 'Style', 'edit', 'Position', [180, 480, 100, 20], 'String', '15');
    hCellCount = uicontrol(f, 'Style', 'edit', 'Position', [180, 450, 100, 20], 'String', '5');
    
    % Control buttons
    uicontrol(f, 'Style', 'pushbutton', 'String', 'Run Simulation', ...
              'Position', [80, 350, 150, 30], 'Callback', @runSimulation);
    uicontrol(f, 'Style', 'pushbutton', 'String', 'Reset', ...
              'Position', [80, 310, 150, 30], 'Callback', @resetSimulation);
    
    % Simulation data storage
    simData = struct();
    
    function runSimulation(~,~)
        % Get all parameters
        channel_length = str2double(get(hLength, 'String')) * 1e-6;
        channel_width = str2double(get(hWidth, 'String')) * 1e-6;
        flow_rate = str2double(get(hFlow, 'String')) * 1e-9/60;
        viscosity = str2double(get(hVisc, 'String')) * 1e-3;
        
        applied_voltage = str2double(get(hVoltage, 'String'));
        pulse_freq = str2double(get(hFreq, 'String'));
        baseline_current = str2double(get(hCurrent, 'String')) * 1e-6;
        
        cell_diameter = str2double(get(hCellDiam, 'String')) * 1e-6;
        cell_count = str2double(get(hCellCount, 'String'));
        
        % Calculate flow parameters
        [y, velocity_profile] = calculate_velocity_rectangular(channel_length, channel_width, flow_rate, viscosity);
        avg_velocity = mean(velocity_profile);
        
        if avg_velocity <= 0 || isnan(avg_velocity) || isinf(avg_velocity)
            errordlg('Average velocity is invalid. Check your input parameters.','Simulation Error');
            return;
        end
        
        transit_time = channel_length / avg_velocity;
        transit_time = min(transit_time, 2); % cap to 2 seconds to avoid large arrays
        
        % Simulate electroporation
        Fs = 10000;
        t = 0:1/Fs:transit_time*1.5;
        V_signal = applied_voltage * square(2 * pi * pulse_freq * t);
        
        try
            [current, impedance, event_times] = simulate_cell_events(t, cell_count, transit_time, baseline_current, applied_voltage);
            simData.current = current;
            simData.impedance = impedance;
            simData.event_times = event_times;
        catch ME
            disp('Error in simulate_cell_events:');
            disp(ME.message);
            
            % Fallback dummy data
            simData.current = baseline_current * ones(size(t));
            simData.impedance = ones(size(t));
            simData.event_times = [];
        end
        
        % Store simulation data
        simData.t = t;
        simData.velocity_profile = velocity_profile;
        simData.y = y;
        simData.V_signal = V_signal;
        
        % Update all plots
        updatePlots();
    end

    function updatePlots()
        if ~isfield(simData, 't'), return; end
        
        % Tab 1: Original Plots + Extensions
        % Flow profile
        plot(hFlowAxes, simData.velocity_profile*1000, simData.y*1e6, 'b', 'LineWidth', 2);
        title(hFlowAxes, 'Microchannel Velocity Profile');
        xlabel(hFlowAxes, 'Velocity (mm/s)');
        ylabel(hFlowAxes, 'Channel Height (µm)');
        grid(hFlowAxes, 'on');
        
        % Electroporation monitoring
        yyaxis(hElectroAxes, 'left');
        plot(hElectroAxes, simData.t, simData.current*1e6, 'b', 'LineWidth', 1.5);
        ylabel(hElectroAxes, 'Current (µA)');
        yyaxis(hElectroAxes, 'right');
        plot(hElectroAxes, simData.t, simData.impedance, 'r', 'LineWidth', 1.5);
        ylabel(hElectroAxes, 'Impedance (Ω)');
        title(hElectroAxes, 'Electroporation Monitoring');
        grid(hElectroAxes, 'on');
        
        % Mark electroporation events
        if ~isempty(simData.event_times)
            for i = 1:length(simData.event_times)
                xline(hElectroAxes, simData.event_times(i), '--r', sprintf('EP %d', i));
            end
        end
        
        % Pulse Response (zoomed)
        plot(hPulseAxes, simData.t(1:1000), simData.V_signal(1:1000), 'k');
        hold(hPulseAxes, 'on');
        plot(hPulseAxes, simData.t(1:1000), simData.current(1:1000)*1e6, 'r');
        title(hPulseAxes, 'Pulse Response (First 100ms)');
        legend(hPulseAxes, {'Voltage', 'Current (µA)'}, 'Location', 'northeast');
        hold(hPulseAxes, 'off');
        
        % Transit Time Histogram
        if length(simData.event_times) > 1
            transit_times = diff(simData.event_times);
            histogram(hTransitAxes, transit_times, 'BinWidth', 0.05);
            title(hTransitAxes, 'Cell Transit Times');
            xlabel(hTransitAxes, 'Time (s)');
            ylabel(hTransitAxes, 'Count');
        end
        
        % Tab 2: Field Analysis
        % Electric Field (simplified approximation)
        [Ex, Ey] = gradient(-simData.V_signal(1:100:end)); % Downsample for performance
        E_mag = sqrt(Ex.^2 + Ey.^2);
        imagesc(hEFieldAxes, E_mag);
        title(hEFieldAxes, 'Electric Field Magnitude (V/m)');
        colorbar(hEFieldAxes);
        
        % FFT Analysis
        L = length(simData.impedance);
        Y = fft(simData.impedance);
        P2 = abs(Y/L);
        P1 = P2(1:L/2+1);
        f = (0:(L/2)) * (1/(simData.t(2)-simData.t(1))) / L;
        plot(hFFTAxes, f, P1);
        title(hFFTAxes, 'Impedance Spectrum');
        xlabel(hFFTAxes, 'Frequency (Hz)');
        ylabel(hFFTAxes, '|Z(f)|');
        
        % Tab 3: Viability (mock data)
        voltages = [1, 3, 5, 10];
        viability = [100, 80, 50, 20]; % Example data
        plot(hViabilityAxes, voltages, viability, '-o');
        title(hViabilityAxes, 'Cell Viability vs Voltage');
        xlabel(hViabilityAxes, 'Voltage (V)');
        ylabel(hViabilityAxes, 'Viability (%)');
        ylim(hViabilityAxes, [0 100]);
    end

    function resetSimulation(~,~)
        % Reset all parameter fields to defaults
        set(hLength, 'String', '1000');
        set(hWidth, 'String', '50');
        set(hFlow, 'String', '10');
        set(hVisc, 'String', '1');
        set(hVoltage, 'String', '3');
        set(hFreq, 'String', '1000');
        set(hCurrent, 'String', '30');
        set(hCellDiam, 'String', '15');
        set(hCellCount, 'String', '5');
        
        % Clear simulation data
        simData = struct();
        
        % Clear all axes
        cla(hFlowAxes);
        cla(hElectroAxes);
        cla(hPulseAxes);
        cla(hTransitAxes);
        cla(hEFieldAxes);
        cla(hFFTAxes);
        cla(hViabilityAxes);
    end

    function [current, impedance, event_times] = simulate_cell_events(t, cell_count, transit_time, baseline_current, applied_voltage)
        % Initialize outputs
        current = baseline_current * ones(size(t));
        impedance = ones(size(t));
        event_times = zeros(1, cell_count);
        
        % Randomize cell arrival times within the transit period
        start_times = sort(rand(1, cell_count) * transit_time);
        event_times = start_times;
        
        % Add current spikes for each cell
        for i = 1:cell_count
            % Find the time index closest to this event
            [~, idx] = min(abs(t - start_times(i)));
            
            % Create a current spike with some random variation
            spike_magnitude = baseline_current * (1 + 0.5*rand()); % 50% variation
            spike_duration = 0.005 + 0.002*randn(); % 5ms ± 2ms
            
            % Apply the spike
            spike_samples = round(spike_duration * (t(2)-t(1))^-1);
            idx_range = idx:min(idx+spike_samples-1, length(t));
            current(idx_range) = spike_magnitude;
            
            % Corresponding impedance drop
            impedance(idx_range) = 0.7 + 0.3*rand(size(idx_range));
        end
        
        % Add some noise to the signals
        current = current + baseline_current * 0.02 * randn(size(t));
        impedance = impedance + 0.02 * randn(size(t));
    end

    function [y, velocity_profile] = calculate_velocity_rectangular(channel_length, channel_width, flow_rate, viscosity)
        % Constants
        h = channel_width / 2;  % Half-height of channel
        y = linspace(-h, h, 100);  % vertical positions from bottom to top
        deltaP = (12 * viscosity * channel_length * flow_rate) / (channel_width^3); % Pressure drop (simplified)
        
        % Velocity profile for pressure-driven flow between parallel plates
        velocity_profile = (1 / (2 * viscosity)) * deltaP * (h^2 - y.^2);
    end
end