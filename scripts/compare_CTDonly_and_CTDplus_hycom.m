%% deep casts only 
addpath('H:\PhD work\Miguel\seawater_ver3_3.1')
rho = 1040; %(kg/m3) 
g = 9.8; %(m/s2) 
% Load CTD data
ctd_file = 'deep_casts_2021.csv';
ctd_data = readtable(ctd_file);

% Convert CTD DateTime to MATLAB datetime format
ctd_data.DateTime = datetime(ctd_data.DateTime, 'InputFormat', 'M/d/yyyy');

% Get unique CTD dates
unique_ctd_dates = unique(ctd_data.DateTime);
%% Initialize results table with an additional 'Source' column
results_table = table('Size', [0 3], 'VariableTypes', {'datetime', 'double', 'string'}, ...
                      'VariableNames', {'DateTime', 'bpga_combine', 'CastType'});

% Loop through each unique date
for i = 1:length(unique_ctd_dates)
    % Filter the data for the current date
    current_date = unique_ctd_dates(i);
    daily_data = ctd_data(ctd_data.DateTime == current_date, :);

    % Extract pressure, salinity, temperature, and source for the current day
    Pressure = daily_data.Pressure;
    Salinity = daily_data.Salinity;
    Temperature = daily_data.Temperature;
 
    % Remove rows with NaN values
    valid_idx = ~isnan(Pressure) & ~isnan(Salinity) & ~isnan(Temperature);
    Pressure = Pressure(valid_idx);
    Salinity = Salinity(valid_idx);
    Temperature = Temperature(valid_idx);

    % Skip if no valid data is available for the current date
    if isempty(Pressure)
        disp(['No valid data for ', datestr(current_date), ', skipping.']);
        continue;
    end

    % Determine cast type based on pressure
    if max(Pressure) > 2000
        cast_type = "Deep Cast";
    else
        cast_type = "Shallow Cast";
    end

    % Calculate geopotential anomaly (ga) and bottom pressure anomaly
    try
        ga = sw_gpan(Salinity, Temperature, Pressure); % Geopotential anomaly
        bpga_combine = rho * g * ga(end) / g; % Bottom pressure anomaly at deepest level
    catch ME
        warning(['Failed to calculate geopotential anomaly for ', datestr(current_date), ': ', ME.message]);
        continue;
    end

    % Add the result to the results table
    new_row = table(current_date, bpga_combine, cast_type, 'VariableNames', {'DateTime', 'bpga_combine', 'CastType'});
    results_table = [results_table; new_row];
end

% Save results to a CSV file
writetable(results_table, 'bpga_CTD_deep_casts_2021.csv');
disp('Results saved to bpga_CTDonly_results_2020_withtag.csv.');

%% 
% Calculate demeaned data
demeaned_data = results_table.bpga_combine - nanmean(results_table.bpga_combine);

% Assign colors based on CastType
deep_cast_idx = results_table.CastType == "Deep Cast";
shallow_cast_idx = results_table.CastType == "Shallow Cast";

% Create scatter plot with different colors for CastType
figure;
hold on;
scatter(results_table.DateTime(deep_cast_idx), demeaned_data(deep_cast_idx) / 10000, ...
        'filled', 'MarkerFaceColor', 'b', 'DisplayName', 'Deep Cast');
scatter(results_table.DateTime(shallow_cast_idx), demeaned_data(shallow_cast_idx) / 10000, ...
        'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'Shallow Cast');
hold off;

% Add title, labels, and legend
title('BPGA Combine Results (CTD Data) Over Time');
xlabel('Date');
ylabel('Pressure (dbar)');
grid on;
legend('show', 'Location', 'best');

% Format the x-axis for better datetime display
datetick('x', 'yyyy-mm-dd', 'keepticks');
xtickangle(45); % Rotate x-axis labels for better readability
%% shallow casts only 
addpath('H:\PhD work\Miguel\seawater_ver3_3.1')
rho = 1040; %(kg/m3) 
g = 9.8; %(m/s2) 
% Load CTD data
ctd_file = 'shallow_casts_2021.csv';
ctd_data = readtable(ctd_file);

% Convert CTD DateTime to MATLAB datetime format
ctd_data.DateTime = datetime(ctd_data.DateTime, 'InputFormat', 'M/d/yyyy');

% Get unique CTD dates
unique_ctd_dates = unique(ctd_data.DateTime);
%% Initialize results table with an additional 'Source' column
results_table_shallow = table('Size', [0 3], 'VariableTypes', {'datetime', 'double', 'string'}, ...
                      'VariableNames', {'DateTime', 'bpga_combine', 'CastType'});

% Loop through each unique date
for i = 1:length(unique_ctd_dates)
    % Filter the data for the current date
    current_date = unique_ctd_dates(i);
    daily_data = ctd_data(ctd_data.DateTime == current_date, :);

    % Extract pressure, salinity, temperature, and source for the current day
    Pressure = daily_data.Pressure;
    Salinity = daily_data.Salinity;
    Temperature = daily_data.Temperature;
 
    % Remove rows with NaN values
    valid_idx = ~isnan(Pressure) & ~isnan(Salinity) & ~isnan(Temperature);
    Pressure = Pressure(valid_idx);
    Salinity = Salinity(valid_idx);
    Temperature = Temperature(valid_idx);

    % Skip if no valid data is available for the current date
    if isempty(Pressure)
        disp(['No valid data for ', datestr(current_date), ', skipping.']);
        continue;
    end

    % Determine cast type based on pressure
    if max(Pressure) > 2000
        cast_type = "Deep Cast";
    else
        cast_type = "Shallow Cast";
    end

    % Calculate geopotential anomaly (ga) and bottom pressure anomaly
    try
        ga = sw_gpan(Salinity, Temperature, Pressure); % Geopotential anomaly
        bpga_combine = rho * g * ga(end) / g; % Bottom pressure anomaly at deepest level
    catch ME
        warning(['Failed to calculate geopotential anomaly for ', datestr(current_date), ': ', ME.message]);
        continue;
    end

    % Add the result to the results table
    new_row = table(current_date, bpga_combine, cast_type, 'VariableNames', {'DateTime', 'bpga_combine', 'CastType'});
    results_table_shallow = [results_table_shallow; new_row];
end

% Save results to a CSV file
writetable(results_table_shallow, 'bpga_CTD_shallow_casts_2021.csv');
disp('Results saved to bpga_CTDonly_results_2020_withtag.csv.');
%% shallow + deep ctd data 
addpath('H:\PhD work\Miguel\seawater_ver3_3.1')
rho = 1040; %(kg/m3) 
g = 9.8; %(m/s2) 
% Load CTD data
ctd_file = 'ctd_salinity_and_temperature_2021_combine_deep_data_filtered_time.csv';
ctd_data = readtable(ctd_file);

% Convert CTD DateTime to MATLAB datetime format
ctd_data.DateTime = datetime(ctd_data.DateTime, 'InputFormat', 'M/d/yyyy');

% Get unique CTD dates
unique_ctd_dates = unique(ctd_data.DateTime);
%% Initialize results table with an additional 'Source' column
results_table_shallow_deep = table('Size', [0 3], 'VariableTypes', {'datetime', 'double', 'string'}, ...
                      'VariableNames', {'DateTime', 'bpga_combine', 'CastType'});

% Loop through each unique date
for i = 1:length(unique_ctd_dates)
    % Filter the data for the current date
    current_date = unique_ctd_dates(i);
    daily_data = ctd_data(ctd_data.DateTime == current_date, :);

    % Extract pressure, salinity, temperature, and source for the current day
    Pressure = daily_data.Pressure;
    Salinity = daily_data.Salinity;
    Temperature = daily_data.Temperature;
 
    % Remove rows with NaN values
    valid_idx = ~isnan(Pressure) & ~isnan(Salinity) & ~isnan(Temperature);
    Pressure = Pressure(valid_idx);
    Salinity = Salinity(valid_idx);
    Temperature = Temperature(valid_idx);

    % Skip if no valid data is available for the current date
    if isempty(Pressure)
        disp(['No valid data for ', datestr(current_date), ', skipping.']);
        continue;
    end

    % Determine cast type based on pressure
    if max(Pressure) > 2000
        cast_type = "Deep Cast";
    else
        cast_type = "Shallow Cast";
    end

    % Calculate geopotential anomaly (ga) and bottom pressure anomaly
    try
        ga = sw_gpan(Salinity, Temperature, Pressure); % Geopotential anomaly
        bpga_combine = rho * g * ga(end) / g; % Bottom pressure anomaly at deepest level
    catch ME
        warning(['Failed to calculate geopotential anomaly for ', datestr(current_date), ': ', ME.message]);
        continue;
    end

    % Add the result to the results table
    new_row = table(current_date, bpga_combine, cast_type, 'VariableNames', {'DateTime', 'bpga_combine', 'CastType'});
    results_table_shallow_deep = [results_table_shallow_deep; new_row];
end

% Save results to a CSV file
writetable(results_table_shallow_deep, 'bpga_CTD_shallow_and_deep_combine_casts_2021.csv');
disp('Results saved to bpga_CTDonly_results_2020_withtag.csv.');

%%
% Calculate demeaned data
demeaned_data_shallow = results_table_shallow.bpga_combine - nanmean(results_table_shallow.bpga_combine);
demeaned_data_deep = results_table.bpga_combine - nanmean(results_table.bpga_combine);
demeaned_data_shallow_deep = results_table_shallow_deep.bpga_combine - nanmean(results_table_shallow_deep.bpga_combine);

% Find the common DateTime across all three datasets
common_dates_12 = intersect(results_table_shallow.DateTime, results_table.DateTime);  % Shallow & Deep
common_dates_all = intersect(common_dates_12, results_table_shallow_deep.DateTime);   % Shallow, Deep & Shallow-Deep

% Find the indices for common dates in each dataset
[~, idx_shallow] = ismember(common_dates_all, results_table_shallow.DateTime);
[~, idx_deep] = ismember(common_dates_all, results_table.DateTime);
[~, idx_shallow_deep] = ismember(common_dates_all, results_table_shallow_deep.DateTime);

% Filter the data based on common dates
common_shallow_data = demeaned_data_shallow(idx_shallow);
common_deep_data = demeaned_data_deep(idx_deep);
common_shallow_deep_data = demeaned_data_shallow_deep(idx_shallow_deep);

% Define a colormap for unique colors for each shared date
num_points = length(common_dates_all);
colors = lines(num_points);  % 'lines' colormap provides distinct colors

% Create scatter plot for common shallow, deep, and shallow-deep data with the same color
figure;
hold on;

for i = 1:num_points
    % Plot Shallow Cast data
    scatter(common_dates_all(i), common_shallow_data(i) / 10000, ...
            60, 'o', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', 'k', ...
            'DisplayName', 'Shallow Cast');
    
    % Plot Deep Cast data
    scatter(common_dates_all(i), common_deep_data(i) / 10000, ...
            60, '^', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', 'k', ...
            'DisplayName', 'Deep Cast');
    
    % Plot Shallow + Deep Cast data
    scatter(common_dates_all(i), common_shallow_deep_data(i) / 10000, ...
            60, 's', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', 'k', ...
            'DisplayName', 'Shallow + Deep Cast');
end

hold off;

% Add title, labels, and legend
title('Common BPGA Combine Results (CTD Data) Over Shared Time', 'FontSize', 14);
xlabel('Date', 'FontSize', 12);
ylabel('Pressure Anomaly (dbar)', 'FontSize', 12);
grid on;

% Combine legend entries for shallow, deep, and shallow+deep casts
legend({'Shallow Cast', 'Deep Cast', 'Shallow + Copy Deep Cast'}, 'Location', 'best');

% Format the x-axis for better datetime display
datetick('x', 'yyyy-mm-dd', 'keepticks');
xtickangle(45);  % Rotate x-axis labels for better readability

%%



%%
ctd_combine_hycom= 'combined_ctd_hycom_profiles_2020.csv';
ctd_combine_hycom_data = readtable(ctd_combine_hycom);
unique_ctd_hycom_dates = unique(ctd_combine_hycom_data.DateTime);
%% Initialize results table with an additional 'Source' column

results_table_hycom = table('Size', [0 3], 'VariableTypes', {'datetime', 'double', 'string'}, ...
                      'VariableNames', {'DateTime', 'bpga_combine_hycom', 'CastType'});

% Loop through each unique date
for i = 1:length(unique_ctd_hycom_dates)
    % Filter the data for the current date
    current_date = unique_ctd_hycom_dates(i);
    daily_data = ctd_combine_hycom_data(ctd_combine_hycom_data.DateTime == current_date, :);

    % Extract pressure, salinity, temperature, and source for the current day
    Pressure = daily_data.Pressure;
    Salinity = daily_data.Salinity;
    Temperature = daily_data.Temperature;
 

    % Remove rows with NaN values
    valid_idx = ~isnan(Pressure) & ~isnan(Salinity) & ~isnan(Temperature);
    Pressure = Pressure(valid_idx);
    Salinity = Salinity(valid_idx);
    Temperature = Temperature(valid_idx);


    % Skip if no valid data is available for the current date
    if isempty(Pressure)
        disp(['No valid data for ', datestr(current_date), ', skipping.']);
        continue;
    end
    % Determine cast type based on pressure
    if max(Pressure) > 2000
        cast_type = "Deep Cast";
    else
        cast_type = "Shallow Cast";
    end

    % Calculate geopotential anomaly (ga) and bottom pressure anomaly
    try
        ga = sw_gpan(Salinity, Temperature, Pressure); % Geopotential anomaly
        bpga_combine_hycom = rho * g * ga(end) / g; % Bottom pressure anomaly at deepest level
    catch ME
        warning(['Failed to calculate geopotential anomaly for ', datestr(current_date), ': ', ME.message]);
        continue;
    end

    % Add the result to the results table (using the first source tag for simplicity)
    new_row = table(current_date, bpga_combine_hycom,cast_type, 'VariableNames', {'DateTime', 'bpga_combine_hycom','CastType'});
    results_table_hycom = [results_table_hycom; new_row];
end
%%
demeaned_data_hycom = results_table_hycom.bpga_combine_hycom-nanmean(results_table_hycom.bpga_combine_hycom);
% Create scatter plot with datetime labels
% Assign colors based on CastType
deep_cast_idx = results_table.CastType == "Deep Cast";
shallow_cast_idx = results_table.CastType == "Shallow Cast";

% Create scatter plot with different colors for CastType
figure;
subplot(211)
hold on;
scatter(results_table_hycom.DateTime(deep_cast_idx), demeaned_data_hycom(deep_cast_idx) / 10000, ...
        'filled', 'MarkerFaceColor', 'b', 'DisplayName', 'Deep Cast');
scatter(results_table_hycom.DateTime(shallow_cast_idx), demeaned_data_hycom(shallow_cast_idx) / 10000, ...
        'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'Shallow Cast');
hold off;

legend('show', 'Location', 'southeast');

title('Demeaned Steric Pressure Results (CTD + HYCOM Data) Over Time');
xlabel('Date');
ylabel('Pressure (dbar)');
grid on;

% Format the x-axis for better datetime display
datetick('x', 'yyyy-mm-dd', 'keepticks');
xtickangle(45); % Rotate x-axis labels for better readability
subplot(212)
hold on;
scatter(results_table.DateTime(deep_cast_idx), demeaned_data(deep_cast_idx) / 10000, ...
        'filled', 'MarkerFaceColor', 'b', 'DisplayName', 'Deep Cast');
scatter(results_table.DateTime(shallow_cast_idx), demeaned_data(shallow_cast_idx) / 10000, ...
        'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'Shallow Cast');
hold off;

% Add title, labels, and legend
title('Demeaned Steric Pressure Results (CTD Data only) Over Time');
xlabel('Date');
ylabel('Pressure (dbar)');
grid on;
legend('show', 'Location', 'southeast');

% Format the x-axis for better datetime display
datetick('x', 'yyyy-mm-dd', 'keepticks');
xtickangle(45); % Rotate x-axis labels for better readability
%% no demean steric pressure data
figure;
subplot(212)
deep_cast_idx = results_table.CastType == "Deep Cast";
shallow_cast_idx = results_table.CastType == "Shallow Cast";
hold on;
scatter(results_table.DateTime(deep_cast_idx), results_table.bpga_combine(deep_cast_idx) / 10000, ...
        'filled', 'MarkerFaceColor', 'b', 'DisplayName', 'Deep Cast');
scatter(results_table.DateTime(shallow_cast_idx), results_table.bpga_combine(shallow_cast_idx) / 10000, ...
        'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'Shallow Cast');
hold off;
title('Steric Pressure Results (CTD Data only) Over Time')
legend('show', 'Location', 'southeast');
xlabel('Date');
ylabel('Pressure (dbar)');
grid on;


deep_cast_idx = results_table.CastType == "Deep Cast";
shallow_cast_idx = results_table.CastType == "Shallow Cast";

% Create scatter plot with different colors for CastType
subplot(211)
hold on;
scatter(results_table_hycom.DateTime(deep_cast_idx), results_table_hycom.bpga_combine_hycom(deep_cast_idx) / 10000, ...
        'filled', 'MarkerFaceColor', 'b', 'DisplayName', 'Deep Cast');
scatter(results_table_hycom.DateTime(shallow_cast_idx), results_table_hycom.bpga_combine_hycom(shallow_cast_idx) / 10000, ...
        'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'Shallow Cast');
hold off;

legend('show', 'Location', 'southeast');
title('Steric Pressure Results (CTD + HYCOM Data) Over Time')
xlabel('Date');
ylabel('Pressure (dbar)');
grid on;

%% plot CTD with HYCOM profile for PTS
% Load the combined CTD and HYCOM data
ctd_combine_hycom = 'combined_ctd_hycom_profiles_2020.csv';
ctd_combine_hycom_data = readtable(ctd_combine_hycom);

% Convert DateTime to MATLAB datetime format
ctd_combine_hycom_data.DateTime = datetime(ctd_combine_hycom_data.DateTime, 'InputFormat', 'M/d/yyyy');

% Separate data based on source
ctd_data = ctd_combine_hycom_data(strcmp(ctd_combine_hycom_data.Source, 'CTD'), :);
hycom_data = ctd_combine_hycom_data(strcmp(ctd_combine_hycom_data.Source, 'HYCOM'), :);

% Plot Pressure vs. Time
figure;
scatter(ctd_data.DateTime, ctd_data.Pressure, 'filled', 'MarkerFaceColor', 'b', 'DisplayName', 'CTD');
hold on;
scatter(hycom_data.DateTime, hycom_data.Pressure, 'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'HYCOM');
title('Pressure vs. Time (CTD and HYCOM)');
xlabel('Time');
ylabel('Pressure (dbar)');
set(gca, 'YDir', 'reverse'); % Reverse Y-axis for oceanographic convention
legend('show', 'Location', 'best');
grid on;

% Plot Temperature vs. Time
figure;
scatter(ctd_data.DateTime, ctd_data.Temperature, 'filled', 'MarkerFaceColor', 'b', 'DisplayName', 'CTD');
hold on;
scatter(hycom_data.DateTime, hycom_data.Temperature, 'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'HYCOM');
title('Temperature vs. Time (CTD and HYCOM)');
xlabel('Time');
ylabel('Temperature (°C)');
legend('show', 'Location', 'best');
grid on;

% Plot Salinity vs. Time
figure;
scatter(ctd_data.DateTime, ctd_data.Salinity, 'filled', 'MarkerFaceColor', 'b', 'DisplayName', 'CTD');
hold on;
scatter(hycom_data.DateTime, hycom_data.Salinity, 'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'HYCOM');
title('Salinity vs. Time (CTD and HYCOM)');
xlabel('Time');
ylabel('Salinity (PSU)');
legend('show', 'Location', 'best');
grid on;
%%
%% Plot Layer-by-Layer Steric Height (bpga) vs Depth for Deep Casts

%% Plot steric height with depth for CTD data
figure;
subplot(121)

% Define color map for distinct profiles
colors = lines(length(unique_ctd_dates));

% Initialize storage for bpga profiles and pressure profiles
bpga_profiles = struct();
pressure_profiles = struct();

% Loop through each unique date
for i = 1:length(unique_ctd_dates)
    current_date = unique_ctd_dates(i);
    
    % Filter data for the current date
    daily_data = ctd_data(ctd_data.DateTime == current_date, :);
    
    % Extract valid data (removing NaNs)
    valid_idx = ~isnan(daily_data.Pressure) & ~isnan(daily_data.Salinity) & ~isnan(daily_data.Temperature);
    Pressure = daily_data.Pressure(valid_idx);
    Salinity = daily_data.Salinity(valid_idx);
    Temperature = daily_data.Temperature(valid_idx);
    
    % Skip if no valid data or if it's not a deep cast
    if isempty(Pressure) || max(Pressure) <= 2000
        continue;  % Skip casts that don't reach deeper than 2000 dbar
    end

    % Calculate geopotential anomaly (ga) for each layer
    try
        ga = sw_gpan(Salinity, Temperature, Pressure); % Geopotential anomaly

        % Calculate steric height difference between consecutive 2-meter layers
        layer_bpga = diff(ga);  % Difference between adjacent layers
        layer_pressure = Pressure(2:end);  % Corresponding depth for the differences

        % Convert to dbar
        layer_bpga_dbar = rho * layer_bpga / 10000;

        % Plot the layer-by-layer bpga profile against depth
        plot(layer_bpga_dbar, -layer_pressure, 'LineWidth', 1.5, 'Color', colors(i, :));  % Negative Pressure for depth
        hold on;

        % Save bpga and pressure profiles with valid field names
        date_str = ['D_', datestr(current_date, 'yyyy_mm_dd')];  % Prefix 'D_' to avoid invalid field names
        bpga_profiles.(date_str) = layer_bpga_dbar;
        pressure_profiles.(date_str) = -layer_pressure;

    catch ME
        warning(['Failed to calculate bpga for ', datestr(current_date), ': ', ME.message]);
        continue;
    end
end

% Customize the plot
xlabel('Layer-by-Layer Geopotential Anomaly (bpga) (dbar)');
ylabel('Depth (m)');
title('Deep Casts: bpga Profile vs Depth');
grid on;
set(gca, 'FontSize', 18); % Set font size for axes

% Add legend for all deep casts (matching plotted lines)
h = findobj(gca, 'Type', 'Line');  % Get all plotted lines
legend_entries = {};
for i = 1:length(unique_ctd_dates)
    current_date = unique_ctd_dates(i);
    daily_data = ctd_data(ctd_data.DateTime == current_date, :);
    
    % Check if the cast is deep (>2000 dbar)
    if ~isempty(daily_data.Pressure) && max(daily_data.Pressure) > 2000
        legend_entries{end+1} = datestr(current_date, 'yyyy-mm-dd');
    end
end
legend(h(end:-1:1), legend_entries, 'Location', 'eastoutside', 'FontSize', 10);  % Reverse 'h' to match plot order

hold off;
%%
figure
% Zoomed-In Plot for the Upper 1000m Depth
subplot(121)

% Loop through each unique date again for the zoomed-in plot
for i = 1:length(unique_ctd_dates)
    current_date = unique_ctd_dates(i);
    
    % Filter data for the current date
    daily_data = ctd_data(ctd_data.DateTime == current_date, :);
    
    % Extract valid data (removing NaNs)
    valid_idx = ~isnan(daily_data.Pressure) & ~isnan(daily_data.Salinity) & ~isnan(daily_data.Temperature);
    Pressure = daily_data.Pressure(valid_idx);
    Salinity = daily_data.Salinity(valid_idx);
    Temperature = daily_data.Temperature(valid_idx);
    
    % Skip if no valid data or if it's not a deep cast
    if isempty(Pressure) || max(Pressure) <= 2000
        continue;  % Skip casts that don't reach deeper than 2000 dbar
    end

    % Calculate geopotential anomaly (ga) for each layer
    try
        ga = sw_gpan(Salinity, Temperature, Pressure); % Geopotential anomaly

        % Calculate steric height difference between consecutive 2-meter layers
        layer_bpga = diff(ga);  % Difference between adjacent layers
        layer_pressure = Pressure(2:end);  % Corresponding depth for the differences

        % Convert to dbar
        layer_bpga_dbar = rho * layer_bpga / 10000;

        % Plot the layer-by-layer bpga profile against depth
        %plot(layer_bpga_dbar, -layer_pressure, 'LineWidth', 1.5, 'Color', colors(i, :));  % Negative Pressure for depth
        plot(layer_bpga_dbar, -layer_pressure, 'LineWidth', 1.5, 'Color', 'r');
        hold on;

    catch ME
        warning(['Failed to calculate bpga for ', datestr(current_date), ': ', ME.message]);
        continue;
    end
end

% Customize the zoomed-in plot
xlabel('Layer-by-Layer Geopotential Anomaly (bpga) (dbar)');
ylabel('Depth (m)');
title('Deep Casts: bpga Profile (0–1000m)');
grid on;
set(gca, 'FontSize', 18); % Set font size for axes
%ylim([-1000 0]);  % Zoom into the upper 1000m

% Add legend for zoomed-in plot
h = findobj(gca, 'Type', 'Line');
%legend(h(end:-1:1), legend_entries, 'Location', 'eastoutside', 'FontSize', 10);

%hold off;

%% Save Data to .mat File
save('bpga_profiles_with_pressure.mat', 'bpga_profiles', 'pressure_profiles', 'unique_ctd_dates');
disp('Layer-by-layer bpga profiles and pressures have been saved.');
%% Plot the same steric height pressure vs depth profile with hycom data
A = 'H:\aloha\ALOHA Cabled Observatory Database\HYCOM\sal\2020_all';
B = 'H:\aloha\ALOHA Cabled Observatory Database\HYCOM\tem\2020_all';
C = 'H:\aloha\ALOHA Cabled Observatory Database\HYCOM\ssh\2020_all';

% Add the required library for seawater calculations
addpath('H:\PhD work\Miguel\seawater_ver3_3.1')

% Initialize parameters
rho = 1040; %(kg/m3) 
g = 9.8; %(m/s2) 
referenceDate = datetime('2000-01-01 00:00:00'); % Reference date for time conversion

% Initialize arrays for concatenating data
S_all = [];
T_all = [];
ssh_all = [];
time_all = [];

% Loop over each day of the year to read and concatenate data
for month = 1:12
    for day = 1:31
        % Ensure that day exists within the month
        if (month == 2 && day > 28) || (ismember(month, [4, 6, 9, 11]) && day > 30)
            continue;
        end
        
        % Format the day into month_day (e.g., 01_01, 02_15, etc.)
        day_str = sprintf('2020_%02d_%02d', month, day);

               % Define file paths for salinity, temperature, and SSH
        S_file = fullfile(A, [day_str, '.nc']);
        T_file = fullfile(B, [day_str, '.nc']);
        ssh_file = fullfile(C, [day_str, '.nc']);
        
        % Check if all files exist for the current day
        if ~exist(S_file, 'file') || ~exist(T_file, 'file') || ~exist(ssh_file, 'file')
            % If any file is missing, skip to the next day
            fprintf('File missing for %s, skipping.\n', day_str);
            continue;
        end
        
        % Read salinity, temperature, SSH, and time data for the specific day
        S = ncread(fullfile(A, [day_str, '.nc']), 'salinity');
        T = ncread(fullfile(B, [day_str, '.nc']), 'water_temp');
        ssh = ncread(fullfile(C, [day_str, '.nc']), 'surf_el');
        timeData = ncread(fullfile(C, [day_str, '.nc']), 'time');
        
        % Convert timeData to datetime format
        Timedate = referenceDate + hours(timeData);
        
        % Concatenate the data along the time dimension (4th dimension for S and T, 3rd for SSH)
        S_all = cat(4, S_all, S);
        T_all = cat(4, T_all, T);
        ssh_all = cat(3, ssh_all, ssh);
        
        % Concatenate time data
        time_all = [time_all; Timedate];
    end
end
%% Constants
time_count = length(time_all);  % Set the time_count as the length of the time_all array
% Initialize storage for bpga profiles and pressure profiles
hycom_bpga_profiles = struct();
pressure_profiles = struct();

% Physical constants
rho = 1040; % Seawater density (kg/m^3)
g = 9.8;    % Gravitational acceleration (m/s^2)

% Read latitude, longitude, and depth from a sample file
Lat = ncread(fullfile(C, '2020_07_01.nc'), 'lat');
Lon = ncread(fullfile(C, '2020_07_01.nc'), 'lon');
Depth = ncread(fullfile(A, '2020_07_01.nc'), 'depth');

% Calculate pressure (dbar)
P = rho * Depth * g / 10^4;

% Target Location
target_lat = 22.76;
target_lon = -158.07167;

% Find the nearest latitude and longitude indices
[~, lat_idx] = min(abs(Lat - target_lat));  % Nearest latitude index
[~, lon_idx] = min(abs(Lon - target_lon));  % Nearest longitude index

% Display the selected grid point for verification
disp(['Nearest Latitude: ', num2str(Lat(lat_idx))]);
disp(['Nearest Longitude: ', num2str(Lon(lon_idx))]);

% Initialize arrays
[m, n, k, time_count] = size(S_all);
hycom_bpga_profiles = struct();
hycom_pressure_profiles = struct();

% Define colors for plotting
colors = lines(time_count);

% Create a new figure
figure;

% Loop over each time step (hour)
for num = 1:time_count
    % Initialize storage for bpga at each depth
    layer_bpga_all = nan(k-1, 1);
    layer_pressure_all = nan(k-1, 1);
    
    % Use the nearest grid point for latitude and longitude
    x = lon_idx;
    y = lat_idx;

    % Extract salinity and temperature profiles at (x,y) and time num
    Salinity = squeeze(S_all(x, y, :, num));
    Temperature = squeeze(T_all(x, y, :, num));
    
    % Remove NaNs
    valid_idx = ~isnan(Salinity) & ~isnan(Temperature);
    Salinity = Salinity(valid_idx);
    Temperature = Temperature(valid_idx);
    Pressure = P(valid_idx);

    % Skip if no valid data
    if isempty(Salinity)
        continue;
    end

    % Define uniform 2-meter depth bins within the valid depth range
    depth_bins = 0:2:max(Pressure);

    % Calculate geopotential anomaly (ga) for each layer
    try
        ga = sw_gpan(Salinity, Temperature, Pressure); % Geopotential anomaly

        % Calculate steric height difference between consecutive layers
        layer_ga_change = diff(ga);          % Difference in GA between layers
        layer_pressure_diff = diff(Pressure); % Depth difference between layers

        % Normalize GA change over each depth interval
        ga_change_per_meter = layer_ga_change  ./ layer_pressure_diff;

        % Map GA change into uniform 2m depth bins
        ga_per_bin = zeros(size(depth_bins));

        for i = 1:length(layer_pressure_diff)
            % Depth range for this layer
            depth_start = Pressure(i);
            depth_end = Pressure(i+1);

            % Find 2m bins within this range
            bin_indices = find(depth_bins >= depth_start & depth_bins < depth_end);

            % Assign the normalized GA change to these bins
            ga_per_bin(bin_indices) = ga_change_per_meter(i);
        end
        bpga_dbar = rho * ga_per_bin / 10000;
        % Plot the layer-by-layer bpga profile
        %plot(2*ga_change_per_meter*rho/10000,-Depth(2:end-1), 'LineWidth', 1.5, 'Color', colors(num, :));
        plot(2*ga_change_per_meter*rho/10000,-Depth(2:end-1), 'LineWidth', 1.5, 'Color', 'b');
        hold on;

        % Save bpga and pressure profiles
        date_str = ['D_', num2str(num)];
        hycom_bpga_profiles.(date_str) = 2*ga_change_per_meter*rho/10000;
        hycom_pressure_profiles.(date_str) = -Depth(2:end-1);

    catch ME
        warning(['Failed to calculate bpga for time step ', num2str(num), ': ', ME.message]);
        continue;
    end
end

% Customize the plot
xlabel('Layer-by-Layer Geopotential Anomaly Change (bpga) (dbar)');
ylabel('Depth (m)');
title('Model Output: Layer-by-Layer bpga Change in 2m Intervals');
grid on;
set(gca, 'FontSize', 18); % Set font size for axes

% Add legend for all time steps
legend_entries = {};
for num = 1:time_count
    legend_entries{end+1} = ['Hour ', num2str(num)];
end
%legend(legend_entries, 'Location', 'eastoutside', 'FontSize', 10);
ylim([-5000 0])

hold off;
%%
save('hycom_bpga_profiles_with_pressure.mat', 'hycom_bpga_profiles', 'hycom_pressure_profiles', 'time_all');
disp('Layer-by-layer bpga profiles and pressures have been saved.');
