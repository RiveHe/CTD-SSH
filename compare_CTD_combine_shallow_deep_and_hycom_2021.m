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
%%
%% Shallow casts only - Pressure < 500
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

% Initialize results table with an additional 'Source' column
results_table_shallow_500 = table('Size', [0 3], 'VariableTypes', {'datetime', 'double', 'string'}, ...
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

    % Remove rows with NaN values and filter pressure < 500
    valid_idx = ~isnan(Pressure) & ~isnan(Salinity) & ~isnan(Temperature) & Pressure <= 500;
    Pressure = Pressure(valid_idx);
    Salinity = Salinity(valid_idx);
    Temperature = Temperature(valid_idx);

    % Skip if no valid data is available for the current date
    if isempty(Pressure)
        disp(['No valid data for ', datestr(current_date), ' with pressure < 500, skipping.']);
        continue;
    end

    % Determine cast type
    cast_type = "Shallow Cast";

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
    results_table_shallow_500 = [results_table_shallow_500; new_row];
end

% Save results to a CSV file
%writetable(results_table_shallow_500, 'bpga_CTD_shallow_casts_pressure_lt_500_2021.csv');
disp('Results saved to bpga_CTD_shallow_casts_pressure_lt_500_2021.csv.');

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


%% shallow + hycom data

addpath('H:\PhD work\Miguel\seawater_ver3_3.1')
rho = 1040; %(kg/m3) 
g = 9.8; %(m/s2) 
% Load CTD data
ctd_file = 'combined_shallow_ctd_hycom_profiles_2021.csv';
ctd_data = readtable(ctd_file);

% Convert CTD DateTime to MATLAB datetime format
ctd_data.DateTime = datetime(ctd_data.DateTime, 'InputFormat', 'M/d/yyyy');

% Get unique CTD dates
unique_ctd_dates = unique(ctd_data.DateTime);
%% Initialize results table with an additional 'Source' column
results_table_shallow_hycom = table('Size', [0 3], 'VariableTypes', {'datetime', 'double', 'string'}, ...
                      'VariableNames', {'DateTime', 'bpga_combine', 'CastType'});
all_bpga_values = [];
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
        all_bpga_values = [all_bpga_values; bpga_combine]; % Collect values for outlier detection
    catch ME
        warning(['Failed to calculate geopotential anomaly for ', datestr(current_date), ': ', ME.message]);
        continue;
    end

    % Add the result to the results table
    new_row = table(current_date, bpga_combine, cast_type, 'VariableNames', {'DateTime', 'bpga_combine', 'CastType'});
    results_table_shallow_hycom = [results_table_shallow_hycom; new_row];
end

% Save results to a CSV file
writetable(results_table_shallow_hycom, 'bpga_CTD_shallow_and_hycom_2021.csv');
disp('Results saved to bpga_CTDonly_results_2020_withtag.csv.');
%%

% Remove outliers from bpga_combine using the IQR method
if ~isempty(all_bpga_values)
    Q1 = prctile(all_bpga_values, 25);
    Q3 = prctile(all_bpga_values, 75);
    IQR = Q3 - Q1;

    % Define the acceptable range
    lower_bound = Q1 - 1.5 * IQR;
    upper_bound = Q3 + 1.5 * IQR;

    % Filter the results table
    results_table_shallow_hycom = results_table_shallow_hycom(...
        results_table_shallow_hycom.bpga_combine >= lower_bound & ...
        results_table_shallow_hycom.bpga_combine <= upper_bound, :);
end

% Save the cleaned results to a CSV file
writetable(results_table_shallow_hycom, 'bpga_CTD_shallow_and_hycom_2021_cleaned.csv');
disp('Results saved to bpga_CTD_shallow_and_hycom_2021_cleaned.csv.');
%% Shallow 500m +hycom


addpath('H:\PhD work\Miguel\seawater_ver3_3.1')
rho = 1040; %(kg/m3) 
g = 9.8; %(m/s2) 
% Load CTD data
ctd_file = 'combined_ctd_hycom_depth_2021_500.csv';
ctd_data = readtable(ctd_file);

% Convert CTD DateTime to MATLAB datetime format
ctd_data.DateTime = datetime(ctd_data.DateTime, 'InputFormat', 'M/d/yyyy');

% Get unique CTD dates
unique_ctd_dates = unique(ctd_data.DateTime);
%% Initialize results table with an additional 'Source' column
results_table_shallow_hycom_500 = table('Size', [0 3], 'VariableTypes', {'datetime', 'double', 'string'}, ...
                      'VariableNames', {'DateTime', 'bpga_combine', 'CastType'});
all_bpga_values = [];
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
        all_bpga_values = [all_bpga_values; bpga_combine]; % Collect values for outlier detection
    catch ME
        warning(['Failed to calculate geopotential anomaly for ', datestr(current_date), ': ', ME.message]);
        continue;
    end

    % Add the result to the results table
    new_row = table(current_date, bpga_combine, cast_type, 'VariableNames', {'DateTime', 'bpga_combine', 'CastType'});
    results_table_shallow_hycom_500 = [results_table_shallow_hycom_500; new_row];
end

% Save results to a CSV file
%writetable(results_table_shallow_hycom, 'bpga_CTD_shallow_and_hycom_2021_500.csv');
disp('Results saved to bpga_CTDonly_results_2020_withtag.csv.');
%%

% Remove outliers from bpga_combine using the IQR method
if ~isempty(all_bpga_values)
    Q1 = prctile(all_bpga_values, 25);
    Q3 = prctile(all_bpga_values, 75);
    IQR = Q3 - Q1;

    % Define the acceptable range
    lower_bound = Q1 - 1.5 * IQR;
    upper_bound = Q3 + 1.5 * IQR;

    % Filter the results table
    results_table_shallow_hycom_500 = results_table_shallow_hycom_500(...
        results_table_shallow_hycom_500.bpga_combine >= lower_bound & ...
        results_table_shallow_hycom_500.bpga_combine <= upper_bound, :);
end

% Save the cleaned results to a CSV file
%writetable(results_table_shallow_hycom_500, 'bpga_CTD_shallow_and_hycom_2021_cleaned.csv');
disp('Results saved to bpga_CTD_shallow_and_hycom_2021_cleaned.csv.');
%%
% Calculate demeaned data
demeaned_data_shallow = results_table_shallow.bpga_combine - nanmean(results_table_shallow.bpga_combine);
demeaned_data_shallow_500 =  results_table_shallow_500.bpga_combine - nanmean(results_table_shallow_500.bpga_combine);
demeaned_data_deep = results_table.bpga_combine - nanmean(results_table.bpga_combine);
demeaned_data_shallow_deep = results_table_shallow_deep.bpga_combine - nanmean(results_table_shallow_deep.bpga_combine);
demeaned_data_shallow_hycom = results_table_shallow_hycom.bpga_combine - nanmean(results_table_shallow_hycom.bpga_combine);
demeaned_data_shallow_hycom_500 = results_table_shallow_hycom_500.bpga_combine - nanmean(results_table_shallow_hycom_500.bpga_combine);
% Find the common DateTime across all three datasets
common_dates_12 = intersect(results_table_shallow.DateTime, results_table.DateTime);  % Shallow & Deep
common_dates_123 = intersect(common_dates_12, results_table_shallow_deep.DateTime); 
common_dates_1234 = intersect(common_dates_123, results_table_shallow_hycom_500.DateTime);
common_dates_12345 = intersect(common_dates_1234,results_table_shallow_500.DateTime);
common_dates_all = intersect(common_dates_12345, results_table_shallow_hycom.DateTime); 


% Shallow, Deep & Shallow-Deep

% Find the indices for common dates in each dataset
[~, idx_shallow] = ismember(common_dates_all, results_table_shallow.DateTime);
[~, idx_deep] = ismember(common_dates_all, results_table.DateTime);
[~, idx_shallow_deep] = ismember(common_dates_all, results_table_shallow_deep.DateTime);
[~, idx_shallow_deep_hycom] = ismember(common_dates_all, results_table_shallow_hycom.DateTime);
[~, idx_shallow_deep_hycom_500] = ismember(common_dates_all, results_table_shallow_hycom_500.DateTime);
[~, idx_shallow_500] = ismember(common_dates_all, results_table_shallow_500.DateTime);

% Filter the data based on common dates
common_shallow_data = demeaned_data_shallow(idx_shallow);
common_deep_data = demeaned_data_deep(idx_deep);
common_shallow_deep_data = demeaned_data_shallow_deep(idx_shallow_deep);
common_shallow_hycom_data = demeaned_data_shallow_hycom(idx_shallow_deep_hycom);
common_shallow_hycom_data_500 = demeaned_data_shallow_hycom_500(idx_shallow_deep_hycom_500);
common_shallow_data_500 = demeaned_data_shallow_500(idx_shallow_500);

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

    scatter(common_dates_all(i), common_shallow_data_500(i) / 10000, ...
            60, 'pentagram', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', 'k', ...
            'DisplayName', 'Shallow Cast 500m');

    scatter(common_dates_all(i), common_shallow_hycom_data_500(i) / 10000, ...
            60, '+', 'MarkerFaceColor', colors(i,:),'MarkerEdgeColor', colors(i,:), 'LineWidth', 1.5,...
            'DisplayName', 'Shallow Cast 500m plus hycom');
    
    % Plot Shallow + Deep Cast data
    scatter(common_dates_all(i), common_shallow_deep_data(i) / 10000, ...
            60, 's', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', 'k', ...
            'DisplayName', 'Shallow + Deep Cast');

    % Plot Shallow + hycom data
    scatter(common_dates_all(i), common_shallow_hycom_data(i) / 10000, ...
            60, 'diamond', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', 'k', ...
            'DisplayName', 'Shallow + Deep Cast');

     % Plot Deep Cast data
    scatter(common_dates_all(i), common_deep_data(i) / 10000, ...
            60, '^', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', 'k', ...
            'DisplayName', 'Deep Cast');
    
end
plot(time_all, filtered_eff_bpga_nearest / 100, 'r', 'LineWidth', 1.5, ...
     'DisplayName', 'HYCOM Steric Height Pressure');

hold off;

% Add title, labels, and legend
title('Steric Height Pressure Combine Results Over Shared Time', 'FontSize', 14);
xlabel('Date', 'FontSize', 12);
ylabel('Pressure Anomaly (dbar)', 'FontSize', 12);
grid on;

% Combine legend entries for shallow, deep, and shallow+deep casts
legend({'Shallow Cast 1000m', 'Shallow Cast 500m', 'Shallow Cast 500m + HYCOM', ...
         'Shallow 1000m + Copy Deep Cast', 'Shallow 1000m + HYCOM', 'Deep Cast',...
        'HYCOM Steric Height Pressure'}, 'Location', 'best');


% Format the x-axis for better datetime display
datetick('x', 'mm-dd', 'keepticks');
%ylim([-0.5 0.5]);
xtickangle(0);  % Rotate x-axis labels for better readability
ax = gca; 
ax.FontSize = 16;