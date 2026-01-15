%% Add the required library for seawater calculations
addpath('H:\PhD work\Miguel\seawater_ver3_3.1')

% Initialize parameters
target_lat = 22.75; % Latitude of the target location
target_lon = -158.0062 + 360; % Longitude of the target location
rho = 1025; % Approximate density of seawater (kg/m³)
g = 9.8; % Gravitational acceleration (m/s²)
%output_csv = 'H:\aloha\ALOHA Cabled Observatory Database\ACO\hot\combined_ctd_hycom_results.csv';

% Define the time range
daily_time_range = datetime(2020, 1, 1):datetime(2020, 1, 31);

% Define the CTD base directory
base_dir = 'H:\aloha\ALOHA Cabled Observatory Database\ACO\hot';

% Create an empty array to store CTD data
ctd_data_all = [];

%% Process CTD Data
% Initialize a table to store the processed CTD data
ctd_salinity_and_temperature_2020 = table('Size', [0 4], ...
    'VariableTypes', {'double', 'double', 'double', 'string'}, ...
    'VariableNames', {'Pressure', 'Salinity', 'Temperature', 'FileName'});

% Process CTD Data
for folder_num = 326:334
    folder_name = sprintf('hot-%d', folder_num); % Construct the folder name
    folder_path = fullfile(base_dir, folder_name); % Full path to the folder

    for file_num = 201:222
        file_name = sprintf('h%da0%03d.ctd', folder_num, file_num); % Construct the file name
        file_path = fullfile(folder_path, file_name); % Full path to the file

        if exist(file_path, 'file')
            disp(['Processing ', file_name]);

            % Define import options for reading the CTD file
            opts = fixedWidthImportOptions('NumVariables', 8);
            opts.DataLines = [7 Inf]; % Skip the first 6 rows
            opts.VariableWidths = [8, 8, 9, 8, 8, 8, 8, 8]; % Column widths
            opts.VariableNames = {'Pressure', 'Temperature', 'Salinity', 'Oxygen', ...
                                  'Transmission_Nitrate', 'Chloropigments', 'Observations', 'Quality'};

            % Read the CTD data
            ctd_data = readtable(file_path, opts);

            % Convert columns to numeric arrays
            Pressure = str2double(string(ctd_data.Pressure));
            Salinity = str2double(string(ctd_data.Salinity));
            Temperature = str2double(string(ctd_data.Temperature));

            % Remove NaN values
            valid_idx = ~isnan(Pressure) & ~isnan(Salinity) & ~isnan(Temperature);
            Pressure = Pressure(valid_idx);
            Salinity = Salinity(valid_idx);
            Temperature = Temperature(valid_idx);

            % Append to the CTD data table if not empty
            if ~isempty(Pressure)
                new_data = table(Pressure, Salinity, Temperature, repmat(string(file_name), length(Pressure), 1), ...
                                 'VariableNames', {'Pressure', 'Salinity', 'Temperature', 'FileName'});
                ctd_salinity_and_temperature_2020 = [ctd_salinity_and_temperature_2020; new_data];
            end
        else
            disp(['File not found: ', file_name, ', skipping.']);
        end
    end
end

% Display a message indicating the process is complete
disp('CTD data processing complete.');


% Optional: Save the table to a CSV file for easier viewing
writetable(ctd_salinity_and_temperature_2020, 'ctd_salinity_and_temperature_2021_combine_deep_data.csv');
disp('CTD data saved to ctd_salinity_and_temperature.csv.');

%% cat HYCOM to CTD
% Directories for salinity, temperature, and SSH data
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

%% Downsample HYCOM to daily data 
% Define time data
unique_days = unique(dateshift(time_all, 'start', 'day')); % Unique daily timestamps

% Initialize arrays for daily downsampled data
daily_S_all = NaN(size(S_all, 1), size(S_all, 2), size(S_all, 3), length(unique_days));
daily_T_all = NaN(size(T_all, 1), size(T_all, 2), size(T_all, 3), length(unique_days));

% Loop through each unique day and compute the daily averages
for d = 1:length(unique_days)
    % Find the indices corresponding to the current day
    day_indices = find(dateshift(time_all, 'start', 'day') == unique_days(d));
    
    % Compute the mean for the current day across the time dimension
    daily_S_all(:, :, :, d) = mean(S_all(:, :, :, day_indices), 4, 'omitnan');
    daily_T_all(:, :, :, d) = mean(T_all(:, :, :, day_indices), 4, 'omitnan');
end

% Confirm results
disp(['Downsampled S_all size: ', mat2str(size(daily_S_all))]);
disp(['Downsampled T_all size: ', mat2str(size(daily_T_all))]);

% Extract unique daily timestamps
downsampled_time_all = unique(dateshift(time_all, 'start', 'day'));

% Confirm the size of downsampled time array
disp(['Downsampled time_all size: ', num2str(length(downsampled_time_all))]);
%% Find the target point in HYCOM data
Lat = ncread(fullfile(C, '2020_07_01.nc'), 'lat');
Lon = ncread(fullfile(C, '2020_07_01.nc'), 'lon');
Depth = ncread(fullfile(A, '2020_07_01.nc'), 'depth');
% Target location
target_lat = 22.75;  % Latitude of the target location
target_lon = -158.0062 + 360;  % Longitude of the target location

% Find the nearest latitude and longitude indices
[~, lat_idx] = min(abs(Lat - target_lat));
[~, lon_idx] = min(abs(Lon - target_lon));

% Extract the salinity and temperature data at the target location
target_S = squeeze(daily_S_all(lon_idx, lat_idx, :, :)); % [Depth x Time]
target_T = squeeze(daily_T_all(lon_idx, lat_idx, :, :)); % [Depth x Time]

% Confirm the size of the extracted data
disp(['Target S size: ', mat2str(size(target_S))]); % Should be [40 x 365]
disp(['Target T size: ', mat2str(size(target_T))]); % Should be [40 x 365]


%% Load CTD data
ctd_file = 'ctd_salinity_and_temperature_2020_time.csv';
ctd_data = readtable(ctd_file);

% Convert CTD DateTime to MATLAB datetime format
ctd_data.DateTime = datetime(ctd_data.DateTime, 'InputFormat', 'M/d/yyyy');

% Get unique CTD dates
unique_ctd_dates = unique(ctd_data.DateTime);

%%
% Initialize arrays for storing matched salinity and temperature
matched_S = [];
matched_T = [];
matched_dates = [];

% Loop through each unique CTD date
for i = 1:length(unique_ctd_dates)
    % Find the index of the matching date in downsampled_time_all
    match_idx = find(downsampled_time_all == unique_ctd_dates(i), 1);
    
    if ~isempty(match_idx)
        % If a match is found, extract the corresponding S and T data
        matched_S = [matched_S, target_S(:, match_idx)]; % Append salinity for the date
        matched_T = [matched_T, target_T(:, match_idx)]; % Append temperature for the date
        matched_dates = [matched_dates; unique_ctd_dates(i)]; % Keep track of matched dates
    end
end

% Display results
disp(['Matched dates count: ', num2str(length(matched_dates))]);
disp('Matched salinity (S) matrix size:');
disp(size(matched_S)); % [Depth x Matched Dates]
disp('Matched temperature (T) matrix size:');
disp(size(matched_T)); % [Depth x Matched Dates]
%% Ensure consistent column structure
ctd_data.Source = repmat("CTD", height(ctd_data), 1); % Add Source column to CTD data

%% Find maximum pressure for each date in CTD data
max_ctd_pressure_by_date = varfun(@max, ctd_data, 'InputVariables', 'Pressure', ...
    'GroupingVariables', 'DateTime');
max_ctd_pressure_by_date.Properties.VariableNames{'GroupCount'} = 'Count';

% Display the max pressure for each date
disp(max_ctd_pressure_by_date);





%% Loop through each date to compare and concatenate HYCOM data
final_ctd_data = ctd_data; % Initialize with original CTD data

for i = 1:height(max_ctd_pressure_by_date)
    % Current date and its max pressure
    current_date = max_ctd_pressure_by_date.DateTime(i);
    max_ctd_pressure = max_ctd_pressure_by_date.max_Pressure(i);
    
    % Find the index of the current date in downsampled_time_all
    hycom_idx = find(downsampled_time_all == current_date, 1);
    
    if isempty(hycom_idx)
        % Skip if no matching HYCOM data for the current date
        disp(['No HYCOM data found for ', datestr(current_date)]);
        continue;
    end
    
    % Get HYCOM data for the current date
    hycom_depth = Depth(:); % Convert HYCOM depth to column vector
    hycom_salinity = target_S(:, hycom_idx);
    hycom_temperature = target_T(:, hycom_idx);
    
    % Filter HYCOM data for depths greater than max CTD pressure
    valid_hycom_idx = hycom_depth > max_ctd_pressure;
    
    if any(valid_hycom_idx)
        % Create a table for the filtered HYCOM data
        hycom_table = table(... 
            repmat(current_date, sum(valid_hycom_idx), 1), ... % Repeat the date
            hycom_depth(valid_hycom_idx), ... % Depths greater than max pressure
            hycom_salinity(valid_hycom_idx), ... % Salinity at those depths
            hycom_temperature(valid_hycom_idx), ... % Temperature at those depths
            repmat("HYCOM", sum(valid_hycom_idx), 1), ... % Source identifier
            'VariableNames', {'DateTime', 'Pressure', 'Salinity', 'Temperature', 'Source'}); 
        
        % Add missing columns in hycom_table
        missing_columns = setdiff(final_ctd_data.Properties.VariableNames, hycom_table.Properties.VariableNames);
        for col = missing_columns
            if strcmp(col{1}, 'FileName')
                % Add 'FileName' column with empty cells for HYCOM data
                hycom_table.(col{1}) = repmat({''}, height(hycom_table), 1);
            else
                % Add other missing columns with NaN
                hycom_table.(col{1}) = repmat(NaN, height(hycom_table), 1);
            end
        end
        
        % Reorder columns to match final_ctd_data structure
        hycom_table = hycom_table(:, final_ctd_data.Properties.VariableNames);
        
        % Concatenate the HYCOM data with the CTD data
        final_ctd_data = [final_ctd_data; hycom_table];
    end
end

%% Sort the final data by date and pressure
final_ctd_data = sortrows(final_ctd_data, {'DateTime', 'Pressure'});

%% Save the combined data to a CSV file
writetable(final_ctd_data, 'combined_ctd_hycom_profiles_2020.csv');
disp('Combined CTD and HYCOM data saved to combined_ctd_hycom_profiles.csv.');

%% %% Define the target date range
start_date = datetime(2020, 1, 1);
end_date = datetime(2020, 12, 31);
all_dates = (start_date:end_date)'; % Generate a full-year date range

% Initialize the final interpolated table
interpolated_ctd_data = table('Size', [0 6], ...
    'VariableTypes', {'datetime', 'double', 'double', 'double', 'cell', 'string'}, ...
    'VariableNames', {'DateTime', 'Pressure', 'Salinity', 'Temperature', 'FileName', 'Source'});

%% Loop through each day in the target date range
for i = 1:length(all_dates)
    current_date = all_dates(i);

    % Filter data for the current date
    daily_data = final_ctd_data(final_ctd_data.DateTime == current_date, :);

    if ~isempty(daily_data)
        % If data exists for the current date, append it to the interpolated table
        interpolated_ctd_data = [interpolated_ctd_data; daily_data];
    else
        % If data is missing for the current date, interpolate using surrounding data
        % Get unique pressures from existing data
        unique_pressures = unique(final_ctd_data.Pressure);

        % Preallocate arrays for interpolation
        interpolated_salinity = NaN(size(unique_pressures));
        interpolated_temperature = NaN(size(unique_pressures));

        % Loop through each pressure level to interpolate salinity and temperature
        for j = 1:length(unique_pressures)
            pressure_level = unique_pressures(j);

            % Filter data for the current pressure level
            pressure_data = final_ctd_data(final_ctd_data.Pressure == pressure_level, :);

            if ~isempty(pressure_data)
                % Interpolate salinity and temperature for the current pressure level
                interpolated_salinity(j) = interp1(datenum(pressure_data.DateTime), ...
                    pressure_data.Salinity, datenum(current_date), 'linear', 'extrap');
                interpolated_temperature(j) = interp1(datenum(pressure_data.DateTime), ...
                    pressure_data.Temperature, datenum(current_date), 'linear', 'extrap');
            end
        end

        % Create a new table for the interpolated data
        interpolated_table = table(...
            repmat(current_date, length(unique_pressures), 1), ... % Date
            unique_pressures, ... % Pressure levels
            interpolated_salinity, ... % Interpolated salinity
            interpolated_temperature, ... % Interpolated temperature
            repmat({''}, length(unique_pressures), 1), ... % FileName
            repmat("Interpolated", length(unique_pressures), 1), ... % Source
            'VariableNames', {'DateTime', 'Pressure', 'Salinity', 'Temperature', 'FileName', 'Source'});

        % Append the interpolated data to the final table
        interpolated_ctd_data = [interpolated_ctd_data; interpolated_table];
    end
end

%% Sort the final data by date and pressure
interpolated_ctd_data = sortrows(interpolated_ctd_data, {'DateTime', 'Pressure'});

%% Save the interpolated data to a CSV file
writetable(interpolated_ctd_data, 'interpolated_ctd_hycom_profiles_2020.csv');
disp('Interpolated CTD and HYCOM data saved to interpolated_ctd_hycom_profiles.csv.');

%%
%% Load the interpolated CTD and HYCOM combined data
data_file = 'interpolated_ctd_hycom_profiles_2020.csv';
combined_data = readtable(data_file);

% Extract unique dates in the dataset
unique_dates = unique(combined_data.DateTime);

% Initialize results table
results_table = table('Size', [0 2], 'VariableTypes', {'datetime', 'double'}, ...
                      'VariableNames', {'DateTime', 'bpga_combine'});

% Add required seawater library path
addpath('H:\PhD work\Miguel\seawater_ver3_3.1');

% Physical constants
rho = 1040; % Seawater density (kg/m^3)
g = 9.8; % Gravitational acceleration (m/s^2)
%%
% Read combined_data table
data = readtable('interpolated_ctd_hycom_profiles_2020.csv');

% Extract unique sources
unique_sources = unique(data.Source);

% Define colors for each source
colors = lines(length(unique_sources)); % Generates distinct colors

% Create a figure
figure;
hold on;

% Loop through each unique source and plot
for i = 1:length(unique_sources)
    % Get data for the current source
    source_mask = strcmp(data.Source, unique_sources{i});
    source_dates = data.DateTime(source_mask);
    source_pressure = data.Temperature(source_mask);
    
    % Plot with a unique color
    scatter(source_dates, source_pressure, 2, colors(i, :), 'filled', 'DisplayName', unique_sources{i});
end

% Add title, labels, and legend
title('Combined Data - Temperature vs DateTime');
xlabel('DateTime');
ylabel('Temperature (c)');
legend('Location', 'best');
grid on;
hold off;

%% Initialize results table with an additional 'Source' column
results_table = table('Size', [0 3], 'VariableTypes', {'datetime', 'double', 'string'}, ...
                      'VariableNames', {'DateTime', 'bpga_combine', 'Source'});

% Loop through each unique date
for i = 1:length(unique_dates)
    % Filter the data for the current date
    current_date = unique_dates(i);
    daily_data = combined_data(combined_data.DateTime == current_date, :);

    % Extract pressure, salinity, temperature, and source for the current day
    Pressure = daily_data.Pressure;
    Salinity = daily_data.Salinity;
    Temperature = daily_data.Temperature;
    Source = daily_data.Source; % Extract source

    % Remove rows with NaN values
    valid_idx = ~isnan(Pressure) & ~isnan(Salinity) & ~isnan(Temperature);
    Pressure = Pressure(valid_idx);
    Salinity = Salinity(valid_idx);
    Temperature = Temperature(valid_idx);
    Source = Source(valid_idx);

    % Skip if no valid data is available for the current date
    if isempty(Pressure)
        disp(['No valid data for ', datestr(current_date), ', skipping.']);
        continue;
    end

    % Calculate geopotential anomaly (ga) and bottom pressure anomaly
    try
        ga = sw_gpan(Salinity, Temperature, Pressure); % Geopotential anomaly
        bpga_combine = rho * g * ga(end) / g; % Bottom pressure anomaly at deepest level
    catch ME
        warning(['Failed to calculate geopotential anomaly for ', datestr(current_date), ': ', ME.message]);
        continue;
    end

    % Add the result to the results table (using the first source tag for simplicity)
    new_row = table(current_date, bpga_combine, string(Source{1}), 'VariableNames', {'DateTime', 'bpga_combine', 'Source'});
    results_table = [results_table; new_row];
end

% Save results to a CSV file
writetable(results_table, 'bpga_combine_results_2020_withtag.csv');
disp('Results saved to results_with_tags.csv.');


%% Save the results to a CSV file
output_csv = 'bpga_combine_results_2020.csv';
writetable(results_table, output_csv);
disp(['Geopotential anomalies and bottom pressure results saved to ', output_csv]);
%%
%% Load the bpga_combine results
data_file = 'bpga_combine_results_2020_withtag.csv';
bpga_data = readtable(data_file);

%% Demean bpga_combine
bpga_combine_demeaned = bpga_data.bpga_combine - nanmean(bpga_data.bpga_combine);

%% Add Demeaned Values to the Table
bpga_data.bpga_combine_demeaned = bpga_combine_demeaned;

%% Save the Updated Table with Tags
output_file = 'bpga_combine_results_2020_withtag_demeaned.csv';
writetable(bpga_data, output_file);

disp(['Demeaned bpga_combine saved to ', output_file]);


%% Save the demeaned results to a table
output_table = table(bpga_data.DateTime, bpga_combine_demeaned, ...
                     'VariableNames', {'DateTime', 'bpga_combine_demeaned'});

% Define the output file path for saving the results
output_csv_file_path = 'bpga_combine_results_demeaned_2020.csv';

% Save the table to the new CSV file
writetable(output_table, output_csv_file_path);

% Display a message indicating the file has been saved
disp(['Demeaned bpga_combine results saved to ', output_csv_file_path]);

%% Plot bpga_combine_demeaned vs DateTime
% Group data by source
sources = unique(bpga_data.Source);

% Create a figure
figure;
hold on;

% Assign colors or markers for each source
colors = lines(length(sources)); % Generate distinct colors

for i = 1:length(sources)
    % Filter data by the current source
    current_source = sources{i};
    source_data = bpga_data(strcmp(bpga_data.Source, current_source), :);
    
    % Plot the data for the current source
    scatter(source_data.DateTime, source_data.bpga_combine / 10000, ...
        50, 'MarkerFaceColor', colors(i, :), 'DisplayName', current_source, ...
        'MarkerEdgeColor', 'k');
end

% Add labels and title
title('Demeaned Geopotential Anomaly (steric height) over Time by Source');
xlabel('Time');
ylabel('Geopotential Anomaly (bpga\_combine) (m^2/s^2)');
legend('Location', 'best'); % Add a legend
grid on;
hold off;

%%
% Define the file path for the CSV file containing SSH Pressure Anomaly data
%https://data.marine.copernicus.eu/products
ssh_file_path = 'H:\aloha\ALOHA Cabled Observatory Database\ACO\ssh\2020\ssh_pressure_anomaly_series_2020.csv';

% Read the CSV file into a table
data_ssh = readtable(ssh_file_path);

% Convert the 'Time' column to datetime, assuming it's in ISO 8601 format
data_ssh.Time = datetime(data_ssh.Time, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSSSSS');

% Plot the time series (Time vs. ssh_Pressure_Anomaly)
figure;
plot(data_ssh.Time, data_ssh.ssh_Pressure_Anomaly, 'b', 'LineWidth', 1.5);
title('SSH Pressure Anomaly over Time');
xlabel('Time');
ylabel('SSH Pressure Anomaly (dbar)');
%%
% Plot the time series (Time vs. ssh_Pressure_Anomaly)
sources = unique(results_table.Source);

figure;
subplot(211)
scatter(data_ssh.Time, data_ssh.ssh_Pressure_Anomaly-bpga_combine_demeaned / 10000);
hold on 
plot(data_ssh.Time, data_ssh.ssh_Pressure_Anomaly, 'b', 'LineWidth', 1.5);
plot(data_prs.Time_UTC,data_prs.Filtered_t_tide_detided_hourly,'r','LineWidth', 2)
legend('SSH - CTD pressure', 'SSH pressure', 'BPR data','Location','southeast')
title('SSH - CTD Pressure Anomaly over Time');
xlabel('Time');
ylabel('SSH - CTD Pressure Anomaly (dbar)');
subplot(212)
% Assign colors or markers for each source
colors = lines(length(sources)); % Generate distinct colors
hold on
for i = 1:length(sources)
    % Filter data by the current source
    current_source = sources{i};
    source_data = results_table(strcmp(results_table.Source, current_source), :);
    
    % Plot the data for the current source
    scatter(source_data.DateTime, source_data.bpga_combine / 10000, ...
        50, 'MarkerFaceColor', colors(i, :), 'DisplayName', current_source, ...
        'MarkerEdgeColor', 'k');
end
title('Demeaned Geopotential Anomaly (steric height) over Time by Source');
xlabel('Time');
ylabel('Geopotential Anomaly (bpga\_combine) (m^2/s^2)');
legend('Location', 'southeast'); % Add a legend
grid on;


%%
% Define a finer time grid for interpolation
fine_time_grid = linspace(1, length(data_ssh.Time), 10000); % 1000 points for smoothness

% Interpolate data onto the finer grid
interp_anomaly = interp1(1:length(data_ssh.Time), ...
    data_ssh.ssh_Pressure_Anomaly - bpga_combine_demeaned / 10000, ...
    fine_time_grid, 'spline');

% Plot the interpolated data
figure;
subplot(211)
plot(fine_time_grid, interp_anomaly, 'b', 'LineWidth', 1.5);
subplot(212)
plot(data_prs.Time_UTC,data_prs.Filtered_t_tide_detided_hourly,'r','LineWidth', 2)
hold on 
plot(data_ssh.Time, data_ssh.ssh_Pressure_Anomaly, 'b', 'LineWidth', 2);
title('SSH Pressure Anomaly over Time - Smoothed');
xlabel('Time');
ylabel('SSH Pressure Anomaly (dbar)');
grid on;
%% 
% Define the standard deviation and create a Gaussian kernel
sigma = 50; % Adjust this value as needed
kernel_size = 6 * sigma; % Kernel size (large enough to cover most weights)
x = -kernel_size:kernel_size;
gaussian_kernel = exp(-x.^2 / (2 * sigma^2));
gaussian_kernel = gaussian_kernel / sum(gaussian_kernel); % Normalize the kernel

% Apply Gaussian smoothing using convolution
smoothed_data = conv(interp_anomaly, gaussian_kernel, 'same');

% Plot the smoothed data
figure;
subplot(211)
plot(fine_time_grid, smoothed_data, 'b', 'LineWidth', 1.5);
subplot(212)
plot(data_prs.Filtered_t_tide_detided_hourly,'r','LineWidth', 2)
title('SSH Pressure Anomaly - Smoothed (Gaussian Filter - Manual)');
xlabel('Time');
ylabel('SSH Pressure Anomaly (dbar)');
grid on;
%%
% Plot the smoothed data
% Define fine_time_grid with corresponding datetime values
start_date = datetime(2020, 1, 1);
end_date = datetime(2020, 12, 31);
fine_time_grid = linspace(datenum(start_date), datenum(end_date), length(smoothed_data)); % Ensure the length matches smoothed_data
fine_time_grid = datetime(fine_time_grid, 'ConvertFrom', 'datenum'); % Convert to datetime

% Plot with annotated time
figure;

% Subplot 1: Smoothed Data
subplot(211);
plot(fine_time_grid, smoothed_data, 'b', 'LineWidth', 1.5);
title('Smoothed Pressure - Gaussian Filter');
xlabel('Time');
ylabel('SSH - CTD Pressure Anomaly (dbar)');
grid on;

% Subplot 2: Original Data
subplot(212);
plot(data_prs.Time_UTC, data_prs.Filtered_t_tide_detided_hourly, 'r', 'LineWidth', 2);
title('Detided Hourly Data');
xlabel('Time');
ylabel('BPR Pressure Anomaly (dbar)');
grid on;

% Adjust x-axis limits and format for clarity
xlim([start_date end_date]);
datetick('x', 'yyyy-mm', 'keeplimits'); % Optional: Format datetime tick labels
%% interpolate ctd steric pressure and plot as Matt requested


filePath = 'H:\aloha\ALOHA Cabled Observatory Database\ACO\prs_bsp4\prs_pressure_downsample_detide_filtered_detrended_2020.csv';

opts = detectImportOptions(filePath, 'Delimiter', ',', 'NumHeaderLines', 0);
data_prs = readtable(filePath, opts);

% Convert Time_UTC to datetime format if not already
if ~isdatetime(data_prs.Time_UTC)
    data_prs.Time_UTC = datetime(data_prs.Time_UTC, 'InputFormat', 'yyyy-MM-dd HH:mm:ss', 'TimeZone', 'UTC');
end
%% Define the file path for the CSV file
csv_file_path = "H:\aloha\ALOHA Cabled Observatory Database\ACO\hot\avg_CTD_bpga_2020.xlsx";

% Read the CSV file into a table
data_ctd_diff_avg = readtable(csv_file_path);

data_ctd_diff_avg.Date = datetime(data_ctd_diff_avg.DateTime, 'InputFormat', 'M/d/yyyy');
% Define the start and end date of the year
start_date = datetime(2020, 1, 1);
end_date = datetime(2020, 12, 31);

% Create a daily time series for the entire year
year_time_series = start_date:end_date;

% Ensure the input date and diff arrays are sorted by date
[data_ctd_diff_avg.Date, sortIdx] = sort(data_ctd_diff_avg.Date);
data_ctd_diff_avg.diff = data_ctd_diff_avg.diff(sortIdx);

% Interpolate the diff values to the complete time series
interpolated_diff = interp1(data_ctd_diff_avg.Date, data_ctd_diff_avg.diff, year_time_series, 'spline', 'extrap');

% Plot the original and interpolated data
figure;

% Scatter plot of original data
scatter(data_ctd_diff_avg.Date, data_ctd_diff_avg.diff, 'c', '^', 'LineWidth', 1.5);
hold on;

% Line plot of interpolated data
plot(year_time_series, interpolated_diff, 'b', 'LineWidth', 1.5);

% Add labels and legend
title('Interpolated Time Series (2020)');
xlabel('Date');
ylabel('Diff');
legend('Original Data', 'Interpolated Data', 'Location', 'best');
grid on;
%%
% Define the file path for the CSV file containing SSH Pressure Anomaly data
%https://data.marine.copernicus.eu/products
ssh_file_path = 'H:\aloha\ALOHA Cabled Observatory Database\ACO\ssh\2020\ssh_pressure_anomaly_series_2020.csv';

% Read the CSV file into a table
data_ssh = readtable(ssh_file_path);

% Convert the 'Time' column to datetime, assuming it's in ISO 8601 format
data_ssh.Time = datetime(data_ssh.Time, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSSSSS');

% Plot the time series (Time vs. ssh_Pressure_Anomaly)
figure;
plot(data_ssh.Time, data_ssh.ssh_Pressure_Anomaly, 'b', 'LineWidth', 1.5);
title('SSH Pressure Anomaly over Time');
xlabel('Time');
ylabel('SSH Pressure Anomaly (dbar)');
%%
figure;
plot(data_ssh.Time, interpolated_diff', 'b', 'LineWidth', 1.5);

xlim([datetime(2020, 1, 7), datetime(2020, 12, 20)]);
hold on
plot(data_ssh.Time, data_ssh.ssh_Pressure_Anomaly, 'r', 'LineWidth', 1.5);
plot(data_prs.Time_UTC,data_prs.Filtered_t_tide_detided_hourly,'r','LineWidth', 2)
scatter(year_time_series, interpolated_diff, 'b', 'LineWidth', 1.5);
scatter(data_ctd_diff_avg.Date, data_ctd_diff_avg.diff, 'c', '^', 'LineWidth', 1.5);
legend('SSH-CTD pressure','SSH pressure','BPR hourly data','Interplated points for steric pressure','CTD data points')
title('Steric pressure interpolation points');
xlabel('Time');
ylabel('Pressure Anomaly (dbar)');
