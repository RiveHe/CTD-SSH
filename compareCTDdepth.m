%% Define CTD Files for Comparison
ctd_files = {'shallow_casts_2021.csv', 'deep_casts_2021.csv', 'combined_shallow_ctd_hycom_profiles_2021.csv'};
ctd_labels = {'Shallow CTD', 'Deep CTD', 'Shallow + HYCOM'};

% Define colors for plotting
ctd_colors = {'r', 'g', 'm'}; % Red, Green, Magenta

% Define the target date range
start_date = datetime(2021, 1, 1);
end_date = datetime(2021, 12, 31);
all_dates = (start_date:end_date)'; % Generate a full-year date range

%% Load SSH (Altimetry) Data
ssh_file_path = 'H:\aloha\ALOHA Cabled Observatory Database\ACO\ssh\2021\ssh_pressure_anomaly_series_2021.csv';
data_ssh = readtable(ssh_file_path);
data_ssh.Time = datetime(data_ssh.Time, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSSSSS');

%% Load BPR Data
%% deal with bpr data and outliers
filePath = 'H:\aloha\ALOHA Cabled Observatory Database\ACO\prs_bsp2\prsbsp2_pressure_downsample_detide_filtered_detrended_before_2021_06_02.csv';

opts = detectImportOptions(filePath, 'Delimiter', ',', 'NumHeaderLines', 0);
data = readtable(filePath, opts);


filePath = 'H:\aloha\ALOHA Cabled Observatory Database\ACO\prs_bsp2\prsbsp2_pressure_downsample_detide_filtered_detrended_after_2021_06_02.csv';

opts = detectImportOptions(filePath, 'Delimiter', ',', 'NumHeaderLines', 0);
data1 = readtable(filePath, opts);

% Extract data
time_series = data1.Time_UTC;
pressure_series = data1.Filtered_t_tide_detided_detrended;

% Compute Interquartile Range (IQR)
Q1 = quantile(pressure_series, 0.25); % 25th percentile
Q3 = quantile(pressure_series, 0.75); % 75th percentile
IQR_value = Q3 - Q1;

% Define outlier thresholds
lower_bound = Q1 - 1.5 * IQR_value;
upper_bound = Q3 + 1.5 * IQR_value;

% Remove outliers
valid_idx = (pressure_series >= lower_bound) & (pressure_series <= upper_bound);
filtered_time = time_series(valid_idx);
filtered_pressure = pressure_series(valid_idx);

% Merge timestamps
combined_BPR_time = [data.Time_UTC; filtered_time];

% Merge pressure data
combined_BPR_pressure = [data.Filtered_t_tide_detided_detrended; filtered_pressure];

% Sort to ensure chronological order
[combined_BPR_time, sort_idx] = sort(combined_BPR_time);
combined_BPR_pressure = combined_BPR_pressure(sort_idx);

%% Initialize Storage for Results
results_ctd = cell(length(ctd_files), 1);
trend_ctd = zeros(length(ctd_files), 1);

%% Process Each CTD File
for idx = 1:length(ctd_files)
    ctd_file = ctd_files{idx};
    final_ctd_data = readtable(ctd_file);

    % Get unique pressure levels
    unique_pressures = unique(final_ctd_data.Pressure);

    % Initialize Interpolated Table
    interpolated_ctd_data = table('Size', [0 5], ...
        'VariableTypes', {'datetime', 'double', 'double', 'double', 'cell'}, ...
        'VariableNames', {'DateTime', 'Pressure', 'Salinity', 'Temperature', 'FileName'});

    % Interpolate across time for each pressure level
    for j = 1:length(unique_pressures)
        pressure_level = unique_pressures(j);
        pressure_data = final_ctd_data(final_ctd_data.Pressure == pressure_level, :);

        if height(pressure_data) < 2
            continue;
        end

        % Use cubic interpolation
        interpolated_salinity = interp1(pressure_data.DateTime, pressure_data.Salinity, all_dates, 'linear', 'extrap');
        interpolated_temperature = interp1(pressure_data.DateTime, pressure_data.Temperature, all_dates, 'linear', 'extrap');

        % Store Interpolated Data
        interpolated_table = table(all_dates, ...
            repmat(pressure_level, length(all_dates), 1), ...
            interpolated_salinity, interpolated_temperature, ...
            repmat({''}, length(all_dates), 1), ...
            'VariableNames', {'DateTime', 'Pressure', 'Salinity', 'Temperature', 'FileName'});

        interpolated_ctd_data = [interpolated_ctd_data; interpolated_table];
    end

    % Sort Data
    interpolated_ctd_data = sortrows(interpolated_ctd_data, {'DateTime', 'Pressure'});

    % Save to File
    %writetable(interpolated_ctd_data, ['interpolated_', ctd_file]);

    % Compute Bottom Pressure Anomaly
    addpath('H:\PhD work\Miguel\seawater_ver3_3.1');
    rho = 1040; % Seawater density (kg/m^3)
    g = 9.8; % Gravitational acceleration (m/s^2)

    % Extract unique dates
    unique_dates = unique(interpolated_ctd_data.DateTime);
    results_table = table('Size', [0 2], 'VariableTypes', {'datetime', 'double'}, ...
        'VariableNames', {'DateTime', 'bpga_ctd'});

    for i = 1:length(unique_dates)
        current_date = unique_dates(i);
        daily_data = interpolated_ctd_data(interpolated_ctd_data.DateTime == current_date, :);

        % Extract relevant data
        Pressure = daily_data.Pressure;
        Salinity = daily_data.Salinity;
        Temperature = daily_data.Temperature;

        % Remove NaNs
        valid_idx = ~isnan(Pressure) & ~isnan(Salinity) & ~isnan(Temperature);
        Pressure = Pressure(valid_idx);
        Salinity = Salinity(valid_idx);
        Temperature = Temperature(valid_idx);

        if isempty(Pressure)
            continue;
        end

        % Compute Geopotential Anomaly and Bottom Pressure Anomaly
        try
            ga = sw_gpan(Salinity, Temperature, Pressure);
            bpga_ctd = rho * g * ga(end) / g;
        catch
            continue;
        end

        % Store results
        new_row = table(current_date, bpga_ctd, 'VariableNames', {'DateTime', 'bpga_ctd'});
        results_table = [results_table; new_row];
    end

    % Apply Low-Pass Butterworth Filter
    [b, a] = butter(4, 0.05, 'low');
    results_table.bpga_ctd = filtfilt(b, a, results_table.bpga_ctd);

    % Save results
    %writetable(results_table, ['bpga_results_', ctd_file]);
    results_ctd{idx} = results_table;
end
%%
%% Define Time Limits for Trend Analysis
start_date = datetime('2021-01-17');
end_date = datetime('2021-12-01');

%% Filter Data to the Specified Time Period

% Filter SSH Data
ssh_mask = (data_ssh.Time >= start_date) & (data_ssh.Time <= end_date);
filtered_ssh_time = data_ssh.Time(ssh_mask);
filtered_ssh_pressure = data_ssh.ssh_Pressure_Anomaly(ssh_mask);

% Initialize storage for SSH-CTD adjusted pressure anomaly
filtered_ssh_ctd_pressure = cell(length(ctd_files), 1);
demeanedresult = cell(length(ctd_files), 1); % Store demeaned results for each CTD file

% Compute Demeaned BPGA for Each CTD File
for idx = 1:length(ctd_files)
    demeanedresult{idx} = results_ctd{idx}.bpga_ctd - nanmean(results_ctd{idx}.bpga_ctd);
    
    % Compute SSH-CTD adjusted pressure anomaly
    adjusted_ssh_ctd = data_ssh.ssh_Pressure_Anomaly - demeanedresult{idx} / 10000;
    
    % Apply time filtering
    filtered_ssh_ctd_pressure{idx} = adjusted_ssh_ctd(ssh_mask);
end

% Filter BPR Data
bpr_mask = (combined_BPR_time >= start_date) & (combined_BPR_time <= end_date);
filtered_bpr_time = combined_BPR_time(bpr_mask);
filtered_bpr_pressure = combined_BPR_pressure(bpr_mask);

% Filter HYCOM Data
hycom_mask = (time_all >= start_date) & (time_all <= end_date);
filtered_hycom_time = time_all(hycom_mask);
filtered_hycom_pressure = filtered_eff_nearest(hycom_mask) / 100;


%% Plot Comparison of All Data: SSH, BPR, and CTD
figure;
hold on;

% Plot SSH (Altimetry) Data
plot(filtered_ssh_time, filtered_ssh_pressure, '--b', 'LineWidth', 1.5);

% Plot BPR Data
plot(filtered_bpr_time, filtered_bpr_pressure, 'k', 'LineWidth', 1.5);
plot(filtered_hycom_time, filtered_hycom_pressure, 'r', 'LineWidth', 1.5);

custom_ctd_colors = [ 
    0 1 1;   % Red for Shallow CTD
    1 0 1;   % Green for Deep CTD
    0 0 1;   % Blue for Combined Shallow + HYCOM
];

% Loop Through Each CTD File
for idx = 1:length(ctd_files)
     ctd_data = readtable(ctd_files{idx});  % Load the raw data
    
    % Extract the original CTD timestamps (before interpolation)
    ctd_times = unique(ctd_data.DateTime); % Ensure only unique timestamps
    % Load the original CTD file (before interpolation)
  
    % valid_ctd_mask = (ctd_times >= start_date) & (ctd_times <= end_date);
    % filtered_ctd_times = ctd_times(valid_ctd_mask); % Apply filter

    % Find the nearest SSH-CTD pressure values at these original times
    [~, nearest_idx] = min(abs(filtered_ssh_time - ctd_times'), [], 1); % Find closest indices
    sampled_ssh_ctd_pressure = filtered_ssh_ctd_pressure{idx}(nearest_idx); % Extract values
    % sampled_ssh_ctd_pressure =  sampled_ssh_ctd_pressure(valid_ctd_mask);

    % Plot smoothed CTD-derived SSH-adjusted pressures
    plot(filtered_ssh_time, filtered_ssh_ctd_pressure{idx}, 'Color', custom_ctd_colors(idx, :), 'LineWidth', 1.5);
end

for idx = 1:length(ctd_files)
      ctd_data = readtable(ctd_files{idx});  % Load the raw data
    
    % Extract the original CTD timestamps (before interpolation)
    ctd_times = unique(ctd_data.DateTime); % Ensure only unique timestamps
    [~, nearest_idx] = min(abs(filtered_ssh_time - ctd_times'), [], 1); % Find closest indices
    sampled_ssh_ctd_pressure = filtered_ssh_ctd_pressure{idx}(nearest_idx); % Extract values
     % Scatter original CTD data points on the smoothed line
    scatter(ctd_times, sampled_ssh_ctd_pressure, 30, custom_ctd_colors(idx, :), 'filled', 'MarkerEdgeColor', 'k');
end

text(1, text_y-0.25, trend_text, 'FontSize', 14, 'FontWeight', 'bold', 'BackgroundColor', 'w', 'EdgeColor', 'k');
xlabel('Time');
ylabel('Pressure (dbar)');
title('Comparison of SSH, BPR, and CTD-Derived Bottom Pressure Anomalies in Hawaii');
legend({'Altimetry Data', 'BPR Data', 'HYCOM Data','SSH-Shallow CTD', 'SSH-Deep CTD', 'SSH-Shallow + HYCOM','CTD Samples'}, 'Location', 'southeast');
ax = gca;
ax.FontSize = 18;
grid on;
hold off;
%% plot trend change the adjusted_ssh_ctd = data_ssh.ssh_Pressure_Anomaly - demeanedresult{idx} / 10000; to demeanedresult{idx} / 10000;
figure;
hold on;

% Convert datetime to numeric if necessary
if isdatetime(filtered_ssh_time)
    filtered_ssh_time = datenum(filtered_ssh_time);
end
if isdatetime(filtered_bpr_time)
    filtered_bpr_time = datenum(filtered_bpr_time);
end
if isdatetime(filtered_hycom_time)
    filtered_hycom_time = datenum(filtered_hycom_time);
end

% Ensure pressure values are numeric
filtered_ssh_pressure = double(filtered_ssh_pressure);
filtered_bpr_pressure = double(filtered_bpr_pressure);
filtered_hycom_pressure = double(filtered_hycom_pressure);

% Plot SSH (Altimetry) Data and Trend
plot(filtered_ssh_time, filtered_ssh_pressure, '--b', 'LineWidth', 1.5);
ssh_trend = polyfit(filtered_ssh_time, filtered_ssh_pressure, 1);
plot(filtered_ssh_time, polyval(ssh_trend, filtered_ssh_time), '--b', 'LineWidth', 2);

% Plot BPR Data and Trend
plot(filtered_bpr_time, filtered_bpr_pressure, 'k', 'LineWidth', 1.5);
bpr_trend = polyfit(filtered_bpr_time, filtered_bpr_pressure, 1);
plot(filtered_bpr_time, polyval(bpr_trend, filtered_bpr_time), '--k', 'LineWidth', 2);

% Plot HYCOM Data and Trend
plot(filtered_hycom_time, filtered_hycom_pressure, 'r', 'LineWidth', 1.5);

% Define Custom Colors for CTD Data
custom_ctd_colors = [ 
    0 1 1;   % Cyan for Shallow CTD
    1 0 1;   % Magenta for Deep CTD
    0 0 1;   % Blue for Combined Shallow + HYCOM
];

% Loop Through Each CTD File for Smoothed Data and Trend Computation
for idx = 1:length(ctd_files)
    ctd_data = readtable(ctd_files{idx});  
    ctd_times = unique(ctd_data.DateTime);  

    % Ensure CTD time is numeric
    if isdatetime(ctd_times)
        ctd_times = datenum(ctd_times);
    end

    % Ensure CTD pressure is numeric
    filtered_ssh_ctd_pressure{idx} = double(filtered_ssh_ctd_pressure{idx});

    % Find the nearest SSH-CTD pressure values at original times
    [~, nearest_idx] = min(abs(filtered_ssh_time - ctd_times'), [], 1);
    sampled_ssh_ctd_pressure = filtered_ssh_ctd_pressure{idx}(nearest_idx);

    % Plot smoothed CTD-derived SSH-adjusted pressures
    plot(filtered_ssh_time, filtered_ssh_ctd_pressure{idx}, 'Color', custom_ctd_colors(idx, :), 'LineWidth', 1.5);

    % Compute and plot trend for CTD data
    ctd_trend = polyfit(filtered_ssh_time, filtered_ssh_ctd_pressure{idx}, 1);
    plot(filtered_ssh_time, polyval(ctd_trend, filtered_ssh_time), '--', 'Color', custom_ctd_colors(idx, :), 'LineWidth', 2);
    ctd_trends(idx) = ctd_trend(1) * 365;
end

% Scatter Original CTD Data Points
for idx = 1:length(ctd_files)
    ctd_data = readtable(ctd_files{idx});
    
    % Ensure time column is numeric
    if isdatetime(ctd_data.DateTime)
        ctd_times = datenum(ctd_data.DateTime);
    else
        ctd_times = double(ctd_data.DateTime);
    end
    
    % Extract unique timestamps
    ctd_times = unique(ctd_times);
    [~, nearest_idx] = min(abs(filtered_ssh_time - ctd_times'), [], 1);
    sampled_ssh_ctd_pressure = filtered_ssh_ctd_pressure{idx}(nearest_idx);
    
    % Scatter original CTD data points
    scatter(ctd_times, sampled_ssh_ctd_pressure, 30, custom_ctd_colors(idx, :), 'filled', 'MarkerEdgeColor', 'k');
end

% Convert Trends to dbar per Year and Display in Upper-Left Corner
ssh_trend_yr = ssh_trend(1) * 365; % Convert from dbar/day to dbar/year
bpr_trend_yr = bpr_trend(1) * 365;
%hycom_trend_yr = hycom_trend(1) * 365;

% Define position for text annotation
text_x = min(filtered_ssh_time) + (max(filtered_ssh_time) - min(filtered_ssh_time)) * 0.02;
text_y = max(filtered_ssh_pressure) - (max(filtered_ssh_pressure) - min(filtered_ssh_pressure)) * 0.05;

% Create formatted text with trend values
trend_text = sprintf(['Trends (dbar/year):\n' ...
    'SSH: %.3f\nBPR: %.3f\nShallow CTD: %.3f\nDeep CTD: %.3f\nShallow+HYCOM: %.3f'], ...
    ssh_trend_yr, bpr_trend_yr, ctd_trends(1), ctd_trends(2), ctd_trends(3));

% Display text on plot
text(text_x, text_y, trend_text, 'FontSize', 14, 'FontWeight', 'bold', 'BackgroundColor', 'w', 'EdgeColor', 'k');
% Labels and Formatting
xlabel('Time');
ylabel('Pressure (dbar)');
title('Comparison of SSH, BPR, and CTD-Derived Bottom Pressure Anomalies in Hawaii');
% legend({'Altimetry Data', 'Altimetry Trend', ...
%         'BPR Data', 'BPR Trend', ...
%         'HYCOM Data', 'HYCOM Trend', ...
%         'Shallow CTD', 'Shallow CTD Trend', ...
%         'Deep CTD', 'Deep CTD Trend', ...
%         'Shallow + HYCOM', 'Shallow + HYCOM Trend', ...
%         'CTD Samples'}, 'Location', 'southeast');

ax = gca;
ax.FontSize = 18;
grid on;
hold off;
