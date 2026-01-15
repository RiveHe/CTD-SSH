% Define the file path for the CSV file containing SSH Pressure Anomaly data
%https://data.marine.copernicus.eu/products
ssh_file_path = 'F:\aloha\ALOHA Cabled Observatory Database\ACO\ssh\2020\ssh_pressure_anomaly_series_2020.csv';

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

%% Add the required library for seawater calculations (replace with your library path)
addpath('F:\PhD work\Miguel\seawater_ver3_3.1')

% Initialize parameters
rho = 1040; % (kg/m^3)
g = 9.8; % (m/s^2)

% Base directory where the folders are located
base_dir = 'F:\aloha\ALOHA Cabled Observatory Database\ACO\hot\original';

% Output CSV file path
output_csv = 'F:\aloha\ALOHA Cabled Observatory Database\ACO\hot\bpga_ctd_results_2020_original.csv';

% Create an empty table to store results
results_table = table('Size', [0 2], 'VariableTypes', {'string', 'double'}, ...
                      'VariableNames', {'FileName', 'bpga_ctd'});

%% Loop through folders (hot-318 to hot-325) and files (h318a0201 to h325a0222)
for folder_num = 326:334
    folder_name = sprintf('hot-%d', folder_num);  % Construct the folder name
    folder_path = fullfile(base_dir, folder_name);  % Full path to the folder

    for file_num = 201:222
        file_name = sprintf('h%da0%03d.ctd', folder_num, file_num);  % Correct file naming
        file_path = fullfile(folder_path, file_name);  % Full path to the file

        % Check if the file exists
        if exist(file_path, 'file')
            disp(['Processing ', file_name]);

            % Define import options for reading the fixed-width text file
            opts = fixedWidthImportOptions('NumVariables', 8);
            opts.DataLines = [7 Inf];  % Skip the first 6 rows (header)
            opts.VariableWidths = [8, 8, 9, 8, 8, 8, 8, 8];  % Column widths
            opts.VariableNames = {'Pressure', 'Temperature', 'Salinity', 'Oxygen', ...
                                  'Transmission_Nitrate', 'Chloropigments', 'Observations', 'Quality'};

            % Read the CTD data from the file
            ctd_data = readtable(file_path, opts);

            % Assuming pressure, temperature, and salinity columns are in the files
            Pressure = ctd_data.Pressure; % Pressure in dbar
            Salinity = ctd_data.Salinity; % Salinity in PSU (PSS-78)
            Temperature = ctd_data.Temperature; % Temperature in degrees C (ITS-90)

            % Convert from cell arrays to numeric arrays, handling non-numeric cases
            if iscell(Pressure)
                Pressure = cellfun(@(x) str2double(x), Pressure, 'UniformOutput', true);
            end
            if iscell(Salinity)
                Salinity = cellfun(@(x) str2double(x), Salinity, 'UniformOutput', true);
            end
            if iscell(Temperature)
                Temperature = cellfun(@(x) str2double(x), Temperature, 'UniformOutput', true);
            end

            % Handle cases where the data may still have NaN values due to non-numeric cells
            Pressure(isnan(Pressure)) = NaN;
            Salinity(isnan(Salinity)) = NaN;
            Temperature(isnan(Temperature)) = NaN;

            % Now S_all, T_all, and P_all contain the data for the file
            time_count = length(Pressure);  % Set the time_count as the length of the Salinity data array

            % Initialize geopotential anomaly array
            bpga_ctd = nan;  % Single value for geopotential anomaly at the deepest level

            % Loop through data points and calculate geopotential anomaly (ga)
            if ~isempty(Pressure) && ~isnan(Salinity(1)) && ~isnan(Temperature(1)) && ~isnan(Pressure(1))
                try
                    ga = sw_gpan(Salinity(1:time_count), Temperature(1:time_count), Pressure(1:time_count));
                    bpga_ctd = rho * g * ga(end) / g;  % Calculate geopotential anomaly at the deepest level
                catch ME
                    warning(['Failed to calculate ga for ', file_name, ': ', ME.message]);
                    continue;
                end
            end

            % Store the result in the table
            new_row = table({file_name}, bpga_ctd, 'VariableNames', {'FileName', 'bpga_ctd'});
            results_table = [results_table; new_row];
        else
            disp(['File not found: ', file_name, ', skipping.']);
        end
    end
end

% Save the results to a CSV file
writetable(results_table, output_csv);
disp(['Results saved to ', output_csv]);
%%

% Define the file path for the CSV file
csv_file_path = 'F:\aloha\ALOHA Cabled Observatory Database\ACO\hot\bpga_ctd_results_2021.csv';
% Read the CSV file into a table
data_ctd = readtable(csv_file_path);

% Convert 'Date' column to datetime if it's not already
if iscell(data_ctd.Date) || ischar(data_ctd.Date)
    data_ctd.Date = datetime(data_ctd.Date, 'InputFormat', 'M/d/yyyy');  % Adjust format based on your data
end

% Convert 'Time' column to datetime if it's not already, using the format 'HH:mm'
if iscell(data_ctd.Time) || ischar(data_ctd.Time)
    data_ctd.Time = datetime(data_ctd.Time, 'InputFormat', 'HH:mm');
end

% Combine 'Date' and 'Time' into a single 'DateTime' column
data_ctd.DateTime = data_ctd.Date + timeofday(data_ctd.Time);  % Combine Date and Time
%%
bpga_ctd2 = nan(size(data_ctd.bpga_ctd)); 

bpga_ctd2=data_ctd.bpga_ctd(1:end)-nanmean(data_ctd.bpga_ctd(1:end));
output_table = table(data_ctd.DateTime, bpga_ctd2, 'VariableNames', {'DateTime', 'bpga_ctd2'});

%Define the output file path for saving the results
output_csv_file_path = 'F:\aloha\ALOHA Cabled Observatory Database\ACO\hot\bpga_ctd_results_with_demeaned_2021.csv';

%Save the table to the new CSV file
writetable(output_table, output_csv_file_path);

%Display a message indicating the file has been saved
disp(['Results saved to ', output_csv_file_path]);
%%

% Plot bpga_ctd vs DateTime
figure;
scatter(data_ctd.Time, bpga_ctd2/10000, 'b', 'LineWidth', 1.5);
title('Geopotential Anomaly (bpga\_ctd) over Time');
xlabel('Time');
ylabel('Geopotential Anomaly (bpga\_ctd) (m^2/s^2)');
grid on;

% Format the x-axis as year-month-day
datetick('x', 'yyyy-mm-dd HH:MM', 'keeplimits');

% Rotate the x-axis labels for better readability
xtickangle(45);
