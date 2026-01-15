% Assuming your table is named 'data'
% Load your table if not already loaded
data = readtable('ctd_salinity_and_temperature_2021_combine_deep_data.csv'); % Uncomment if reading from a CSV

%% Remove rows where FileName is 'h3260201.ctd'
data(strcmp(data.FileName, 'h334a0210.ctd'), :) = [];

% Optional: Save the updated table to a CSV file
writetable(data, 'ctd_salinity_and_temperature_2021_combine_deep_data_filtered.csv');
data_filtered = readtable('ctd_salinity_and_temperature_2021_combine_deep_data_filtered.csv')
% Display confirmation
disp('Rows with FileName "h326a0201.ctd" have been removed.');
%%
% Assuming your table is named 'data'
% Define the DateTime you want to assign
new_datetime = datetime('2021-12-07', 'InputFormat', 'yyyy-MM-dd');

% Check if 'DateTime' column exists; if not, create it
if ~ismember('DateTime', data.Properties.VariableNames)
    data.DateTime = NaT(height(data), 1); % Initialize with NaT (Not-a-Time)
end

% Assign the new DateTime to rows where FileName is 'h3260201.ctd'
data.DateTime(strcmp(data.FileName, 'h334a0211.ctd')) = new_datetime;

% Optional: Save the updated table to a CSV file
writetable(data, 'ctd_salinity_and_temperature_2021_combine_deep_data_filtered.csv');

% Display confirmation
disp('DateTime has been added to rows with FileName "h3260201.ctd".');
%% check the date and filenames are right
% Load the data from the CSV file
data_edited = readtable('ctd_salinity_and_temperature_2021_combine_deep_data_filtered_time.csv');

% Identify rows where DateTime is missing or invalid (NaT)
invalid_datetime_idx = isnat(data_edited.DateTime);

% Extract unique FileNames associated with invalid/missing DateTime
invalid_filenames = unique(data_edited.FileName(invalid_datetime_idx));

% Display the unique filenames with invalid/missing DateTime
disp('FileNames with Invalid or Missing DateTime:');
for i = 1:length(invalid_filenames)
    fprintf('FileName: %s\n', invalid_filenames{i});
end


%%
% Extract unique combinations of DateTime and FileName
unique_pairs = unique(data_edited(:, {'DateTime', 'FileName'}));

% Display each unique DateTime with its corresponding FileName
disp('Unique DateTime and associated FileName(s):');
for i = 1:height(unique_pairs)
    % Check if DateTime is valid (not missing or NaT)
    if ~isnat(unique_pairs.DateTime(i))
        % Print valid DateTime
        fprintf('DateTime: %s, FileName: %s\n', datestr(unique_pairs.DateTime(i), 'yyyy-mm-dd'), unique_pairs.FileName{i});
    else
        % Handle invalid or missing DateTime
        fprintf('DateTime: [Invalid or Missing], FileName: %s\n', unique_pairs.FileName{i});
    end
end

%%
% Load the data from the CSV file
data_edited = readtable('ctd_salinity_and_temperature_2021_time.csv');

% Get unique FileNames
unique_files = unique(data_edited.FileName);

% Initialize tables for deep and shallow casts
deep_casts = [];
shallow_casts = [];

% Loop through each unique FileName
for i = 1:length(unique_files)
    % Extract data for the current FileName
    current_file = unique_files{i};
    file_data = data_edited(strcmp(data_edited.FileName, current_file), :);
    
    % Find the maximum pressure for this FileName
    max_pressure = max(file_data.Pressure);
    
    % Classify based on the maximum pressure
    if max_pressure >= 2000
        % Deep cast (Pressure >= 2000)
        deep_casts = [deep_casts; file_data];
    else
        % Shallow cast (Pressure < 2000)
        shallow_casts = [shallow_casts; file_data];
    end
end

% Save the deep casts to a CSV file
writetable(deep_casts, 'deep_casts_2021.csv');
disp('Deep casts saved to deep_casts.csv');

% Save the shallow casts to a CSV file
writetable(shallow_casts, 'shallow_casts_2021.csv');
disp('Shallow casts saved to shallow_casts.csv');
%%
new_datetime = datetime('2021-10-29', 'InputFormat', 'yyyy-MM-dd');

% Check if 'DateTime' column exists; if not, create it
if ~ismember('DateTime', data.Properties.VariableNames)
    data.DateTime = NaT(height(data), 1); % Initialize with NaT (Not-a-Time)
end

% Assign the new DateTime to rows where FileName is 'h3260201.ctd'
data.DateTime(strcmp(data.FileName, 'h333a0205.ctd')) = new_datetime;

% Optional: Save the updated table to a CSV file
writetable(data, 'ctd_salinity_and_temperature_2021_time_shallow.csv');
%% Load the existing shallow_casts_2021.csv data if it exists
if isfile('shallow_casts_2021.csv')
    existing_shallow_casts = readtable('shallow_casts_2021.csv');
else
    existing_shallow_casts = table();  % Create an empty table if the file doesn't exist
end

% Load the new data table
new_data = readtable('ctd_salinity_and_temperature_2021_time_shallow.csv');  % Replace with the actual CSV file

% Define the target FileName you want to filter
target_FileName = 'h333a0205.ctd';  % Replace this with the specific FileName you want

% Filter rows where FileName matches the target
filtered_data = new_data(strcmp(new_data.FileName, target_FileName), :);

% Append the new filtered data to the existing shallow_casts data
combined_shallow_casts = [existing_shallow_casts; filtered_data];

% Remove duplicate rows if any
combined_shallow_casts = unique(combined_shallow_casts, 'rows');

% Save the updated combined data back to shallow_casts_2021.csv
writetable(combined_shallow_casts, 'shallow_casts_2021.csv');

disp(['Data for ', target_FileName, ' has been appended to shallow_casts_2021.csv']);
%%
% Load the shallow casts data
shallow_data = readtable('shallow_casts_2021.csv');

% Load the deep data that needs DateTime assignment
deep_data = readtable('ctd_salinity_and_temperature_2021_combine_deep_data_filtered.csv');

% Check if 'DateTime' column exists in deep_data; if not, create it
if ~ismember('DateTime', deep_data.Properties.VariableNames)
    deep_data.DateTime = NaT(height(deep_data), 1);  % Initialize with NaT
end

% Get the unique FileNames from the shallow data
unique_files = unique(shallow_data.FileName);

% Loop through each unique FileName to assign corresponding DateTime
for i = 1:length(unique_files)
    current_file = unique_files{i};  % Current FileName
    
    % Get the unique DateTime associated with the current FileName
    date_for_file = unique(shallow_data.DateTime(strcmp(shallow_data.FileName, current_file)));
    
    % Assign the DateTime to matching FileName in the deep_data
    deep_data.DateTime(strcmp(deep_data.FileName, current_file)) = date_for_file;
end

% Save the updated deep_data table with the new DateTime assignments
writetable(deep_data, 'ctd_salinity_and_temperature_2021_combine_deep_data_filtered_time.csv');

% Display confirmation
disp('DateTime has been updated in the deep data file based on shallow_casts_2021.csv.');
