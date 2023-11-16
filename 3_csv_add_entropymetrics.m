%% Add entropy metrics to data_long.csv
% Copywright David Reydellet - 27/09/2023


%% 1. time_series_complexity

%% Add a participantID column to the entropy metric tables
% LSD

% Load tables from the .mat files
directory = '/student/davidreydellet/lsd-basel/data/derivative/analysis/entropy/';
tableName = 'tbl_lsd_time_series_complexity.mat';
load([directory tableName]);

% Read the participant CSV file
csvfile = '/student/davidreydellet/lsd-basel/data/derivative/analysis/participant_list.csv';
participantList = readtable(csvfile);

% Merge the tables on the 'subject'-'number' columns and add the 'participant' column to the original table
tbl_lsd_time_series_complexity = join(tbl_lsd_time_series_complexity, participantList, 'LeftKey', 'subject', 'RightKey', 'number');

% Save the modified table back to .mat file
save([directory tableName], 'tbl_lsd_time_series_complexity');


% PLC
% Load tables from the .mat files
directory = '/student/davidreydellet/lsd-basel/data/derivative/analysis/entropy/';
tableName = 'tbl_plc_time_series_complexity.mat';
load([directory tableName]);

% Read the participant CSV file
csvfile = '/student/davidreydellet/lsd-basel/data/derivative/analysis/participant_list.csv';
participantList = readtable(csvfile);

% Merge the tables on the 'subject'-'number' columns and add the 'participant' column to the original table
tbl_plc_time_series_complexity = join(tbl_plc_time_series_complexity, participantList, 'LeftKey', 'subject', 'RightKey', 'number');

% Save the modified table back to .mat file
save([directory tableName], 'tbl_plc_time_series_complexity');


%% Add a column "time_series_complexity" to data_long.csv
% Load CSV file
csv_file_path = '/student/davidreydellet/lsd-basel/data/derivative/analysis/data_long.csv';
csv_data = readtable(csv_file_path);

% Load .mat files
lsd_file_path = '/student/davidreydellet/lsd-basel/data/derivative/analysis/entropy/tbl_lsd_time_series_complexity.mat';
plc_file_path = '/student/davidreydellet/lsd-basel/data/derivative/analysis/entropy/tbl_plc_time_series_complexity.mat';
load(lsd_file_path);
load(plc_file_path); 

% Initialize the new column
csv_data.time_series_complexity = nan(height(csv_data), 1);

% Initialize the new column
csv_data.time_series_complexity = nan(height(csv_data), 1);

% Loop through each row of the CSV file and assign the value of time_series_complexity
for i = 1:height(csv_data)
    participant = csv_data.participant{i}; % Adjusted assuming participant is a cell array of character vectors
    session = csv_data.session{i}; % Adjusted assuming session is a cell array of character vectors
    
    switch session
        case 'ses-lsd'
            % Find the corresponding row in the tbl_lsd_time_series_complexity
            idx = find(strcmp(tbl_lsd_time_series_complexity.participant, participant), 1);
            if ~isempty(idx)
                csv_data.time_series_complexity(i) = tbl_lsd_time_series_complexity.time_series_complexity(idx);
            end
            
        case 'ses-plc'
            % Find the corresponding row in the tbl_plc_time_series_complexity
            idx = find(strcmp(tbl_plc_time_series_complexity.participant, participant), 1);
            if ~isempty(idx)
                csv_data.time_series_complexity(i) = tbl_plc_time_series_complexity.time_series_complexity(idx);
            end
    end
end


% Write the updated table back to CSV
writetable(csv_data, csv_file_path);











%% 2. geodesic_entropy (path length distribution entropy)

%% Add a participantID column to the entropy metric tables
% LSD

% Load tables from the .mat files
directory = '/student/davidreydellet/lsd-basel/data/derivative/analysis/entropy/';
tableName = 'tbl_lsd_geodesic_entropy.mat';
load([directory tableName]);

% Read the participant CSV file
csvfile = '/student/davidreydellet/lsd-basel/data/derivative/analysis/participant_list.csv';
participantList = readtable(csvfile);

% Merge the tables on the 'subject'-'number' columns and add the 'participant' column to the original table
tbl_lsd_geodesic_entropy = join(tbl_lsd_geodesic_entropy, participantList, 'LeftKey', 'subject', 'RightKey', 'number');

% Save the modified table back to .mat file
save([directory tableName], 'tbl_lsd_geodesic_entropy');


% PLC
% Load tables from the .mat files
directory = '/student/davidreydellet/lsd-basel/data/derivative/analysis/entropy/';
tableName = 'tbl_plc_geodesic_entropy.mat';
load([directory tableName]);

% Read the participant CSV file
csvfile = '/student/davidreydellet/lsd-basel/data/derivative/analysis/participant_list.csv';
participantList = readtable(csvfile);

% Merge the tables on the 'subject'-'number' columns and add the 'participant' column to the original table
tbl_plc_geodesic_entropy = join(tbl_plc_geodesic_entropy, participantList, 'LeftKey', 'subject', 'RightKey', 'number');

% Save the modified table back to .mat file
save([directory tableName], 'tbl_plc_geodesic_entropy');


%% Add a column "geodesic_entropy" to data_long.csv
% Load CSV file
csv_file_path = '/student/davidreydellet/lsd-basel/data/derivative/analysis/data_long.csv';
csv_data = readtable(csv_file_path);

% Load .mat files
lsd_file_path = '/student/davidreydellet/lsd-basel/data/derivative/analysis/entropy/tbl_lsd_geodesic_entropy.mat';
plc_file_path = '/student/davidreydellet/lsd-basel/data/derivative/analysis/entropy/tbl_plc_geodesic_entropy.mat';
load(lsd_file_path);
load(plc_file_path); 

% Initialize the new column
csv_data.geodesic_entropy = nan(height(csv_data), 1);

% Initialize the new column
csv_data.geodesic_entropy = nan(height(csv_data), 1);

% Loop through each row of the CSV file and assign the value of geodesic_entropy
for i = 1:height(csv_data)
    participant = csv_data.participant{i}; % Adjusted assuming participant is a cell array of character vectors
    session = csv_data.session{i}; % Adjusted assuming session is a cell array of character vectors
    
    switch session
        case 'ses-lsd'
            % Find the corresponding row in the tbl_lsd_geodesic_entropy
            idx = find(strcmp(tbl_lsd_geodesic_entropy.participant, participant), 1);
            if ~isempty(idx)
                csv_data.geodesic_entropy(i) = tbl_lsd_geodesic_entropy.geodesic_entropy{idx}(27); %27: cf paper Navigating the chaos
            end
            
        case 'ses-plc'
            % Find the corresponding row in the tbl_plc_geodesic_entropy
            idx = find(strcmp(tbl_plc_geodesic_entropy.participant, participant), 1);
            if ~isempty(idx)
                csv_data.geodesic_entropy(i) = tbl_plc_geodesic_entropy.geodesic_entropy{idx}(27);
            end
    end
end


% Write the updated table back to CSV
writetable(csv_data, csv_file_path);





%% 3. meta_state_complexity

%% Add a participantID column to the entropy metric tables
% LSD

% Load tables from the .mat files
directory = '/student/davidreydellet/lsd-basel/data/derivative/analysis/entropy/';
tableName = 'tbl_lsd_meta_state_complexity.mat';
load([directory tableName]);

% Read the participant CSV file
csvfile = '/student/davidreydellet/lsd-basel/data/derivative/analysis/participant_list.csv';
participantList = readtable(csvfile);

% Merge the tables on the 'subject'-'number' columns and add the 'participant' column to the original table
tbl_lsd_meta_state_complexity = join(tbl_lsd_meta_state_complexity, participantList, 'LeftKey', 'subject', 'RightKey', 'number');

% Save the modified table back to .mat file
save([directory tableName], 'tbl_lsd_meta_state_complexity');


% PLC
% Load tables from the .mat files
directory = '/student/davidreydellet/lsd-basel/data/derivative/analysis/entropy/';
tableName = 'tbl_plc_meta_state_complexity.mat';
load([directory tableName]);

% Read the participant CSV file
csvfile = '/student/davidreydellet/lsd-basel/data/derivative/analysis/participant_list.csv';
participantList = readtable(csvfile);

% Merge the tables on the 'subject'-'number' columns and add the 'participant' column to the original table
tbl_plc_meta_state_complexity = join(tbl_plc_meta_state_complexity, participantList, 'LeftKey', 'subject', 'RightKey', 'number');

% Save the modified table back to .mat file
save([directory tableName], 'tbl_plc_meta_state_complexity');


%% Add a column "meta_state_complexity" to data_long.csv
% Load CSV file
csv_file_path = '/student/davidreydellet/lsd-basel/data/derivative/analysis/data_long.csv';
csv_data = readtable(csv_file_path);

% Load .mat files
lsd_file_path = '/student/davidreydellet/lsd-basel/data/derivative/analysis/entropy/tbl_lsd_meta_state_complexity.mat';
plc_file_path = '/student/davidreydellet/lsd-basel/data/derivative/analysis/entropy/tbl_plc_meta_state_complexity.mat';
load(lsd_file_path);
load(plc_file_path); 

% Initialize the new column
csv_data.meta_state_complexity = nan(height(csv_data), 1);

% Initialize the new column
csv_data.meta_state_complexity = nan(height(csv_data), 1);

% Loop through each row of the CSV file and assign the value of meta_state_complexity
for i = 1:height(csv_data)
    participant = csv_data.participant{i}; % Adjusted assuming participant is a cell array of character vectors
    session = csv_data.session{i}; % Adjusted assuming session is a cell array of character vectors
    
    switch session
        case 'ses-lsd'
            % Find the corresponding row in the tbl_lsd_meta_state_complexity
            idx = find(strcmp(tbl_lsd_meta_state_complexity.participant, participant), 1);
            if ~isempty(idx)
                csv_data.meta_state_complexity(i) = tbl_lsd_meta_state_complexity.meta_state_complexity(idx);
            end
            
        case 'ses-plc'
            % Find the corresponding row in the tbl_plc_meta_state_complexity
            idx = find(strcmp(tbl_plc_meta_state_complexity.participant, participant), 1);
            if ~isempty(idx)
                csv_data.meta_state_complexity(i) = tbl_plc_meta_state_complexity.meta_state_complexity(idx);
            end
    end
end


% Write the updated table back to CSV
writetable(csv_data, csv_file_path);







%% 4. sample_entropy



% 1. Read in the CSV files
data = readtable('MSSE.csv');
parcellationData = readtable('/student/davidreydellet/lsd-basel/scripts/CopBET-master/Atlases/shen_268_parcellation_networklabels_updated.csv');

% 1.bis Add NaNs for the missing scans
% Extract unique subject IDs
uniqueSubs = unique(data.sub);

% Create a template for a new row filled with NaNs
newRowTemplate = array2table(nan(1, width(data)));
newRowTemplate.Properties.VariableNames = data.Properties.VariableNames;

% Loop through each subject
for i = 1:length(uniqueSubs)
    sub = uniqueSubs{i}; % Get the subject ID from cell
    
    % Check if ses-plc session exists for this subject
    hasPlc = any(strcmp(data.sub, sub) & strcmp(data.ses, 'ses-plc'));
    
    % Check if ses-lsd session exists for this subject
    hasLsd = any(strcmp(data.sub, sub) & strcmp(data.ses, 'ses-lsd'));
    
    % If ses-plc session doesn't exist, add a new row with NaNs for this subject
    if ~hasPlc
        newRow = newRowTemplate;
        newRow.sub = {sub};
        newRow.ses = {'ses-plc'};
        data = [data; newRow];
    end
    
    % Similarly, if ses-lsd session doesn't exist, add a new row with NaNs for this subject
    if ~hasLsd
        newRow = newRowTemplate;
        newRow.sub = {sub};
        newRow.ses = {'ses-lsd'};
        data = [data; newRow];
    end
end

% 2. Gather regions into networks
numSubjects = height(data);
networkNames = unique(parcellationData.Network_name);

% Preallocate output table
outputData = table();

for i = 1:numSubjects
    for scale = 1:5
        for j = 1:length(networkNames)
            currentNetwork = strrep(networkNames{j}, '"', '');  % Remove double quotes
            regionsInNetwork = parcellationData.Node(strcmp(parcellationData.Network_name, currentNetwork));
            
            % 3. Compute the average MSSE for the regions in the current network for the current scale
            regionColumns = strcat('scale', num2str(scale), '_region', arrayfun(@num2str, regionsInNetwork, 'UniformOutput', false));
            avgMSSE = nanmean(data{i, regionColumns}, 2);  % Use nanmean to handle potential NaNs
            
            % Add results to output table
            tempTable = table(data.sub(i), data.ses(i), scale, string(currentNetwork), avgMSSE, ...
                'VariableNames', {'sub', 'ses', 'scale_MSSE', 'network_shen_MSSE', 'avg_MSSE'});
            outputData = [outputData; tempTable];
        end
    end
end
outputData.network_shen_MSSE = cellstr(outputData.network_shen_MSSE);

MSSE_network_long = outputData;

% 4. Write to new CSV file
cd('/student/davidreydellet/lsd-basel/data/derivative/analysis/entropy/')
writetable(outputData, 'MSSE_network_long.csv');







% Merge with data_long.csv

% 1. Read the two CSV files
MSSE_data = readtable('/student/davidreydellet/lsd-basel/data/derivative/analysis/entropy/MSSE_network_long.csv');
data_long = readtable('/student/davidreydellet/lsd-basel/data/derivative/analysis/data_long.csv');

% 2. Rename columns in MSSE_data for consistency in merging
MSSE_data.Properties.VariableNames{'sub'} = 'participant';
MSSE_data.Properties.VariableNames{'ses'} = 'session';

% 3. Merge the two datasets on participant and session
mergedData = outerjoin(MSSE_data, data_long, 'Keys', {'participant', 'session'}, 'MergeKeys', true);

% Save the merged data to a new CSV file
writetable(mergedData, '/student/davidreydellet/lsd-basel/data/derivative/analysis/data_long_withMSSE.csv');








%% 5. DCC

%% Add a participantID column to the entropy metric tables
% LSD

% Load tables from the .mat files
directory = '/student/davidreydellet/lsd-basel/data/derivative/analysis/entropy/';
tableName = 'tbl_lsd_DCC.mat';
load([directory tableName]);

% Read the participant CSV file
csvfile = '/student/davidreydellet/lsd-basel/data/derivative/analysis/participant_list.csv';
participantList = readtable(csvfile);

% Merge the tables on the 'subject'-'number' columns and add the 'participant' column to the original table
tbl_lsd_DCC = join(tbl_lsd_DCC, participantList, 'LeftKey', 'subject', 'RightKey', 'number');

% Save the modified table back to .mat file
save([directory tableName], 'tbl_lsd_DCC');


% PLC
% Load tables from the .mat files
directory = '/student/davidreydellet/lsd-basel/data/derivative/analysis/entropy/';
tableName = 'tbl_plc_DCC.mat';
load([directory tableName]);

% Read the participant CSV file
csvfile = '/student/davidreydellet/lsd-basel/data/derivative/analysis/participant_list.csv';
participantList = readtable(csvfile);

% Merge the tables on the 'subject'-'number' columns and add the 'participant' column to the original table
tbl_plc_DCC = join(tbl_plc_DCC, participantList, 'LeftKey', 'subject', 'RightKey', 'number');

% Save the modified table back to .mat file
save([directory tableName], 'tbl_plc_DCC');


%% Load tables and networks
cd('/student/davidreydellet/lsd-basel/data/derivative/analysis/entropy/')
load('tbl_lsd_DCC.mat');
load('tbl_plc_DCC.mat');

% Read CSV conversion file into a table
T = readtable('/student/davidreydellet/lsd-basel/scripts/CopBET-master/Atlases/shen_268_parcellation_networklabels_updated.csv');

% Extract the 'Network_name' column into a cell array
network_labels = table2cell(T(:, 'Network_name'));


%% Create the data_long_withDCC.csv

participants_lsd = unique(tbl_lsd_DCC.participant);
participants_plc = unique(tbl_plc_DCC.participant);
missing_lsd = setdiff(participants_plc, participants_lsd);

result = [];

% Process both tables
for session = {'ses-lsd', 'ses-plc'}
    if strcmp(session, 'ses-lsd')
        table = tbl_lsd_DCC;
    else
        table = tbl_plc_DCC;
    end

    unique_networks = unique(network_labels);
    
    for row = 1:height(table)
        entropy_matrix = table.entropy{row};
        network_matrix = zeros(length(unique_networks));

        for i = 1:length(unique_networks)
            for j = i:length(unique_networks)  % Start from i to include diagonal
                rows = strcmp(network_labels, unique_networks{i});
                cols = strcmp(network_labels, unique_networks{j});
                sub_matrix = entropy_matrix(rows, cols);
                avg_value = mean(sub_matrix(:));
                network_matrix(i, j) = avg_value;
                network_matrix(j, i) = avg_value;
            end
        end
        
        [row_ids, col_ids] = find(triu(network_matrix));  % Include diagonal
        
        for k = 1:length(row_ids)
            i = row_ids(k);
            j = col_ids(k);

            entry = [];
            entry.session = session{1};
            entry.participant = table.participant{row};
            entry.DCC = network_matrix(i,j);
            entry.network_shen = strcat(unique_networks{i}, '_', unique_networks{j});
            result = [result; entry];
        end
    end
end

% Add entries for participants missing ses-lsd with NaN DCC
for i = 1:length(missing_lsd)
    for j = 1:length(unique_networks)
        for k = j:length(unique_networks)  % Start from j instead of j+1
            entry = [];
            entry.session = 'ses-lsd';
            entry.participant = missing_lsd{i};
            entry.DCC = NaN;
            entry.network_shen = strcat(unique_networks{j}, '_', unique_networks{k});
            result = [result; entry];
        end
    end
end

% Convert result to table and save to CSV
result_table = struct2table(result);
dcc_long = result_table;
cd('/student/davidreydellet/lsd-basel/data/derivative/analysis/entropy/')
writetable(dcc_long, 'dcc_long.csv');


% Merge with data_long.csv
% Load the data from data_long.csv
data_long = readtable('/student/davidreydellet/lsd-basel/data/derivative/analysis/data_long.csv');

% Merge the tables based on 'participant' and 'session' columns
data_long_withDCC = innerjoin(dcc_long, data_long, 'Keys', {'participant', 'session'});

% If you want to save the merged table to a new CSV:
writetable(data_long_withDCC, '/student/davidreydellet/lsd-basel/data/derivative/analysis/data_long_withDCC.csv');