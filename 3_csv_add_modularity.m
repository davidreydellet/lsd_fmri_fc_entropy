%% Script to add modularity to data_long.csv

% Load data_long.csv file
cd('/student/davidreydellet/lsd-basel/data/derivative/analysis/')
data = readtable('data_long.csv');

% In the csv, create a participant index column and a session index column to map participant IDs and session names to indices
participantIDs = unique(data.participant);
[~, participant_idx] = ismember(data.participant, participantIDs);
data.participant_idx = participant_idx;
sessionIDs = unique(data.session);
[~, session_idx] = ismember(data.session, sessionIDs);
data.session_idx = session_idx;
clear sessioNIDs participant IDs participant_idx session_idx

% Find graph theory files
filename = '/student/davidreydellet/lsd-basel/data/derivative/analysis/graph_theory/positive_onlythreshold-0.3.mat';
cd('/student/davidreydellet/lsd-basel/data/derivative/analysis/graph_theory/')

load(filename);     % Load .mat file
varname = 'modularity_norm';
matvar = modularity_norm;

% Create new columns for each array in the data_long.csv file
for i = 1:size(matvar, 1) % loop over sessions
    for j = 1:size(matvar, 2) % loop over participants

        % Find rows in the data table corresponding to the current session and participant
        rows = (data.participant_idx == j) & (data.session_idx == i);

        % Extract the value from the current array for the current session and participant
        value = matvar(i,j);

        % Add a new column to the data table with the extracted value for the corresponding rows
        data.('normalized_modularity')(rows) = value;

    end
end

data = removevars(data, {'participant_idx', 'session_idx'});

% Write the updated data table to a new CSV file
cd('/student/davidreydellet/lsd-basel/data/derivative/analysis/')
writetable(data, 'data_long.csv');



% Done