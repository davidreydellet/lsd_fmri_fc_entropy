%% FROM CONN TO CSV SPACE
% Script by David Reydellet 28/03/2023
% 1. Applies connReader (P.Fisher script) to extract a summary of Conn analyses
%    applied on a combination of two atlases
%    OUTPUT: "out_connReader_shen.mat" here: /student/davidreydellet/lsd-basel/data/derivative/analysis/network_conn_2atlases/
% 2. Transform this output and the 2-atlases matrix into a csv
%    OUTPUT: "conn.mat" and "conn.csv" here: /student/davidreydellet/lsd-basel/data/derivative/analysis/network_conn_2atlases/
% 3. Merge existing data (5D-ASC and demographics) to this csv
%    OUTPUT: "lsd-basel.csv" here: /student/davidreydellet/lsd-basel/data/derivative/analysis/


%% 1. connReader applied on shen atlas
% addpath('/users/patrick/github/connReader/matlab/')
addpath('/student/davidreydellet/lsd-basel/scripts/connReader/matlab')
infile = '/student/davidreydellet/lsd-basel/data/derivative/conn_project/20230926_conn_project_shen.mat';
atlas = 'shen_268';
out_connReader_shen = connExtract(infile, atlas);
cd /student/davidreydellet/lsd-basel/data/derivative/analysis/output_conn/
save out_connReader_shen out_connReader_shen




%% 2. Transform the 2-atlases matrix into a csv
cd /student/davidreydellet/lsd-basel/data/derivative/analysis/output_conn/
load('out_connReader_shen.mat')


% Create a 1D-array with connectivity values + a 1D-array with network pairs (strings)
conn_shen.connectivity_shen = [];
conn_shen.network = {};
conn_shen.participant = {};
conn_shen.session = {};

for ses = 1:length(out_connReader_shen.netconn)
    for sub = 1:length(out_connReader_shen.netconn{ses})
        
        % Initialize empty arrays
        unique_connectivity_shen = [];
        network_pairs = {};
        
        % Loop through the connectivity matrix
        for i = 1:size(out_connReader_shen.netconn{ses}{sub}.mean, 1)
            for j = i:size(out_connReader_shen.netconn{ses}{sub}.mean, 2) 
                
                % Get the connectivity value and network names for this pair
                connectivity_shen_value = out_connReader_shen.netconn{ses}{sub}.mean(i,j);
                network1_name = out_connReader_shen.netconn{ses}{sub}.x{i};
                network2_name = out_connReader_shen.netconn{ses}{sub}.x{j};

                % Add the connectivity_shen value to the unique_connectivity array if
                % it hasn't been added before
                if ~ismember(connectivity_shen_value, unique_connectivity_shen)
                    unique_connectivity_shen = [unique_connectivity_shen; connectivity_shen_value];
                end

                % Add the network pair to the network_pairs cell array
                network_pair_name = [network1_name, '_', network2_name];
                network_pairs{end+1,1} = network_pair_name;
                                
            end
        end
        
        conn_shen.connectivity_shen = [conn_shen.connectivity_shen; unique_connectivity_shen];
        conn_shen.network = [conn_shen.network; network_pairs];

        
        % Now we add participant and session field
        participant = {};
        session = {};
        
        str_array = out_connReader_shen.files{ses}{sub};
        if strcmp(str_array,'NA')            % some participants don't have LSD scans
            
            % Check the string of the previous subject, retrieve the
            % session and add +1 to the sub id
            str_array = out_connReader_shen.files{ses}{sub-1};
            sections = strsplit(str_array, '/');            
            participant_str = sections{8};
            session_str = sections{9};
            var_part = str2double(participant_str(end-1:end));
            new_var_part = var_part + 1;
            participant_str = [participant_str(1:end-2) sprintf('%02d', new_var_part)];

            % Then do the normal stuff
            for i = 1:length(unique_connectivity_shen)
                participant{i,1} = participant_str;
                session{i,1} = session_str;
            end
            conn_shen.participant = [conn_shen.participant; participant];
            conn_shen.session = [conn_shen.session; session];
        else
            sections = strsplit(str_array, '/');
            participant_str = sections{8};
            session_str = sections{9};
            for i = 1:length(unique_connectivity_shen)
                participant{i,1} = participant_str;
                session{i,1} = session_str;
            end
            conn_shen.participant = [conn_shen.participant; participant];
            conn_shen.session = [conn_shen.session; session];
        end
    end
end
        
save conn_shen conn_shen
writetable(struct2table(conn_shen), 'conn_shen.csv')



%% 3. Merge conn_shen with data_long

% 1. Load the CSV data into MATLAB
data_long = readtable('/student/davidreydellet/lsd-basel/data/derivative/analysis/data_long.csv');

% 2. Convert the conn_shen structure into a table
conn_shen_table = struct2table(conn_shen);

% 3. Rename the column network to network_shen
conn_shen_table.Properties.VariableNames{'network'} = 'network_shen';

% 4. Merge the tables based on the participant and session columns
merged_data = outerjoin(data_long, conn_shen_table, 'Keys', {'participant', 'session'}, 'MergeKeys', true);

data_long_withShen = merged_data;

% Save
writetable(merged_data, '/student/davidreydellet/lsd-basel/data/derivative/analysis/data_long_withShen.csv');