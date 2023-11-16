%% Computation of entropy metrics with CopBET toolbox
% David Reydellet - 25/09/2023
% 1. Compute LZ TIME COMPLEXITY + DCC + META-STATE COMPLEXITY   with ROI denoised timeseries as input data
% 2. Compute PATH LENGTH DISTRIBUTION                           with ROI-ROI connectivity matrix as input data
% 3. Compute SAMPLE ENTROPY                                     with Voxel denoised timeseries as input data

% Initialization
addpath(genpath('/student/davidreydellet/lsd-basel/scripts/CopBET-master/'))



%% 1. Compute LZ TIME COMPLEXITY + DCC with ROI denoised timeseries as input data

% Location of the ROI data
datafile = '/student/davidreydellet/lsd-basel/data/derivative/conn_project/20230926_conn_project_shen/results/preprocessing/';

% Find the names of all scans
filenames = dir([datafile,'/ROI_Subject*']);
scan_names_full = {filenames.name};     

% Organize these names in a cell scan_names_byCondition, one row per condition
scan_names_byCondition = cell(2,1);
condition_names = unique(cellfun(@(x) x(end-15:end), scan_names_full, 'UniformOutput', false));
condition_names = condition_names(2:3); % There is a condition000 in index1 to remove here
num_conditions = numel(condition_names);

% Loop through each condition
for c = 1:num_conditions
    % Find indices of strings that contain current condition
    condition_indices = endsWith(scan_names_full, condition_names{c});

    % Extract those strings into separate array
    condition_names_array = scan_names_full(condition_indices);

    % Add array to cell array
    scan_names_byCondition{c} = condition_names_array;
end



% Start by looping through conditions
for cond = 1:num_conditions
    scan_names = scan_names_byCondition{cond};
    subject_nb{cond} = cellfun(@str2double, unique(cellfun(@(x) x(end-19:end-17), scan_names, 'UniformOutput', false)));

    % Loop through scans/subjects for this condition
    for sub = 1:length(scan_names)
    
        scan = [datafile,'/',scan_names{sub}];
        load(scan);
        scan_names{sub};
        
        % remove ROIs which are not part of Schen
        data = data(:,4:end-3); 
        
        % In data, both conditions are piled up (condition feature in conn
        % seems broken), so we manually split it
        % Moreover, some subjects have got a missing scan for the lsd condition 
        % but they still have both condition files so I manually deal with
        % them
        missingSub = [30,60,64,71];
        numrois = size(data,2);

        % Convert cell array data into a dataMatrix (300 time x 268 ROI)
        if cond == 1 % lsd
            if ismember(sub, missingSub)
                dataMatrix = [];
            else
                for i = 1:numrois
                    dataMatrix(:, i) = data{i}(1:300);
                end
            end
        elseif cond == 2 % plc
            for i = 1:numrois
                if ismember(sub, missingSub)
                    dataMatrix(:, i) = data{i}(1:300); 
                else
                	dataMatrix(:, i) = data{i}(301:end);
                end
            end
        end
        
        dataCell{sub,cond} = dataMatrix;
        dataMatrix = [];
        
    end
end

% Create a 83*1 cell per condition
timeseries_lsd = dataCell(:, 1); % cond1 = lsd
timeseries_plc = dataCell(:, 2); % cond2 = plc

% Transform it into a table
tbl_plc = table(timeseries_plc);
tbl_plc.subject = subject_nb{2}.';

tbl_lsd = table(timeseries_lsd);
tbl_lsd.subject = subject_nb{1}.';
emptyRows = cellfun(@isempty, tbl_lsd.timeseries_lsd); % Find rows where the 'timeseries_lsd' variable is empty
tbl_lsd(emptyRows, :) = []; % Remove those rows from the table

% Good to keep: ROI names
ROInames = names(:,4:end-3);



% Compute LZ TIME COMPLEXITY + META-STATE COMPLEXITY + DCC

% Time series complexity
tbl_lsd = CopBET_time_series_complexity(tbl_lsd, 'LZ76temporal','keepdata',true,'parallel',true);
tbl_lsd.time_series_complexity = tbl_lsd.entropy;
tbl_plc = CopBET_time_series_complexity(tbl_plc, 'LZ76temporal', 'keepdata',true,'parallel',true);
tbl_plc.time_series_complexity = tbl_plc.entropy;

% Save tables
cd('/student/davidreydellet/lsd-basel/data/derivative/analysis/entropy/')
tbl_lsd_time_series_complexity = tbl_lsd;
tbl_plc_time_series_complexity = tbl_plc;
save('tbl_lsd_time_series_complexity.mat', 'tbl_lsd_time_series_complexity');
save('tbl_plc_time_series_complexity.mat', 'tbl_plc_time_series_complexity');



% Meta-state complexity
tbl_lsd = CopBET_metastate_series_complexity(tbl_lsd,'keepdata',true,'parallel',true);
tbl_lsd.meta_state_complexity = tbl_lsd.entropy;
tbl_plc = CopBET_metastate_series_complexity(tbl_plc, 'keepdata',true,'parallel',true);
tbl_plc.meta_state_complexity = tbl_plc.entropy;

% Save tables
cd('/student/davidreydellet/lsd-basel/data/derivative/analysis/entropy/')
tbl_lsd_meta_state_complexity = tbl_lsd;
tbl_plc_meta_state_complexity = tbl_plc;
save('tbl_lsd_meta_state_complexity.mat', 'tbl_lsd_meta_state_complexity');
save('tbl_plc_meta_state_complexity.mat', 'tbl_plc_meta_state_complexity');



% DCC entropy
cd('/student/davidreydellet/lsd-basel/data/derivative/analysis/entropy/')
tbl_lsd = CopBET_DCC_entropy(tbl_lsd,true,'keepdata',true,'parallel',true);
tbl_lsd.DCC_entropy = tbl_lsd.entropy;
tbl_lsd_DCC = tbl_lsd;
save('tbl_lsd_DCC.mat', 'tbl_lsd_DCC');

tbl_plc = CopBET_DCC_entropy(tbl_plc,true,'keepdata',true,'parallel',true);
tbl_plc.DCC_entropy = tbl_plc.entropy;
tbl_plc_DCC = tbl_plc;
save('tbl_plc_DCC.mat', 'tbl_plc_DCC');




%% 2. Compute PATH LENGTH DISTRIBUTION ENTROPY with ROI-ROI connectivity matrix as input data

% Location of the ROI data
datafile = '/student/davidreydellet/lsd-basel/data/derivative/conn_project/20230926_conn_project_shen/results/firstlevel/SBC_01';

% Find the names of all scans
filenames = dir([datafile,'/resultsROI_Subject*']);
scan_names_full = {filenames.name};     

% Organize these names in a cell scan_names_byCondition, one row per condition
scan_names_byCondition = cell(2,1);
condition_names = unique(cellfun(@(x) x(end-15:end-4), scan_names_full, 'UniformOutput', false));
num_conditions = numel(condition_names);

% Loop through each condition
for c = 1:num_conditions
    % Find indices of strings that contain current condition
    condition_indices = contains(scan_names_full, condition_names{c});

    % Extract those strings into separate array
    condition_names_array = scan_names_full(condition_indices);

    % Add array to cell array
    scan_names_byCondition{c} = condition_names_array;
end


% Start by looping through conditions
for cond = 1:num_conditions
    scan_names = scan_names_byCondition{cond};
    subject_nb{cond} = cellfun(@str2double, unique(cellfun(@(x) x(end-19:end-17), scan_names, 'UniformOutput', false)));

    % Loop through scans for this condition
    for sub = 1:length(scan_names)
    
        scan = [datafile,'/',scan_names{sub}];
        load(scan);
        scan_names{sub};
        
        numrois = size(Z,1);
        Z = Z(:,1:end-3); % remove ROIs which are not from Schen268
        
        if all(all(isnan(Z)))
            dataMatrix = [];
        else
            dataMatrix = Z;
        end 
        dataCell{sub,cond} = dataMatrix;
  
    end
end

% Create a 83*1 cell per condition
roiconnectivity_lsd = dataCell(:, 1);
roiconnectivity_plc = dataCell(:, 2);

% Transform it into a table
tbl_plc = table(roiconnectivity_plc);
tbl_plc.subject = subject_nb{2}.';

tbl_lsd = table(roiconnectivity_lsd);
tbl_lsd.subject = subject_nb{1}.';
emptyRows = cellfun(@isempty, tbl_lsd.roiconnectivity_lsd); % Find rows where the 'roiconnectivity_lsd' variable is empty
tbl_lsd(emptyRows, :) = []; % Remove those rows from the table


% Replace the NaN in the diagonals of the connectivity matrices with 1
numRows = height(tbl_lsd);
for i = 1:numRows
    % Extract the matrix from the cell in the timeseries_lsd column of the table
    matrix = tbl_lsd.roiconnectivity_lsd{i};
    % Replace NaN values in the diagonal with 1
    for j = 1:268
        if isnan(matrix(j,j)) 
            matrix(j,j) = 1;
        end
    end
    % Store the modified matrix back in the table
    tbl_lsd.roiconnectivity_lsd{i} = matrix;
end
numRows = height(tbl_plc);
for i = 1:numRows
    % Extract the matrix from the cell in the timeseries_plc column of the table
    matrix = tbl_plc.roiconnectivity_plc{i};
    % Replace NaN values in the diagonal with 1
    for j = 1:268
        if isnan(matrix(j,j)) 
            matrix(j,j) = 1;
        end
    end
    % Store the modified matrix back in the table
    tbl_plc.roiconnectivity_plc{i} = matrix;
end



% Compute Path length distribution (from ROI-ROI connectivity matrices)
% Geodosic entropy - path length entropy
tbl_lsd = CopBET_geodesic_entropy(tbl_lsd, 'keepdata',true,'parallel',true);
tbl_lsd.geodesic_entropy = tbl_lsd.entropy; 
tbl_lsd.entropy = []; 
tbl_plc = CopBET_geodesic_entropy(tbl_plc, 'keepdata',true,'parallel',true);
tbl_plc.geodesic_entropy = tbl_plc.entropy; 
tbl_plc.entropy = []; 

% Save tables
cd('/student/davidreydellet/lsd-basel/data/derivative/analysis/entropy/')
tbl_lsd_geodesic_entropy = tbl_lsd;
tbl_plc_geodesic_entropy = tbl_plc;
save('tbl_lsd_geodesic_entropy.mat', 'tbl_lsd_geodesic_entropy');
save('tbl_plc_geodesic_entropy.mat', 'tbl_plc_geodesic_entropy');







%% 3. Compute SAMPLE ENTROPY with Voxel denoised timeseries as input data

% Define atlas
atlas = niftiread('/student/davidreydellet/lsd-basel/scripts/CopBET-master/Atlases/Shen268_2mm.nii');

% Define the base directory where subject folders are located
baseDir = "/student/davidreydellet/lsd-basel/data/derivative/fmriprep-fastsurfer/";

% Use dir to list all sub-XXX folders
subFolders = dir(fullfile(baseDir, 'sub-*'));

% Filter out non-directories and potential *.html or other files
subjects = {subFolders([subFolders.isdir]).name};

% Pre-allocate cell arrays to hold the paths for each session type
paths_lsd = cell(length(subjects), 1);
paths_plc = cell(length(subjects), 1);

% Subjects without an LSD scan
missingSub = [30,60,64,71];

% Generate the paths
for s = 1:length(subjects)
    % Common path for both sessions for the current subject
    commonPart = fullfile(baseDir, subjects{s});
    
    if ismember(s,missingSub)
        % For ses-plc
        paths_plc{s} = fullfile(commonPart, 'ses-plc', 'func', ...
        sprintf('d%s_ses-plc_task-rest_acq-ep2d_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.nii', subjects{s}));
    
        % For ses-lsd
        paths_lsd{s} = [];
        
    else
        % For ses-plc
        paths_plc{s} = fullfile(commonPart, 'ses-plc', 'func', ...
        sprintf('d%s_ses-plc_task-rest_acq-ep2d_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.nii', subjects{s}));
        
        % For ses-lsd
        paths_lsd{s} = fullfile(commonPart, 'ses-lsd', 'func', ...
        sprintf('d%s_ses-lsd_task-rest_acq-ep2d_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.nii', subjects{s}));
    end
end

% Create the tables
tbl_lsd_sample = table(paths_lsd, 'VariableNames', {'FilePath_DenoisedVoxelTimeSeries'});
tbl_plc_sample = table(paths_plc, 'VariableNames', {'FilePath_DenoisedVoxelTimeSeries'});

emptyRows = cellfun(@isempty, tbl_lsd_sample.FilePath_DenoisedVoxelTimeSeries); % Find rows where the 'roiconnectivity_lsd' variable is empty
tbl_lsd_sample(emptyRows, :) = []; % Remove those rows from the table

% Retrieve participant and subject IDs
entropyDir = '/student/davidreydellet/lsd-basel/data/derivative/analysis/entropy/';

% Load the tables
loadedData = load(fullfile(entropyDir, 'tbl_lsd_time_series_complexity.mat'));
tbl_lsd_time_series_complexity = loadedData.tbl_lsd_time_series_complexity;

loadedData = load(fullfile(entropyDir, 'tbl_plc_time_series_complexity.mat'));
tbl_plc_time_series_complexity = loadedData.tbl_plc_time_series_complexity;

% Copy the required columns
tbl_lsd_sample.subject = tbl_lsd_time_series_complexity.subject;
tbl_lsd_sample.participant = tbl_lsd_time_series_complexity.participant;

tbl_plc_sample.subject = tbl_plc_time_series_complexity.subject;
tbl_plc_sample.participant = tbl_plc_time_series_complexity.participant;


% Compute sample entropy
tbl_lsd_sample = CopBET_sample_entropy(tbl_lsd_sample, atlas, true, 'keepdata',true,'parallel',true);
tbl_plc_sample = CopBET_sample_entropy(tbl_plc_sample, atlas, true, 'keepdata',true,'parallel',true);

% Save tables
cd('/student/davidreydellet/lsd-basel/data/derivative/analysis/entropy/')
save('tbl_lsd_sample.mat', 'tbl_lsd_sample');
save('tbl_plc_sample.mat', 'tbl_plc_sample');




