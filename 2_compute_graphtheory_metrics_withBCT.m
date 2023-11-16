addpath('/student/davidreydellet/lsd-basel/scripts/brain_connectivity_toolbox/2019_03_03_BCT/');
addpath('/student/davidreydellet/lsd-basel/scripts/SmallWorldNess-master/');
% thresholds = linspace(0.2, 0.5, 4);
thresholds = 0.3;
        
for tt = 1:length(thresholds)

    threshold = thresholds(tt)

    datafile = '/student/davidreydellet/lsd-basel/data/derivative/conn_project/conn_project01/results/firstlevel/SBC_01/';


    %Find the names of all scans
    filenames = dir([datafile,'/resultsROI_Subject*']);
    scan_names_full = {filenames.name};     

    % Organize these names in a cell scan_names_byCondition, one row per condition
    scan_names_byCondition = cell(2,1);
    condition_names = unique(cellfun(@(x) x(end-15:end), scan_names_full, 'UniformOutput', false));
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



    % Load the session connectivity matrix and calculate all metrics

    % Start by looping through conditions
    for iii = 1:num_conditions
        scan_names = scan_names_byCondition{iii};

        % Loop through scans for this condition
        for ii = 1:length(scan_names)

            scan = [datafile,'/',scan_names{ii}];
            load(scan);
            scan_names{ii}

            if ~all(isnan(Z(:))) % we don't analyse empty scans (matrix filled with NaN only)

                z_mat = Z(1:100,1:100); % we take only schaefer100parcels

                numrois = length(z_mat);

                %Make the diagonal zero
                for i = 1:length(z_mat)
                    z_mat(i,i) = 0;
                end

                
                %Calculate Global Functional Connectivity (GFC)
                lower_triangle = tril(z_mat, -1); % Extract the lower triangle elements
                gfc{iii,ii} = mean(lower_triangle(:)); % Calculate the mean 
                

                %Calculate clustering coefficient per ROI
                clustcoeff = clustering_coef_wu(z_mat)';
                %Calculate the mean clustering coefficient for all ROIs
                meancc{iii,ii} = mean(clustcoeff);

                %Threshold the matrix
                mat_thresh = z_mat;
                for i = 1:length(z_mat)
                    for j = 1:length(z_mat)
                        if z_mat(i,j) < threshold
                            mat_thresh(i,j) = 0;
                        end
                    end
                end

                %Calculate mean degree
                degree = degrees_und(mat_thresh)';
                meandegree{iii,ii} = mean(degree);

                %Create distance matrix
                mat_length = z_mat;
                for i = 1:length(z_mat)
                    for j = 1:length(z_mat)
                        if z_mat(i,j) ~= 0
                            mat_length(i,j) = 1/z_mat(i,j);
                        end
                    end
                end

                %calculate path lengths
                cpl = distance_wei(mat_length);

                %calculate efficiencies
                eff_mat = cpl;
                for i = 1:length(z_mat)
                    for j = 1:length(z_mat)
                        if cpl(i,j) ~= 0
                            eff_mat(i,j) = 1/cpl(i,j);
                        end
                    end
                end

                %find the mean efficiencies (ignoring the diagonal)
                glob_eff{iii,ii} = mean(sum(eff_mat)/(numrois-1));

                %find the mean path lengths (ignoring the diagonal)
                meancpl{iii,ii} = mean(sum(cpl)/(numrois-1));

                %calculate local efficiency
                %the 2 denotes local efficiency, revised version
                mean_local_eff{iii,ii} = mean(efficiency_wei(z_mat,2));


                %Normalised modularity
                n_rep = 100;
                n_null = 100;
                Qb = 0; %initialize Qb for this subject

                for rep = 1:n_rep %repeat the community detection algorithm n_rep times

                    A = z_mat;%create a positive-values-only matrix based on
                    %the input matrix M for this subject

                    [~, Qt] = community_louvain(A, 1, [], 'negative_sym'); %apply the community detection
                    %algorithm

                    if Qt > Qb %if the current modularity value is higher than the previous
                        %highest value

                        Qb = Qt;%then set Qb to the current modularity value, to maximize
                        %modularity

                    end

                end

                q_max = Qb; %store the highest modularity value in the q variable for this
                %subject


                q_null = zeros(1, n_null); %initialize the null modularity to zero


                if sum(isnan(A(:))) > 0 %check if there are any NaN values in the
                    %adjacency matrix

                    error('Adjacency matrix contains NaNs') %set the null modularity value
                    %to zero if there are any NaN values
                else

                    for null = 1 : n_null %repeat the null model creation n_null times

                        B = randmio_und(A, 1); %create a null model on
                        %the binary adjacency matrix A

                        Qb = 0; %initialize the highest modularity value to zero

                        for rep = 1: n_rep %repeat the community detection
                            %algorithm n_rep times on the null model

                            [~, Qt] = community_louvain(B, 1, [], 'negative_sym');%apply the
                            %community detection algorithm and get the
                            %modularity value Qt

                            if Qt > Qb %if the current modularity value is
                                %higher than the previous highest value

                                Qb = Qt;%set Qb to be the current modularity value
                            end
                        end

                        q_null(1, null) = Qb;%store the highest modularity
                        %value in the null modularity value array
                    end
                end

                q_null_mean = mean(q_null); %take the mean of the null modularity
                %values and store it in the q_null_mean matrix

                q_norm = q_max / q_null_mean; %calclualtes the normalized modularity values
                %by dividing the mean modularity values by the mean null modularity values

                modularity_norm{iii,ii} = q_norm;
                modularity{iii,ii} = q_max;


                %Smallworldness

                % analysis parameters
                Num_ER_repeats = 100;  % to estimate C and L numerically for E-R random graph
                Num_S_repeats = 1000; % to get P-value for S; min P = 0.001 for 1000 samples
                I = 0.95;

                FLAG_Cws = 1;
                FLAG_Ctransitive = 2;

                A = mat_thresh;

                % get its basic properties
                n = size(A,1);  % number of nodes
                k = sum(A);  % degree distribution of undirected network
                m = sum(k)/2;
                K = mean(k); % mean degree of network

                [expectedC,expectedL] = ER_Expected_L_C(K,n);  % L_rand and C_rand

                [S_ws,C_ws,L] = small_world_ness(A,expectedL,expectedC,FLAG_Cws);  % Using WS clustering coefficient

                smallworldness{iii,ii} = S_ws;


            else % if no scan, ie matrix is filled with NaN
                gfc{iii,ii} = NaN;
                glob_eff{iii,ii} = NaN;
                mean_local_eff{iii,ii} = NaN;
                meancc{iii,ii} = NaN;
                meancpl{iii,ii} = NaN;
                meandegree{iii,ii} = NaN;
                modularity_norm{iii,ii} = NaN;
                modularity{iii,ii} = NaN;
                smallworldness{iii,ii} = NaN;
            end

        end
    end

    save(['/student/davidreydellet/lsd-basel/data/derivative/analysis/graph_theory/posneg_threshold-',num2str(threshold),'.mat'], ...
        'threshold', 'gfc', 'meancc','meandegree', 'glob_eff', 'meancpl', 'mean_local_eff',...
        'modularity_norm', 'modularity', 'smallworldness')

    clearvars gfc meancc meandegree glob_eff meancpl mean_local_eff modularity_norm modularity smallworldness

end





