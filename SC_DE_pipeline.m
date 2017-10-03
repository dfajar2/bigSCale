% SC_DE_pipeline ************************************************************************
% MAIN PIPELINE FOR DIFFERENTIAL EXPRESSION WITH bigSCale
% GIOVANNI IACONO, CNAG, 16/08/2017

% Requires to have:
% total_data: expression matrix (UMI or reads)
% idx_samples: cell array with the indexes of the conditions of the dataset (example, wild-type, ko)
% idx_pools: cell array with the indexes of the pools of the dataset (needed to perform batch correction)

% INITIATE THE VARIABLE EDGES

% For dataset with extremely high depth. This was used for the convoluted
% 10X genomics dataset which had an average library size of 238500 UMIs per cell 
% edges=[ 0 [0.000001:1:10] [11:2:26] [28:4:50] [55:5:100] [110:10:200] [220:20:580] [600:100:1000] [1200:200:2800] [3000:10000:100000] Inf ];

% For dataset with  good depth. This was used for the Zeisel brain dataset
% of 3005 cells
% edges=[ 0 [0.000001:1:10] [11:2:26] [28:4:50] [55:5:100] [110:10:200] [220:20:580] [600:100:1000] [1200:200:2800] [3000:10000:100000] Inf ];


% For sparse datasets. This was used for the 10x genomics brain dataset
% which had an average of 4800 UMIs per cells.
% of 3005 cells
% edges=[ 0 [0.000001:1:10] [11:2:26] [28:4:50] [55:5:100] [110:10:200] [220:20:580] [600:100:1000] [1200:200:2800] [3000:10000:100000] Inf ];

% OPTIONAL, remove batch effects. It is better no to do this step for DE
% only if the design is equilibrated (same proportion of batches in each group).
%total_data_batch=SC_remove_batches( total_data,idx_samples, idx_pools);

% CALCULATE THE MODEL
% Option 1, for datasets < 5000 cells
[~, N_pct] = SC_new_algorithm( total_data , edges, 0');

% Option 1, for big datasets > 5000 cells up to millions, different formnat
% for data
[~, N_pct ] = SC_new_algorithm_bigdata_v2( indices, indptr, data , edges)
close all

% PERFORM DE
D_scores=SC_DE_matlab_MEX_v2( total_data, N_pct , edges,sum(total_data),idx_samples, [],[],[]);
% ***************************************************************************************************************************************************





