% SC_pipeline ************************************************************************
% MAIN PIPELINE FOR CELL CLUSTERING AND MARKERS EXTRACTION WITH bigSCale
% GIOVANNI IACONO, CNAG, 16/08/2017

% Requires to have:
% total_data: expression matrix (UMI or reads)
% idx_samples: cell array with the indexes of the conditions of the dataset (example, wild-type, ko)
% idx_pools: cell array with the indexes of the pools of the dataset (needed to perform batch correction)

% INITIATE THE VARIABLE EDGES

% For dataset with extremely high depth. This was used for the convoluted
% 10X genomics dataset which had an average library size of 238500 UMIs per cell 
% edges=[ 0 [0.000001:1:10] [11:2:26] [28:4:50] [55:5:100] [110:10:200] [220:20:580] [600:100:1000] [1200:200:2800] [3000:10000:100000] Inf ];

% For dataset with  good depth. This was used for the Zeisel brain dataset of 3005 cells
% edges=[ 0 [0.000001:1:10] [11:2:26] [28:4:50] [55:5:100] [110:10:200] [220:20:580] [600:100:1000] [1200:200:2800] [3000:10000:100000] Inf ];


% For sparse datasets. This was used for the 10x genomics brain dataset
% which had an average of 4800 UMIs per cells
% edges=[ 0 [0.000001:1:10] [11:2:26] [28:4:50] [55:5:100] [110:10:200] [220:20:580] [600:100:1000] [1200:200:2800] [3000:10000:100000] Inf ];


% OPTIONAL, remove batch effects
total_data_batch=SC_remove_batches( total_data,idx_samples, idx_pools);

% CALCULATE THE MODEL
% Option 1, for datasets < 5000 cells
[~, N_batch_pct] = SC_new_algorithm( total_data_batch , edges, 0);

% Option 1, for big datasets > 5000 cells up to millions, differnt formnat
% for data
[~, N_batch_pct ] = SC_new_algorithm_bigdata_v2( indices, indptr, data , edges)
close all
   
% CALCULATE OVERDISPERSED GENES
% option 1: calculate simply overdispersed genes
[driving_genes, ~ ,score]= population_calling_v2( total_data_batch,[], 0 );
% option 2: returns only the overdispersed genes displaying a sufficent number of correlated genes
[driving_genes, genes_related, score]= population_calling_v2( total_data_batch,[],1 );

% CALCULATE THE DISTANCE MATRIX OVER THE OVERDISPERSED GENES
close all
[ D_ls_mex ]=SC_1vs1_MEX( total_data_batch(driving_genes,:), N_batch_pct , edges,sum(total_data_batch));
close all
 
% Depth at which Dendrogram is cut, the lower the numner the higher the
% cluster numbers
Depth=0.20;
 
% CALCULATE THE DENDROGRAM AND THE CLUSTERS
Z = linkage(D_ls,'ward'); 
[~,~,outperm]=dendrogram(Z,Inf,'ColorThreshold', Depth*max(Z(:,3)));
T = cluster(Z,'cutoff',Depth*max(Z(:,3)),'criterion','distance');
indices={};
names={};
count=1;
for k=1:max(T)
    if (nnz(T==k)>10)
        indices{count}=find(T==k);
        names{count}=sprintf('Cluster_%g',count);
        count=count+1;
    end
end
ix =  SC_info_dendro( outperm,indices );
indices=indices(ix);
close all

% Optional
% CALCULATE AND WRITES TO DISK THE MARKERS OF LV.1 (markers_hard) AND LV.2 (markers_soft)
[markers_hard markers_soft scores_hard scores_soft I_scores common]= SC_calculate_markers( total_data_batch, indices, N_batch_pct , sum(total_data_batch), 4 );


% CALCULATES AND WRITES TO DISK THE COMPETE STRUCTURE OF HIERARCHICAL  MARKERS
all_markers=SC_bool_v2( I_scores , kg_or_ens, 6,'active');

% ***************************************************************************************************************************************************





