function [ matrix common lists] = SC_exprs_boolean( I_scores , kg_or_ens, treshold,option,depth)
% SC_exprs_boolean ************************************************************************
% GIOVANNI IACONO, CNAG, 16/08/2017
% writes to disk an excel file with the hierarchal markers of level "depth"
% INPUT
% I_scores: output of SC_calcola_markers
% kg_or_ens: annotatation file with gene names
% treshold: Treshold for significant changes. Default is 6.
% option: can be 
% -housekeeping to return the housekeeping genes (stably expressing genes in all cells)
% -active to return the markers expressed in cell types
% -silenced to return the genes silenced in cell types
% OUTPUT (optional)
% matrix: list of all significant markers
% common: list of all significant markers
% lists: list of all significant markers


if isequal(option,'housekeeping')
   disp('IGNORING THE depth parameter') 
else
    tresh_min=length(I_scores)-depth;
end


matrix=zeros( length(kg_or_ens),length(I_scores));

% Run trought the population of cells one by one
% every gene which is up-regulated in a certain population  compared to at least "total_population-depth" others becomes a marker 
for j=1:length(I_scores)
    

    block_values=[];
    for k=1:length(I_scores)
        if ~(k==j)    
            block_values=[block_values I_scores{j,k} ];
        end
    end
    block_values=-block_values;
    
    if isequal(option,'active')
        result=sum(block_values'>treshold);
    end
    if isequal(option,'silenced')
        result=sum(block_values'<-treshold);
        block_values=-block_values;
    end
     if isequal(option,'housekeeping')
        result=sum(abs(block_values')<treshold);
        tresh_min=length(block_values(1,:));
    end   
    [block_values_sorted] = sort(block_values');
    
    if depth<length(block_values_sorted(:,1))
    	dummy_scores_soft=min(block_values_sorted(depth:end,:))';
    else
        dummy_scores_soft=block_values_sorted(depth:end,:);
    end

    matrix(find(result>=tresh_min),j)=dummy_scores_soft(find(result>=tresh_min));
    
end


for j=1:length(I_scores)
    for h=1:length(I_scores)
        common(j,h)=nnz( sum(matrix(:,[j h])')==2);
    end
end

header={'GENE'};
for j=1:length(I_scores)
    lists{j}=find(matrix(:,j));
    header{1,j+1}=sprintf('C_%g',j);
end

% OLD CODE
% dummy=matrix;
% dummy(matrix==0)=NaN;
% order=min(dummy');
% order(isnan(order))=0;

order=sum(matrix');
[c ix]=sort(order,'descend');
matrix=matrix(ix,:);
kg_or_ens=kg_or_ens(ix,:);


% write to disk
xlswrite('./../data/SC_output_matrix_bool/matrix_bool.xlsx',header,sprintf('Lv%g',depth),'A1');
xlswrite('./../data/SC_output_matrix_bool/matrix_bool.xlsx',[kg_or_ens  mio_mat2cell(matrix)],sprintf('Lv%g',depth),'A2');

end


