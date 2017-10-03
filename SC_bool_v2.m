function all_markers=SC_bool_v2( I_scores , kg_or_ens, treshold,option)
% SC_bool_v2 ************************************************************************
% GIOVANNI IACONO, CNAG, 16/08/2017
% writes to disk an excel file with the hierarchal markers
% INPUT
% I_scores: output of SC_calcola_markers
% kg_or_ens: annotatation file with gene names
% SOGLIA,opzione: variable for SC_exprs_boolean
% OUTPUT
% list of all significant markers

% Maxium hierarchical level to look for. Default should be total clusters - 1
max_depth=length(I_scores)-1;

delete('./../data/SC_output_matrix_bool/matrix_bool.xlsx');
copyfile('./../data/SC_output_matrix_bool/BOOL_model.xlsx','./../data/SC_output_matrix_bool/matrix_bool.xlsx')

used=[];
% at every cycle calls SC_exprs_boolean to calculate the markers of level "depth"
for depth=1:4%max_depth
     depth
%     if depth<=1 
%         treshold_use=4;
%     else
%         treshold_use=treshold;
%     end
    [~ , ~, lists] = SC_exprs_boolean( I_scores , kg_or_ens, treshold,option,depth);
    used=[used ; cell2mat(lists')];
end

all_markers=unique(used);



end


