function [  sum_ex ] = SC_bigdata_somma_ex( indices, indptr, data )
% SC_bigdata_somma_ex ************************************************************************
% GIOVANNI IACONO, CNAG, 16/08/2017
% Calculates library size for a matrix saved in sparse format



sum_ex=zeros(1,numel(indptr)-1);

if min(indptr)==0 & min(indices)==0
    disp('Fixing indexes');
    indptr=indptr+1;
    indices=indices+1;
end

for k=1:length(indptr)-1
    sum_ex(k)=sum(data(indptr(k) : indptr(k+1)-1));
end

% Old code for debuggibng
% for k=1:length(indptr)-1
%     factor=somma_ex(k)/mean(somma_ex);
%     data(indptr(k) : indptr(k+1)-1)=data(indptr(k) : indptr(k+1)-1)/factor;
% end
