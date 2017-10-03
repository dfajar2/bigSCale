function [ total_data , indices, indptr, data] = SC_bigdata_reduce( indices, indptr, data , genes )
% SC_bigdata_reduce ************************************************************************
% GIOVANNI IACONO, CNAG, 16/08/2017
% Transforms Expression matrix from a sparse form to a full form.
% INPUT
% indices,indptr, data: sparse representation of the expression matrix
% genes: use only "genes" genes in the full matrix  
% OUTPUT
% total_data: Expression matrix full.

indices = int16(indices);
data = int16(data);
genes = int16(genes);

% Forcing indexes to start from 1 and not 0
if min(indptr)==0 & min(indices)==0
    disp('Fixing indexes');
    indptr=indptr+1;
    indices=indices+1;
end

num_samples=length(indptr)-1



result=ismember(indices,genes);


detected=zeros(1,num_samples);

for k=1:length(indptr)-1
    detected(k)=sum(result(indptr(k) : indptr(k+1)-1));
end

if min(detected==0)
    error('Strange mistake');
end

indices=indices(result);
[~, ~, indices]=unique(indices);
data=data(result);
indptr=cumsum([1 detected]);

total_data=zeros(numel(genes),num_samples ,'single');

% creating the full matrix
for k=1:num_samples
    total_data(indices( indptr(k) : indptr(k+1)-1 ),k)=data( indptr(k) : indptr(k+1)-1 );
end