function [ aassoc_e ] = AverageAssociationEnergy( A, clustering )
%AverageAssociationEnergy computes average assocaiton given affinity A and labeling
%Average association is defined as \sum_k(assoc(S_k,S_k) / size(S_k))
%Inputs:
%   A --- Affinity matrix (sparse or full) of N-by-N
%   clustering --- N-by-1 vector of cluster indicator, start from 1
%Outputs:
%   aassoc_e --- Average Association energy

K = max(clustering); % number of clusters

aassoc_e = 0;
for k=1:K
    % binary indicators (N-by-1) for cluster k
    S_k = double(clustering == k);
    if 0 == sum(clustering==k); continue; end % skip empty cluster
    aassoc_e = aassoc_e + S_k' * A * S_k / sum(S_k);
end

end