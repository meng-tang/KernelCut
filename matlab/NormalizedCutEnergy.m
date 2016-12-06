function [ ncut_e ] = NormalizedCutEnergy( A, clustering )
%NormalizedCutEnergy computes normalized cut given affinity A and labeling
%Inputs:
%   A --- Affinity matrix (sparse or full) of N-by-N
%   clustering --- N-by-1 vector of cluster indicator, start from 1
%Outputs:
%   ncut_e --- Normalized Cut energy

K = max(clustering); % number of clusters
d = sum(A,2); % degree vector
% normalized cut is defined as \sum_k(assoc(S_k,\Omega - S_k) / assoc(S_k,\Omega))
% normalized association is defined as \sum_k(assoc(S_k,S_k) / assoc(S_k,\Omega))
nassoc_e = 0;
for k=1:K
    % binary indicators (N-by-1) for cluster k
    S_k = double(clustering == k);
    if 0 == sum(clustering==k); continue; end % skip empty cluster
    nassoc_e = nassoc_e + S_k' * A * S_k / (d' * S_k);
end
% ncut_e = K - nassoc_e
ncut_e = numel(unique(clustering))- nassoc_e;

end