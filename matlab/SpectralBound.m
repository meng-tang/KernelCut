function [ unaries ] = SpectralBound( A, K, m, current_clustering, energy_type)
%KERNELBOUND derives linear upper bound w.r.t binary indicator variables
%for NC (Normalized Cut) or AA (Average Association). This function first
%finds spectral embedding (for once) and take distances to (weighted) mean
%on such embeddings as unaries.
%Inputs:
%   A --- Affinity matrix (sparse or full) of N-by-N
%   K --- Number of desired clusters, default is 2
%   m --- dimensionality of embedding, typically small
%   current_clustering --- N-by-1 vector of cluster indicator, value from 1
%   to K, default is a random generated vector
%   energy_type --- String 'NC' or 'AA', default is 'NC'
%Outputs:
%   unaries --- N-by-K matrix where unaries(n,k) is the unary (linear) cost
%   of assigning data n to cluster k


N = size(A,1); % number of data points
if nargin < 2
    K = 2;
end
if nargin < 3
    current_clustering = randi([1, K], N, 1);
end
if nargin < 4
    energy_type = 'NC';
end
if ~strcmp(energy_type, 'NC') && ~strcmp(energy_type, 'AA')
    error('Energy type has to be NC (normalized cut) or AA (average association).')
end

% for weighted k-means, only need to be computed ONCE.
persistent embedding weights
if isempty(embedding) || isempty(weights)
    [ embedding, weights ] = SpectralEmbedding( A, m, energy_type);
end

% update cluster centers for (weighted k-means)
centers = zeros(K,m);
weightedembedding = embedding .* repmat(weights,1,m);
for j=1:K
    centers(j,:) = sum( weightedembedding(current_clustering==j,:),1 ) ...
        / sum( weights(current_clustering==j));
end
% distances to centers
disttocenters = zeros(N,K);
for j=1:K
    disttocenters(:,j) = sum( (embedding - repmat(centers(j,:),N,1) ).^2 , 2);
end
% unaries are just distances to centers
unaries = disttocenters;

end