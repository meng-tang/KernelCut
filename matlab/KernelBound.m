function [ unaries ] = KernelBound( A, K, current_clustering, energy_type)
%KERNELBOUND derives linear upper bound w.r.t binary indicator variables
%for NC (Normalized Cut) or AA (Average Association)
%Inputs:
%   A --- Affinity matrix (sparse or full) of N-by-N
%   K --- Number of desired clusters, default is 2
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
    current_labeling = randi([1, K], N, 1);
end
if nargin < 4
    energy_type = 'NC';
end
if ~strcmp(energy_type, 'NC') && ~strcmp(energy_type, 'AA')
    error('Energy type has to be NC (normalized cut) or AA (average association).')
end

unaries = zeros(N, K, 'double');
% degree vector
d = sum(A,2);
for i=1:K
    % current binary indicators (N-by-1) for cluster i
    S_t = double(current_clustering == i);
    % compute gradient as unaries (see kernel bound proposition in the paper)
    if strcmp(energy_type, 'NC') % for normalized cut
        unaries(:,i) = S_t'* A * S_t / (d' * S_t)^2 * d - 2 * A * S_t / (d' * S_t);
    elseif strcmp(energy_type, 'AA') % for average association
        unaries(:,i) = S_t'* A * S_t / (ones(1,N) * S_t)^2 * ones(N,1) - 2 * A * S_t / (ones(1,N) * S_t);
    end
end
end

