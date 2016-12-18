function [ embedding, weights ] = SpectralEmbedding( A, m, energy_type)
%SPECTRALEMBEDDING embeds graph to m-dimension space such that (weighted)
%k-means on such embedding approximates NC (Normalized Cut) or AA (Average
%Association)
%Inputs:
%   A --- Affinity matrix (sparse or full) of N-by-N
%   m --- dimensionality of embedding, typically small
%   energy_type --- String 'NC' or 'AA', default is 'NC'
%Outputs:
%   embedding --- N-by-m matrix n th row is the embedding for n th node
%   weights --- N-by-1 vector denoting weights of each embedding


N = size(A,1); % number of data points
if nargin < 2
    m = 2;
end
if nargin < 3
    energy_type = 'NC';
end
if ~strcmp(energy_type, 'NC') && ~strcmp(energy_type, 'AA')
    error('Energy type has to be NC (normalized cut) or AA (average association).')
end

opts.issym=1;
opts.isreal = 1;
opts.disp=0;
if strcmp(energy_type, 'AA') % for average association
    tic
    [EigVect, EVal] = eigs(A, m, 'lm',opts);
    embedding = EigVect * sqrt(EVal);
    weights = ones(N,1);
    toc
    return;
elseif strcmp(energy_type, 'NC') % for normalized cut
    % degree vector
    d = sum(A,2);
    % degree matrix
    D = sparse(1:N,1:N,d);
    if ~issparse(A); D = full(D); end;
    tic
    [EigVect, EVal] = eigs(D^(-0.5)*A*D^(-0.5), m, 'lm',opts);
    toc
    embedding = EigVect * sqrt(EVal)./repmat(sqrt(d),1,m);
    weights = d;
    return
end
end
