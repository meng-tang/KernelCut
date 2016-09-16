function [ embedding, weights ] = getembeddingandweights( EigVect, EigVal, embeddingmethod, num_vec, D )
%GETEMBEDDING Summary of this function goes here
%   Detailed explanation goes here
% EigVec: eigenvectors of generalized eigenproblem
% EigVal: eigenvalues for generalized eigenproblem, a column vector with
% values from small to large.
% embeddingmethod: eigenvector (just eigenvector) or weightedeigenvector
% (devided by square root of eigen values) or tang or bach
% D: For tang or bach, the degrees must be provided
N = size(EigVect, 1);
EigVect = EigVect ./ repmat(sqrt(sum(EigVect.^2)),N,1);
shift = 0;
if strcmp(embeddingmethod,'tang++')
    embeddingmethod = 'tang';
    shift = 1;
elseif strcmp(embeddingmethod,'tang--')
    embeddingmethod = 'tang';
    shift = -1;
end
if strcmp(embeddingmethod,'eigenvector')
    num_dim = num_vec;
    embedding = EigVect(:, 1:num_vec);
    weights = ones(N,1);
elseif strcmp(embeddingmethod,'weightedeigenvector')
    embedding = EigVect(:, 2:num_vec);
    embedding = embedding ./ repmat((sqrt(EigVal(2:num_vec)))',N,1);
    weights = ones(N,1);
    
    embedding = embedding(:,EigVal(2:num_vec)>eps);
    
elseif strcmp(embeddingmethod,'bach')
    num_dim = num_vec;
    embedding = EigVect(:, 1:num_vec);
    % eigenvector of normalized laplacian
    V = embedding .* repmat(sqrt(D),1,num_dim);
    V = V ./ repmat(sqrt(sum(V.^2)),N,1);
    embedding = V;
    embedding = embedding ./ repmat(sqrt(D),1,num_dim);
    weights = D;
elseif strcmp(embeddingmethod,'tang')
    num_dim = num_vec;
    embedding = EigVect(:, 1:num_vec);
    % eigenvector of normalized laplacian
    V = embedding .* repmat(sqrt(D),1,num_dim);
    V = V ./ repmat(sqrt(sum(V.^2)),N,1);
    embedding = V;
    eig_V = 1 - EigVal(1:num_vec);
    if shift == 1
        embedding = embedding .* repmat(sqrt(eig_V'+1),N,1);
    elseif shift == -1
        embedding = embedding .* repmat(sqrt(eig_V'-min(eig_V)),N,1);
    else
        embedding = embedding .* repmat(sqrt(eig_V'),N,1);
    end
    embedding = embedding ./ repmat(sqrt(D),1,num_dim);
    weights = D;
end
end

