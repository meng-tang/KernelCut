function [ nass_e ] = ncutEnergy( W, labeling )
%NCUTE Summary of this function goes here
%   NCUT energy
[wx, wy] = size(W);
x = 1 : wx;
S = full(sum(W, 1));
D = sparse(x, x, S, wx, wy);
clear S x;
K = max(labeling);
nume = zeros(1,K,'double');
demo = zeros(1,K,'double');
[rows, cols, values] = find(W);
l_rows = labeling(rows);
l_cols = labeling(cols);
sameidx = (l_rows == l_cols);
rows = rows(sameidx);
%cols = cols(diffidx);
values = values(sameidx);
l_rows = labeling(rows);
for k=1:K
    nume(k) = sum(values(l_rows == k));
end
degrees = full(diag(D));
for k=1:K
    demo(k) = sum(degrees(labeling == k));
end
nass_e = - sum( nume ./ (demo + 1e-10));
end