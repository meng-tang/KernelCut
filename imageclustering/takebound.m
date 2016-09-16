function [ dist_to_mu ] = takebound( W, wp,labeling , labelsets)
%TAKEBOUND Summary of this function goes here
%   Detailed explanation goes here
% W is sparse affinity matrix, row-wise
% wp is the degree of points
% labeling is current labeling, vector, row-wise
% take unary bound (dist_to_mu)
N = numel(labeling);
K = max(labeling(:))+1;
dist_to_mu = zeros(N,K,'double');
if nargin == 3
    labelsets = 1:K;
end
for k=labelsets
    s_k = int32(find(labeling==k)); % the set of points labeled k
    linear_sum_k = W(:,s_k);
    linear_sum_k = sum(linear_sum_k,2);
    linear_sum_k = full(linear_sum_k);
    quad_sum_k = sum(linear_sum_k(s_k));
    wp_sum_k = sum(wp(s_k));
    %dist_to_mu(:,k) = (1.0 ./ wp -2 * linear_sum_k / wp_sum_k + quad_sum_k / wp_sum_k / wp_sum_k * wp)*10000;
    dist_to_mu(:,k) = (-2 * linear_sum_k / wp_sum_k + quad_sum_k / wp_sum_k / wp_sum_k * wp);
    if isempty(s_k)
        dist_to_mu(:,k) = 1e+20;
    end
end

end

