function [ energy ] = scutenergy( embedding, Labeling, opt, G, goodgroups)
%SCUTENERGY Summary of this function goes here
%   Detailed explanation goes here
num_cluster = max(Labeling(:));
N = size(embedding,1);
num_dim = size(embedding,2);
numGroups = size(G,2);
if nargin==4
    goodgroups = 1:numGroups
end
% update centers
centers = zeros(num_cluster,num_dim);
for j=1:num_cluster
    centers(j,:) = sum( embedding(Labeling==j,:),1 ) / sum( Labeling==j);
end
% distances to centers
disttocenters = zeros(N,num_cluster);
for j=1:num_cluster
    disttocenters(:,j) = sum( (embedding - repmat(centers(j,:),N,1) ).^2 , 2);
end
kmeans_e = 0;
for j=1:num_cluster
    kmeans_e  = kmeans_e + sum(disttocenters(Labeling==j,j));
end

pnpotts_e = 0;
if opt.pnpottsweight > 1e-10
    for groupid = 1:numGroups
        if numel(find(goodgroups==groupid))~=1; continue; end
        groupindicator = (G(:,groupid)==1);
        grouphist = zeros(1,num_cluster,'int32');
        for j=1:num_cluster
            grouphist(j) = sum(Labeling(groupindicator)==j);
        end
        C = double(sum(groupindicator==1));
        num_dominant = double(max(grouphist));
        num_disagreement = C - num_dominant;
        Q = floor(C * opt.pnpottsratio);
        theta  = opt.pnpottsmaxpen / Q;
        pnpotts_e = pnpotts_e + min(opt.pnpottsmaxpen, theta * num_disagreement);
    end
end

energy = kmeans_e + pnpotts_e*opt.pnpottsweight;
end

