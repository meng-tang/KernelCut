function [ energy ] = kcutenergy( W, Labeling, opt, G, goodgroups)
%SCUTENERGY Summary of this function goes here
%   Detailed explanation goes here
num_cluster = max(Labeling(:));
N = size(W,1);
numGroups = size(G,2);
if nargin==4
    goodgroups = 1:numGroups
end
ncutE = ncutEnergy(W,Labeling);
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

energy = ncutE + pnpotts_e*opt.pnpottsweight;
end

