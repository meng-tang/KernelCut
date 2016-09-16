%add libs
addpath('../libs/Bk_matlab/');
addpath('../libs/Bk_matlab/bin');
load('kernel.mat');
load('kkm_sc.mat');
load('labels.mat');
N = size(W,1);
opt.numClusters = 8;

%Labeling = randi([1,opt.numClusters],N,1); % random init
%load('sc.mat');
itrNum = 0;
nmis = nmi(Labeling,labels);
wp = sum(W,2);
rng(100);

%Labeling = labels;
% for shallow kernel
% opt.pnpottsweight = 1.0e-3;
% opt.pnpottsratio = 0.25;
% opt.pnpottsmaxpen = 10;
% for deep kernel
opt.pnpottsweight = 5e-3;
opt.pnpottsratio = 0.3;
opt.pnpottsmaxpen = 10;
load G19.mat
goodgroups = [1 2 3 7 8 9 10 12 13 14 18 19];
numGroups = size(G,2);

% we sample the grouping information
samplerate = 1.0
sampledgroup = randperm(int32(N),int32(N*(1.0-samplerate)));
G(sampledgroup,:) = zeros(numel(sampledgroup),numGroups);

rng(100);

swap = ones(opt.numClusters,opt.numClusters,'int32'); % set of pair of labels for swapping
swap(logical(eye(opt.numClusters))) = 0;
tic
current_e = 1e+10;
kcutES = [];
itrmNum = 0;
while sum(swap(:))~=0
    itrNum = itrNum + 1;
    for alpha = randperm(opt.numClusters)
        for beta = randperm(opt.numClusters)
            if swap(alpha,beta) == 0 || alpha == beta
                continue;
            end
            vLabeling = Labeling(:);
            s_alpha = find(vLabeling == alpha);
            if isempty(s_alpha)
                swap(alpha, :) = 0;swap(:,alpha) = 0;
                continue;
            end
            s_beta = find(vLabeling == beta);
            if isempty(s_beta)
                swap(beta, :) = 0;swap(:,beta) = 0;
                continue;
            end
            roi = [s_alpha ; s_beta]';
            [ dist_to_mu ] = takebound( W, wp,vLabeling,[alpha beta] );
            dist_to_alpha = dist_to_mu(:,alpha);
            dist_to_beta = dist_to_mu(:,beta);            
            % set of points labeled alpha or beta
            alphabetasize = numel(roi);
            h = BK_Create(alphabetasize+numGroups*2,8*(alphabetasize+numGroups*2));
            dist_to_alpha = dist_to_alpha(roi);
            dist_to_beta = dist_to_beta(roi);
            unaryterm = [dist_to_alpha' ; dist_to_beta'];
            unaryterm = [unaryterm zeros(2,numGroups*2)];
            W_pairwise = [];
            if opt.pnpottsweight>0
                Edges = [];%NumEdges-by-6 matrix [i,j,e00,e01,e10,e11]
                for groupid = 1:numGroups
                    groupindicator = (G(:,groupid)==1);
                    %groupindicator = (labels==groupid);
                    roigroupindicator = groupindicator(roi);
                    C = sum(groupindicator); % group size
                    C_ab = sum(roigroupindicator==1);
                    Q = floor(C * opt.pnpottsratio);
                    theta  = opt.pnpottsmaxpen / Q;
                    if numel(find(goodgroups==groupid))~=1; continue; end
                    if C_ab <= C -Q; continue; end
                    Edgesonegroup = zeros(C_ab*2+1,6);
                    Edgesonegroup(1:C_ab,1) = find(roigroupindicator==1);
                    Edgesonegroup(1:C_ab,2) = alphabetasize + (groupid-1)*2 +1;
                    Edgesonegroup(1:C_ab,3:6) = repmat([0 0 theta 0],C_ab,1);
                    Edgesonegroup((1+C_ab):(2*C_ab),1) = alphabetasize + (groupid-1)*2 +2;
                    Edgesonegroup((1+C_ab):(2*C_ab),2) = find(roigroupindicator==1);
                    Edgesonegroup((1+C_ab):(2*C_ab),3:6) = repmat([0 0 theta 0],C_ab,1);
                    Edgesonegroup(end,1:6) = [alphabetasize + 1 + (groupid-1)*2 alphabetasize + 2 + (groupid-1)*2 ...
                        0 0 opt.pnpottsmaxpen - (C - C_ab)*theta 0];
                    assert(opt.pnpottsmaxpen - (C - C_ab)*theta > 0);
                    Edges = [Edges ; Edgesonegroup];
                    %disp('Edge set!');
                end
                if size(Edges,1)>=1
                    pairwise_i = Edges(:,1);
                    pairwise_j = Edges(:,2);
                    pairwise_v = Edges(:,5);
                    W_pairwise = sparse(pairwise_i,pairwise_j,pairwise_v,alphabetasize+numGroups*2,alphabetasize+numGroups*2);
                    W_pairwise = (W_pairwise + tril(W_pairwise)' );
                    unaryterm(1,Edges(:,1)) = unaryterm(1,Edges(:,1)) - Edges(:,5)'*opt.pnpottsweight;
                    unaryterm(1,Edges(:,2)) = unaryterm(1,Edges(:,2)) + Edges(:,5)'*opt.pnpottsweight;
                end
            end
            
            BK_SetUnary(h,unaryterm*1000000);
            if size(W_pairwise,1) >=1
                BK_SetNeighbors(h,W_pairwise*opt.pnpottsweight*1000000);
            end
            maxflow = BK_Minimize(h);
            binaryLabeling = BK_GetLabeling(h); % 1 for alpha and 2 for beta
            binaryLabeling = binaryLabeling(1:alphabetasize);
            BK_Delete(h);
            if isempty(find(binaryLabeling==1)) || isempty(find(binaryLabeling==2))
                swap(alpha,beta) = 0;
                swap(beta,alpha) = 0;
                continue;
            end
            roiLabeling = binaryLabeling;
            assert(numel(unique(binaryLabeling))==2);
            roiLabeling(binaryLabeling == 1) = alpha;
            roiLabeling(binaryLabeling == 2) = beta;
            vLabeling(roi) = roiLabeling;
            newLabeling = vLabeling;
            % converged?
            diffcount = sum(Labeling(:) ~= newLabeling(:));
            if diffcount > 2 %&& ( current_e - new_e) > 1e-6
                Labeling = newLabeling;
                nmis = [nmis nmi(Labeling,labels)];
                new_e = kcutenergy(W,Labeling, opt, G, goodgroups);
                current_e = new_e;
                disp(['diffcount ' num2str(diffcount)]);
                kcutES = [kcutES current_e];
                disp(['SCUT energy: ' num2str(current_e)]);
                swap(:,[alpha, beta]) = 1;
                swap([alpha, beta],:) = 1;
                swap(alpha,alpha) = 0;
                swap(beta,beta) = 0;
                disp([num2str(alpha) ' and ' num2str(beta) ' swapped!']);
            else
                swap(alpha, beta) = 0;
                swap(beta, alpha) = 0;
                disp([num2str(alpha) ' and ' num2str(beta) ' unable to swap']);
            end
            %toc
        end
    end
    if itrNum > 100
        disp('Max iteration!');
        break;
    end
end


save('kcut.mat','Labeling','nmis');
disp(['nmi of kernel cut: ' num2str(nmi(Labeling,labels))]);
disp(['ncut energy: ' num2str(ncutEnergy( W, Labeling(:) ))]);
return
figure;plot(nmis);xlabel('Iteration');ylabel('NMI');
set(gca,'FontSize',16);
figure;plot(kcutES);xlabel('Iteration');ylabel('Normalized Cut + Group Prior');
set(gca,'FontSize',16);
