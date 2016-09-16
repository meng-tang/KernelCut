% load eigensolution
load eigen.mat
load('kernel.mat');
% kmeans
num_vec = 8;
num_cluster = 8;
[ embedding, weights ] = getembeddingandweights( EigVect, EigVal, 'eigenvector', num_vec, full(sum(W,2)) );
[Labeling,~,~] = kmeans(embedding,num_cluster,...
            'MaxIter',10000,'Start','plus');
Labeling = int32(Labeling);
load('labels.mat');
labels = int32(labels);
save('sc.mat','Labeling');% spectral clustering
v = nmi(Labeling,labels);
disp(['nmi of spectral clustering: ' num2str(v)]);
disp(['ncut energy: ' num2str(ncutEnergy( W, Labeling(:) ))]);