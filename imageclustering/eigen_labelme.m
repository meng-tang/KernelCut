load('kernel.mat');
[wx, wy] = size(W);
%continue;
x = 1 : wx;
S = (sum(W, 1));
D = full(sparse(x, x, S, wx, wy));
clear S x;
opts.issym=1;
opts.isreal = 1;
opts.disp=0;
opts.tol = 1e-9;
tic
[EigVect, EVal] = eigs(D - W, D, 8, 'sm',opts);
toc
clear D opts;

EigVal = diag(EVal);
clear Eval;

EigVal(1:end) = EigVal(end:-1:1);
EigVect(:, 1:end) = EigVect(:, end:-1:1);

% normalzied the eigenvectors
for i=1:8
    EigVect(:,i) = EigVect(:,i) / norm(EigVect(:,i));
end
save(['eigen.mat'],'EigVect','EigVal');