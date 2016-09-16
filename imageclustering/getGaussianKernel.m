function K = getGaussianKernel(X, sigma)
    
%% F: feature matrix: nImg * nFeature
 

%     sigma = 0.5;   
 
    [nImg] = size(X, 1);
    aa = sum(X.*X,2);  
    D = repmat(aa, 1, nImg) + repmat(aa', nImg, 1) - 2*X*X';    
    D = D - diag(diag(D));
    D = abs(D);
    
    K = exp(-D/(2*sigma^2));
    
%     X =X';
%     n1sq = sum(X.^2,1);
%     nImg = size(X,2);
%     D = (ones(nImg,1)*n1sq)' + ones(nImg,1)*n1sq -2*(X'*X);
%     K = exp(-D/(2*sigma^2));
        
    
    % debug
%     imagesc(K); colorbar
end
