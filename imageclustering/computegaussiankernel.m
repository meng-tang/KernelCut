if ~exist('featureOB.mat', 'file')
    param.dataset = 'images/';
    load Dlabelme.mat
    N = size(Dlabelme,2);
    for n=1:N
        filename = Dlabelme(n).annotation.filename;
        param.imageList{n} = filename;
    end
    addpath( 'object-bank/code/partless/')
    addpath(genpath('object-bank/code/partless/code/'));
    addpath(genpath('object-bank/code/partless/code/LSVM'));
    imageFolder = param.dataset;
    imageList = param.imageList;


    % Pre-allocate gist:
    nFeature = 44604 ;
    nImg =  length(imageList);

    featureOB= zeros([nImg nFeature]); 


    for i=1:nImg
       i

        imgPath = [imageFolder, imageList{i} ];
        featureOB(i, :) =   getFeatureOB(imgPath);

    end  


    save featureOB.mat  featureOB

else
   load featureOB.mat
end


paramOB.sigma = 22;
W = getGaussianKernel(featureOB, paramOB.sigma);
save('kernel.mat','W');
