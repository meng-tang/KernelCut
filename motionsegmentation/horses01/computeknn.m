function [] = computeknn(whichimage)% find KNN (K = 400)
if nargin == 0
    whichimage = 'all';
end
Files = dir('images/*.bmp');
LengthFiles = length(Files);
addpath('../../libs/flow-code-matlab');
% relative scale of xy:
positionscale = 0;
% relative scale of motion:
motionscale = 0.1;
targetdir = 'knn';
if 0==exist(targetdir,'dir')
    mkdir(targetdir);
end
for i=1:(LengthFiles-1)
    imgname=Files(i).name;
    imgname = imgname(1:(size(imgname,2)-4));
    if ~strcmp(imgname,whichimage)
        ;%continue;
    end
    % load rgb image
    rgbimg = imread(['images/' imgname '.bmp']);
    [M N C] = size(rgbimg);
    % read subpixel image
    fileID = fopen( ['subpixelimages/' imgname '.bin'],'r');
    labimg = fread(fileID,numel(rgbimg),'double');
    fclose(fileID);
    labimg = reshape(labimg,M,N,C);
    xpos = zeros(M,N,'double');
    for n=1:N
        xpos(:,n) = n;
    end
    ypos = zeros(M,N,'double');
    for m=1:M
        ypos(m,:) = m;
    end
    pos = [reshape(xpos,M*N,1) reshape(ypos,M*N,1)]*positionscale;
    % run knn
    K = 400;
    ldof = readFlowFile(['opticalflow/' imgname 'LDOF.flo']);
    X = [reshape(labimg(:,:,1),M*N,1) reshape(labimg(:,:,2),M*N,1) ...
            reshape(labimg(:,:,3),M*N,1)];
    motion = [reshape(ldof(:,:,1),M*N,1) reshape(ldof(:,:,2),M*N,1)];
    %figure,imagesc(ldof(:,:,1));colorbar
    %figure,imagesc(ldof(:,:,2));colorbar
    motion(:,1) = motion(:,1) *motionscale;
    motion(:,2) = motion(:,2) *motionscale;
    if motionscale > 1e-10
        X = [X pos motion];
    else
        X = [X pos];
    end
    [knnidx, knnd] = knnsearch(X,X,'K',K);
    % save distance to the k th neighbor
    fileID = fopen( [targetdir '/' imgname '.bin'],'w');
    fwrite(fileID,knnidx(:,1:8:end),'int32');
    fclose(fileID);
    disp(['saved into bin files ' targetdir '/' imgname]);
end
end
