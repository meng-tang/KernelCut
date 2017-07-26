close all
clear
Files = dir('opticalflow/*.flo');
LengthFiles = length(Files);
addpath('../../libs/flow-code-matlab');
for i=1:LengthFiles
    imgname=Files(i).name;
    imgname = imgname(1:(size(imgname,2)-4));
    ldof = readFlowFile(['opticalflow/' imgname '.flo']);
    ldof_x = ldof(:,:,1);
    ldof_y = ldof(:,:,2);
    ldof_m = ldof_x.*ldof_x + ldof_y.*ldof_y;
    ldof_m = sqrt(ldof_m);
    %gray = mat2gray(ldof_m);
    gray = (ldof_m)*0.2508;
    X = gray2ind(gray, 128);
    rgb = ind2rgb(X, jet(128));
    disp(['min max flow: ' num2str(min(ldof_m(:))) ' ' num2str(max(ldof(:)))]);
    imwrite(rgb, ['opticalflow/' imgname '.jpg']);
end