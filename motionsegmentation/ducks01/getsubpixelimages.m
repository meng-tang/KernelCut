% add gaussian noise to RGB image
% convert to LAB image
% scale each channel to have unit STD
Files = dir('images/*.bmp');
LengthFiles = length(Files);
mkdir('subpixelimages');
for i=1:LengthFiles
    imgname=Files(i).name;
    rgbimg = imread(['images/' imgname]);
    imgname = imgname(1:(size(imgname,2)-4));
    if ~strcmp(imgname,'book')
        ;%continue;
    end
    [M N C] = size(rgbimg);
    rgbimg = double(rgbimg)/255.0;
    %% add gaussian noise to rgbimg (std = sqrt(0.0002)*255)
    rgbimg = imnoise(rgbimg,'gaussian',0,0.0002);
    %rgbimg = rgbimg+rand(size(rgbimg))*1/255.0;
    %rgbimg = rgbimg / (1.0+1.0/255.0);
    % transform to LAB space
    labimg = rgb2lab(rgbimg);
    lab_l = labimg(:,:,1);lab_a = labimg(:,:,2);lab_b = labimg(:,:,3);
    lab_l = lab_l(:);lab_a = lab_a(:);lab_b = lab_b(:);
    %figure;scatter3(lab_l(1:50:end),lab_a(1:50:end),lab_b(1:50:end),'x');title('LAB');
    %return
    %% normalize lab so that each channel has unit std
    X = [reshape(labimg(:,:,1),M*N,1) reshape(labimg(:,:,2),M*N,1) ...
        reshape(labimg(:,:,3),M*N,1) ];
    x_std = std(X,0,1);
    labimg(:,:,1) = labimg(:,:,1) / x_std(1);
    labimg(:,:,2) = labimg(:,:,2) / x_std(2);
    labimg(:,:,3) = labimg(:,:,3) / x_std(3);
    fileID = fopen( ['subpixelimages/' imgname '.bin'],'w');
    fwrite(fileID,labimg,'double');
    fclose(fileID);
    disp([num2str(i) ' th image: ' imgname]);
end