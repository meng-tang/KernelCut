Files = dir(['images' '/*.bmp']);
LengthFiles = length(Files);
mkdir('opticalflow');
% convert all bmp images to ppm images
for i=1:LengthFiles
    imgname=Files(i).name;
    imgname = imgname(1:(size(imgname,2)-4));
    img = imread(['images/' imgname '.bmp']);
    imwrite(img, ['opticalflow/' imgname '.ppm']);
end
Files = dir(['opticalflow' '/*.ppm']);
LengthFiles = length(Files);
for i=1:(LengthFiles-1)
    imgname=Files(i).name;
    imgname = imgname(1:(size(imgname,2)-4));
    disp(imgname);
    command = ['../../libs/LDOF/ldof ' ' opticalflow/' ...
         Files(i).name ' opticalflow/' Files(i+1).name];
    [status,cmdout] = system(command)
    disp(cmdout)
end
showldof;
