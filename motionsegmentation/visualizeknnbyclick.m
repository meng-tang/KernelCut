db = 'ducks01';
imgname = 'ducks01_0300';
global rgbimg
rgbimg = imread([db '/images/' imgname '.bmp']);
global H
global W
global C
[H, W, C] = size(rgbimg);
iptsetpref('ImshowBorder','tight');
set(0,'DefaultFigureMenu','none');
global f
f= figure;imshow(rgbimg);
fileID = fopen( [db '/knn/' imgname '.bin'],'r');
global knnidx
knnidx = fread(fileID,[H*W, 50],'int32');
fclose(fileID);
set(f,'WindowButtonDownFcn',@mycallback);
