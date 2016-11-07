function [  ] = mycallback( hObject,~)
%MYCALLBACK Summary of this function goes here
%   Detailed explanation goes here
pos = get(hObject,'CurrentPoint');
global H
global W
global rgbimg
global knnidx
x = pos(1);
y = H+1 - pos(2);
disp(['You clicked X:',num2str(x),', Y:',num2str(y)]);
imshow(rgbimg);
hold on
p_idx = (x-1)*H+y-1+1;
q_idx = knnidx(p_idx,:);
q_x = (int32(q_idx-1) - mod(int32(q_idx-1),H) ) / int32(H) + 1;
q_y = mod(int32(q_idx-1),H) + 1;
for i=1:50
    line([x q_x(i)],[y q_y(i)],'Color','red');
end
end
