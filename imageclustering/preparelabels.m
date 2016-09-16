load Dlabelme.mat
N = size(Dlabelme,2);
labels = zeros(N,1,'int32');
for n=1:N
    filename = Dlabelme(n).annotation.filename;
    scenetype = filename(1:(find(filename=='_')-1));
    %disp(scenetype);
    if strcmp(scenetype,'coast')
        sceneid = 1;
    elseif strcmp(scenetype,'mountain')
        sceneid = 2;
    elseif strcmp(scenetype,'forest')
        sceneid = 3;
    elseif strcmp(scenetype,'opencountry')
        sceneid = 4;
    elseif strcmp(scenetype,'street')
        sceneid = 5;
    elseif strcmp(scenetype,'insidecity')
        sceneid = 6;
    elseif strcmp(scenetype,'tallbuilding')
        sceneid = 7;
    elseif strcmp(scenetype,'highway')
        sceneid = 8;
    else
        disp('Unrecognized category!');
        return;
    end
    labels(n) = sceneid;
end
save('labels.mat','labels');
