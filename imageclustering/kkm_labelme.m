load('kernel.mat');
load('labels.mat');
N = size(W,1);
initmode = 'sc'; % rand or sc
if strcmp(initmode,'rand')
    Labeling = randi([1,8],N,1); % random init
elseif strcmp(initmode,'sc')
    load('sc.mat');
end
initLabelingNcutE = ncutEnergy( W, Labeling(:) );
itrNum = 0;
ncutES = initLabelingNcutE;
nmis = nmi(Labeling,labels);
wp = sum(W,2);
while 1
    [ dist_to_mu ] = takebound( W, wp,Labeling(:) );
    [mintestue, newLabeling] = min(dist_to_mu,[],2);
    % converged?
    diffcount = sum(Labeling(:) ~= newLabeling(:));
    if diffcount > 1 % 1 for deep kernel, 5 for shallow kernel
        Labeling = newLabeling;
        disp(['diffcount ' num2str(diffcount)]);
        ncutE = ncutEnergy( W, Labeling(:) );
        ncutES = [ncutES ncutE];
        nmis = [nmis nmi(Labeling,labels)];
        disp(['NCUT energy: ' num2str(ncutE)]);
        if numel(ncutES)>=2 && ( ncutES(end-1)-ncutES(end) ) < 1e-6
            disp('Energy converged!');
            break;
        end
    else
        disp('converged');
        break;
    end
    itrNum = itrNum + 1;
    if itrNum>100
        break;
    end
end
save(['kkm_' initmode '.mat'],'Labeling','ncutES','nmis');


disp(['nmi of kernel k-means: ' num2str(nmi(Labeling,labels))]);
disp(['ncut energy: ' num2str(ncutEnergy( W, Labeling(:) ))]);

figure;plot(ncutES);xlabel('Iteration');ylabel('negative normlized association');
set(gca,'FontSize',16);
figure;plot(nmis);xlabel('Iteration');ylabel('NMI');
set(gca,'FontSize',16);
