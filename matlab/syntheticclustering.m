%Example of normalized cut (NC) or average association (AA) clustering using our kernel cut 
% The objective is NC or AA ONLY. No MRF term here so just use unary bounds.
% generate ring data, kernel and initial clustering
genringkernelandinit;
itrnum = 0;
energies = [];
energy_type = 'AA'; % 'NC' or 'AA'
while 1
    % current energy
    if strcmp(energy_type, 'NC')
        energies = [energies NormalizedCutEnergy( A, clustering )];
    elseif strcmp(energy_type, 'AA')
        % note that we are minimizing negative AA (or maximinzing AA)
        energies = [energies -1*AverageAssociationEnergy( A, clustering )];
    end
    disp(['current energy: ' num2str(energies(end))]);
    % linear bound
    unaries = KernelBound(A, 2, clustering, energy_type);
    % take the label correspond to the minimum unary term for each point
    [~, new_clustering] = min(unaries, [], 2); 
    % check convergence
    if 0 == sum(new_clustering~=clustering)
        disp('converged!');
        break;
    else
        clustering = new_clustering;
    end
    itrnum = itrnum + 1;
    if itrnum > 50
        disp('max iteration');
        break;
    end
end
figure;
plot(X(clustering==2,1),X(clustering==2,2),['r','s'], 'MarkerFaceColor','r','MarkerSize',5); hold on; 
plot(X(clustering==1,1),X(clustering==1,2),['b','s'], 'MarkerFaceColor','b','MarkerSize',5);axis equal
title('converged clustering');
figure,plot(energies,'--ks');title('Energies Convergence');xlabel('Iteration');ylabel('Energy');


