function [Keigs,eq,rel_exact,K_eig_R,K_eig_L]=spec_decomp(K)

%compute all the left and right eigenvectors, eigenvalues, relaxation times
%and equilibrium probability distribution

N=size(K,1);
%calculate equilibrium from spectral decomposition
[eigvec,eigval]=eig(K); % diagonalize K, eigvec stores the eigenvectors, eigval the eigenvalues
[dsorted,index]=sort(real(diag(eigval)),'descend'); % sort the eigenvalues.
Keigs=dsorted;
eigvec=eigvec(:,index);
eq=eigvec(:,1)/sum(eigvec(:,1));
rel_exact=-1/dsorted(2);

for i=2:N
    eigvec(:,i)=eigvec(:,i)./sqrt(sum(eigvec(:,i).^2./eq));
end
eigvec(:,1)=eigvec(:,1)./sum(eigvec(:,1));
K_eig_R=eigvec;

%calculate equilibrium from spectral decomposition
[eigvec,eigval]=eig(K'); % diagonalize K, eigvec stores the eigenvectors, eigval the eigenvalues
[~,index]=sort(diag(eigval),'descend'); % sort the eigenvalues.
eigvec=eigvec(:,index);

for i=2:N
    eigvec(:,i)=eigvec(:,i)./sqrt(sum(eigvec(:,i).^2.*eq));
end
eigvec(:,1)=eigvec(:,1)./sum(eigvec(:,1).*eq);
K_eig_L=eigvec;

for i=1:N
    for j=1:N
        if abs(K_eig_L(i,j)-K_eig_R(i,j)/eq(i))>0.01
            K_eig_L(:,j)=K_eig_L(:,j)*(-1);
        end
    end
end

end