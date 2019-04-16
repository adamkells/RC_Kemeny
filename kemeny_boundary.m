function [kemenyR]=kemeny_boundary(K,INV_K,eq,A,jj,param)

% computing hummer-szabo reduced matrix
if jj==0
    [R,~,~]=hummer_szabo_clustering_A(K,INV_K, eq, A);
elseif jj==1
    R=localeq(K',eq,A);
end

% analyse eigenvalues vectors of clustered matrix R
[Reigs,~,~,~,~]=spec_decomp(R);

% compute variational parameter
if param==0
    kemenyR=sum(-1./Reigs(2:end));
elseif param==1
   kemenyR=sum(-1./Reigs(2)); 
end
%kemenyR=sum(-1./Reigs(2:(end-1)));
end