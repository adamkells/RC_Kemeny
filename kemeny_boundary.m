function [kemenyR]=kemeny_boundary(K,INV_K,eq,boundary,tmp2)

N=size(K,1);
A=zeros(N,3);
A(tmp2(1:boundary(1)),1)=1;
A(tmp2(boundary(1)+1:boundary(2)),2)=1;
A(tmp2(boundary(2)+1:end),3)=1;
[R,P_EQ,Aclus]=hummer_szabo_clustering_A(K',INV_K, eq, A);

% analyse eigenvalues vectors of clustered matrix R
[Reigs,~,rel__R,R_eig_R,R_eig_L]=spec_decomp(R);

% compute reduced kemeny (for 2 state clustering this is just the
% same as the slowest timescale
kemenyR=sum(-1./Reigs(2:end));
kemenyR=sum(-1./Reigs(2));

end