function [R,P_EQ,A]=hummer_szabo_clustering_A(K,INV_K, P_eq, A)

% Variable names in upper case correspond to the reduced system. Variable
% names in lower case correspond to the full system.
[n, ~] = size(K);

num_of_clusters=size(A,2);

P_EQ=zeros(size(A,2),1);
for i=1:size(A,1)
    for j=1:size(A,2)
        if A(i,j)==1
            P_EQ(j,1)=P_EQ(j,1)+P_eq(i);
        end
    end
end
    

ONE_VEC=ones(1,num_of_clusters);
one_vec=ones(1,n);

% Diagonal matrices of probabilities
D_N=diag(P_EQ);
D_n=diag(P_eq);

% A = aggregation matrix
% Creates a three-column matrix of length n: one column for each state.
% Shows how original states are aggregated.

% Calculating reduced matrix for the reduced two-state model
% Hummer-Szabo 2014, equation 12
% keyboard
R=P_EQ*ONE_VEC - D_N*(inv(A'*INV_K*D_n*A));

end