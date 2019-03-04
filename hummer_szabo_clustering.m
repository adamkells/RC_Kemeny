function [R,P_EQ,A]=hummer_szabo_clustering(K, P_eq, committor)

% Variable names in upper case correspond to the reduced system. Variable
% names in lower case correspond to the full system.
N=size(K,1);
n_cluster = 2;
A=zeros(N,n_cluster);
for i = 1:N
    if committor(i)<0.5
        A(i,1) = 1;
    elseif committor(i)>0.5
        A(i,2) = 1;
    else
        A(i,rand(n_cluster))=1;
    end
end

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
R=P_EQ*ONE_VEC - D_N*(inv(transpose(A)*(inv(P_eq*one_vec-K))*D_n*A));

end