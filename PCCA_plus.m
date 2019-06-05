function [chi]=PCCA_plus(P,n)

%n=3; % number of clusters
m = size(P,1); % number of states

%% Different possibilities to define the initial distribution
%sd=rand(1,m); fprintf('\n initial distribution = random');
%sd=ones(1,m); fprintf('\n initial distribution = uniformly');
[sd, onv]=eigs(P',1,'lr'); fprintf('\n initial distribution = stationary');


%% Schur decomposition with correct initial distribution (sd)
Pd=diag(sqrt(sd))*P*diag(1./sqrt(sd));
[Q, R]=schur(Pd);
[Q, R]=SRSchur(Q, R, 1, n);
X=diag(1./sqrt(sd))*Q; X=X/X(1,1);
%% Clustering
chi=optimizeMetastab(X,n);

end

