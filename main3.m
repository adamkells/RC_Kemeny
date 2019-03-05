% Code by Adam Kells to use Kemeny constant for reaction coordinate
% identification
clear all
close all
% Create graph
% choose a number of edges
% Draw a random pair of states and connect them (if they are not already
% connected) until theres sufficient edges

M = 25; % number of edges within a cluster

[A]=make_graph(M); % make random graph
%[A]=make_lin_graph(M); % make linear chain graph

%if i want a multi state-esque network
multi_state = 2; % number of clusters
ncon = 1; % number of intercluster links
for i = 1:multi_state-1
    [A1]=make_graph(M);
    A2 = blkdiag(A,A1);
    for j=1:ncon
        rand_a = randi(size(A,1),1);
        rand_b = randi([size(A,1)+1,size(A2,1)],1);
        if A2(rand_a,rand_b)==0
            A2(rand_a,rand_b)=1;
            A2(rand_b,rand_a)=1;
        end
    end
    A = A2;
end
N=size(A,1); % total number of nodes in the network

% make a figure of the graph we have generated
G = digraph(A); % directed graph from adjacency matrix generated
h=plot(G);

% create a matrix which describes the rate of transition between states
% evenly divide rate among number of states linked to
for i = 1:N
    for j = 1:N
        K(i,j)=A(i,j)/sum(A(:,j));
    end
end

for i = 1:N
    K(i,i)=-sum(K(i,:).*A(i,:));
end

% do spectral decomposition
[Keigs,eq,rel_exact,K_eig_R,K_eig_L]=spec_decomp(K');
kemeny = sum(-1./Keigs(2:end));

% what are the endpoints suggested by the second eigenvector
[a,b1] = min(K_eig_R(:,2));
[a,b2] = max(K_eig_R(:,2));

% I want to choose two states at random and define a reaction coordinate
% based on the commitor probability to reach one state over the other and
% then search for the two states which maximise the kemeny constant
% could maybe generalise to three stable states? i would need to consider
% committor between each pair of states

%%
% lets choose two states
tic
end_points = randi(N,[1,2]);
kem_max=0;
num_clusters = 3;
for i=1:N-1
    for j=i+1:N
        end_points=[i,j];
        [i,j]
        
        % now for each other state I want to compute the commitor probability to
        % reach one state or the other first, details on what i'm doing are here:
        % www.emma-project.org/v2.2.1/api/generated/msmtools.analysis.committor.html
        [committor]=compute_commit(K,end_points);
        
        % the next bit is to now coarse grain the states based on
        % similarity of committor probability
        % do hummer-szabo clustering of rate matrix K
        kem_local = 0;
        [~,tmp2]=sort(committor);
        for i1=1:N-2
            for j1=(i1+1):N-1
                A=zeros(N,num_clusters);
                A(tmp2(1:i1),1)=1;
                A(tmp2(i1+1:j1),2)=1;
                A(tmp2(j1+1:end),3)=1;
                [R,P_EQ,Aclus]=hummer_szabo_clustering_A(K', eq, A);
                
                % analyse eigenvalues vectors of clustered matrix R
                [Reigs,~,rel__R,R_eig_R,R_eig_L]=spec_decomp(R);
                
                % compute reduced kemeny (for 2 state clustering this is just the
                % same as the slowest timescale
                kemenyR=sum(-1./Reigs(2:end));
                
                %  save the info for the iteration which maximises kemeny (this may
                %  not be unique), will only save the first choice found
                if kemenyR>kem_max
                    kem_max = kemenyR;
                    best_wells = end_points;
                    best_split = Aclus;
                end
                if kemenyR>kem_local
                    kem_local = kemenyR;
                    local_split(i,j,:,:)=Aclus;
                end
                
            end
        end
        %save for eigenvector end points
        if i==min(b1,b2)
            if j==max(b1,b2)
                sec_eig_split = best_split;
            end
        end
    end
end
kemeny_time=toc;

% Here I want to do the exhaustive version where I try every combination of
% states to cluster in to.
tic

kem_max2=0;
count=0;
for i=1:(N-2)
    clus1=nchoosek([1:N],i);
    [i]
    for k=1:size(clus1,1)
        clus1_tmp=clus1(k,:);
        for j=1:N-i-1
            clus2=nchoosek(setdiff([1:N],clus1_tmp),j);
            for kk=1:size(clus2,1)
                count=count+1;
                A=zeros(N,3);
                clus2_tmp=clus2(kk,:);
                clus3_tmp=setdiff([1:N],[clus1_tmp,clus2_tmp]);
                A(clus1_tmp,1)=1;
                A(clus2_tmp,2)=1;
                A(clus3_tmp,3)=1;
                [R,P_EQ,Aclus]=hummer_szabo_clustering_A(K', eq, A);
                [Reigs,~,rel__R,R_eig_R,R_eig_L]=spec_decomp(R);
                kemenyR=sum(-1./Reigs(2:end));
                if kemenyR>kem_max2
                    kem_max2 = kemenyR;
                    best_split_exhaus = Aclus;
                end
            end
        end
    end
end
exhaust_time=toc;


% colour based on 'best' splitting
figure()
h = plot(G);
highlight(h,find(best_split(:,1)),'NodeColor','g')
highlight(h,find(best_split(:,2)),'NodeColor','r')
title('Best splitting')

% colour based on exhaustive
figure()
h = plot(G);
highlight(h,find(best_split_exhaus(:,1)),'NodeColor','g')
highlight(h,find(best_split_exhaus(:,2)),'NodeColor','r')
title('Exhaustive splitting')

% colour based on second eigenvector end points
figure()
h = plot(G);
highlight(h,find(sec_eig_split(:,1)),'NodeColor','g')
highlight(h,find(sec_eig_split(:,2)),'NodeColor','r')
title('Second eigenvector end point splitting')