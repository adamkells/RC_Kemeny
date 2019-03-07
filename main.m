% Code by Adam Kells to use Kemeny constant for reaction coordinate
% identification
clear all
close all
% Create graph
% choose a number of edges
% Draw a random pair of states and connect them (if they are not already
% connected) until theres sufficient edges

M = 5; % number of edges within a cluster

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

G = digraph(A); % directed graph from adjacency matrix generated

tmp=reshape(K',[N^2,1]);
tmp2=tmp(tmp>0);
G.Edges.Weight = tmp2;
G.Edges.LWidths = 3*G.Edges.Weight/max(G.Edges.Weight);

% make a figure of the graph we have generated
h=plot(G);
h.LineWidth = G.Edges.LWidths;
% layout(h,'force3')
% view(3)
% keyboard

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
end_points = randi(N,[1,2]);
kem_max=0;
for i=1:N-1
    for j=i+1:N
        end_points=[i,j];
        
        % now for each other state I want to compute the commitor probability to
        % reach one state or the other first, details on what i'm doing are here:
        % www.emma-project.org/v2.2.1/api/generated/msmtools.analysis.committor.html
        [committor]=compute_commit(K,end_points);
        
        % the next bit is to now coarse grain the states based on
        % similarity of committor probability
        % do hummer-szabo clustering of rate matrix K
        [R,P_EQ,Aclus]=hummer_szabo_clustering(K', eq, committor);
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
        
        %save for eigenvector end points
        if i==min(b1,b2)
            if j==max(b1,b2)
                sec_eig_split = Aclus;
            end
        end
    end
end

% colour based on 'best' splitting
figure()
h = plot(G);
h.LineWidth = G.Edges.LWidths;
highlight(h,find(best_split(:,1)),'NodeColor','g')
title('Best splitting')

% colour based on second eigenvector sign
figure()
h = plot(G);
h.LineWidth = G.Edges.LWidths;
highlight(h,find(sign(K_eig_R(:,2))+1),'NodeColor','g')
title('Second eigenvector sign splitting')

% colour based on second eigenvector end points
figure()
h = plot(G);
h.LineWidth = G.Edges.LWidths;
highlight(h,find(sec_eig_split(:,1)),'NodeColor','g')
title('Second eigenvector end point splitting')