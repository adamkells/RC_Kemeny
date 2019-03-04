% Code by Adam Kells to use Kemeny constant for reaction coordinate
% identification
clear all
% Create graph
% choose a number of edges
% Draw a random pair of states and connect them (if they are not already
% connected) until theres sufficient edges

M = 60; % number of edges

[A]=make_graph(M);
%[A]=make_lin_graph(M);

%if i want a two state-esque network
multi_state = 5;
ncon = 5;
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



N=size(A,1);

% make a figure of the graph we have generated
G = digraph(A); % directed graph from adjacency matrix generated
h=plot(G)
keyboard
% create a matrix which describes the rate of transition between states
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

% I want to choose two states at random and define a reaction coordinate
% based on the commitor probability to reach one state over the other and
% then search for the two states which maximise the kemeny constant
% could maybe generalise to three stable states? how would that work?

%%
% lets choose two states
end_points = randi(N,[1,2]);
kem_max=0;
for i=1:N-1
    for j=i+1:N
        %end_points=[2,55];
        end_points=[i,j];
        %end_points = [1,N];
        
        % now for each other state I want to compute the commitor probability to
        % reach one state or the other first
        % details on what i'm doing are here:
        % www.emma-project.org/v2.2.1/api/generated/msmtools.analysis.committor.html
        [committor]=compute_commit(K,end_points);
        
        % the next bit is to now coarse grain the states based on
        % similarity of committor probability
        
        % now we just need to decide upon a number of clusters
        
        [R,P_EQ,Aclus]=hummer_szabo_clustering(K', eq, committor);
        
        [Reigs,~,rel__R,R_eig_R,R_eig_L]=spec_decomp(R);
        kemenyR=sum(-1./Reigs(2:end));
        if kemenyR>kem_max
            kem_max = kemenyR;
            best_wells = end_points;
            best_split = Aclus;
            h = plot(G)
            highlight(h,find(Aclus(:,1)),'NodeColor','g')
        end
    end
end